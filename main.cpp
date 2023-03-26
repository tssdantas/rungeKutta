#include <iostream>
#include <cmath>
#include <vector>

//fonte: https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
 
using namespace std;

// Function to define the differential equation dy/dx
vector<double> f(double t0, double h, vector<double> xDot, int k, vector<double> _modelInputs, vector<double> lastK = {0, 0, 0})
{
    if (xDot.size() != lastK.size())
    {
        throw ("Erro de dimensionamento dos vetores");
    }

    vector<double> newK;
    newK.resize(xDot.size());
    
    double xDot0, xDot1, xDot2; 

    switch (k)
    {     
        case 1:
            xDot0 = xDot.at(0);
            xDot1 = xDot.at(1);
            xDot2 = xDot.at(2);
            break;
        
        case 2:
            xDot0 = xDot.at(0) + h/2;
            xDot1 = xDot.at(1) + lastK.at(1)/2;
            xDot2 = xDot.at(2) + lastK.at(2)/2;
            break;

        case 3:
            xDot0 = xDot.at(0) + h/2;
            xDot1 = xDot.at(1) + lastK.at(1)/2;
            xDot2 = xDot.at(2) + lastK.at(2)/2;
            break;

        case 4:
            xDot0 = xDot.at(0) + h;
            xDot1 = xDot.at(1) + lastK.at(1);
            xDot2 = xDot.at(2) + lastK.at(2);
            break;
        
        default:
            throw ("Erro de chamada a funcao");
            break;
    }


    vector<double> yDot(3,0);
    yDot.at(0) = xDot0;
    yDot.at(1) = xDot2;
    yDot.at(2) = -2 * xDot2 - xDot1 + sin(xDot0);

    //multipica tudo por h e retorna para o u
    for (unsigned int i = 0; i < yDot.size(); i++ )
    {
        yDot[i] = h*yDot[i];
    }
    return yDot;
    //return t*y; // Example function: y' = xy
}

// Function to define the differential equation dy/dx
vector<double> g(double t0, double h, vector<double> xDot, int k, vector<double> _modelInputs, vector<double> lastK = {0, 0, 0})
{

    if (_modelInputs.size() != 4)
    {
        throw ("Erro de dimensionamento dos vetores: _modelInputs");
    }

    // Inputs
    // 0 - QH
    // 1 - Ca0
    // 2 - T0  (temp reator)
    // 4 - u (F/V) taxa de diluicao
    double QH = _modelInputs[0];
    double CA0 = _modelInputs[1];
    double T0 = _modelInputs[2];
    double u =  _modelInputs[3];

    if (xDot.size() != lastK.size())
    {
        throw ("Erro de dimensionamento dos vetores");
    }

    // state variables (xDot)
    // xDot[0] - Ca
    // xDot[1] - Cb
    // xDot[2] - Temp

    vector<double> newK;
    newK.resize(xDot.size());
    
    double xDot0, xDot1, xDot2; 

    switch (k)
    {     
        case 1:
            xDot0 = xDot.at(0);
            xDot1 = xDot.at(1);
            xDot2 = xDot.at(2);
            break;
        
        case 2:
            xDot0 = xDot.at(0) + h/2;
            xDot1 = xDot.at(1) + lastK.at(1)/2;
            xDot2 = xDot.at(2) + lastK.at(2)/2;
            break;

        case 3:
            xDot0 = xDot.at(0) + h/2;
            xDot1 = xDot.at(1) + lastK.at(1)/2;
            xDot2 = xDot.at(2) + lastK.at(2)/2;
            break;

        case 4:
            xDot0 = xDot.at(0) + h;
            xDot1 = xDot.at(1) + lastK.at(1);
            xDot2 = xDot.at(2) + lastK.at(2);
            break;
        
        default:
            throw ("Erro de chamada a funcao");
            break;
    }

    // PARAMETROS FIXOS
    long double k10 = 1.287*pow(10,12);
    long double k20 = 1.287*pow(10,12);
    long double k30 = 9.043*pow(10,9) ;

    double Er1 = 9758.3;
    double Er2 = 9758.3;
    double Er3 = 8560.0;

    double DeltaH1 = 4.2;
    double DeltaH2 = -11;
    double DeltaH3 = -41.85; 

    double rho =  0.9342 ;
    double cp  =  3.01   ;

    // constante cinetica !! nao sao iguais aos k1, k2, k3 do metodo RK !!!
    double k1 = k10*exp(-Er1/xDot[2]);
    double k2 = k20*exp(-Er2/xDot[2]);
    double k3 = k30*exp(-Er3/xDot[2]);  

    //modelo em estado estacionario
    vector<double> yDot(3,0);
    yDot.at(0) = -k1*xDot0  - k3*pow(xDot0, 2) + (CA0 - xDot0)*u;
    yDot.at(1) = k1*xDot0 - k2*xDot1 - xDot1*u;
    yDot.at(2) = (T0 - xDot2)*u + ( (-DeltaH1*k1*xDot0) + (-DeltaH2*k2*xDot1) + (-DeltaH3*k3*pow(xDot0, 2)) + QH )/(rho*cp);


    //multipica tudo por h e retorna para o u
    for (unsigned int i = 0; i < yDot.size(); i++ )
    {
        yDot[i] = h*yDot[i];
    }
    return yDot;
    //return t*y; // Example function: y' = xy
}

// Function to implement the Runge-Kutta 4th Order Method
vector<double> rungeKutta(double t0, double tf, double h, vector<double> y0, vector<double> _modelInputs)
{
    // Count number of iterations using step size or
    // step height h
    int n = (int) ((tf - t0) / h);

    vector<double> k1, k2, k3, k4;
    k1.resize(y0.size());
    k2.resize(y0.size());
    k3.resize(y0.size());
    k4.resize(y0.size());

    double t = t0;

    // Iterate for number of iterations
    vector<double> y = y0;
    for (int i = 0; i <= n; i++) {

        // Apply Runge-Kutta Formulas to find
        // next value of y
        k1 = g(t, h, y, 1, _modelInputs);
        k2 = g(t, h, y, 2, _modelInputs, k1);
        k3 = g(t, h, y, 3, _modelInputs, k2);
        k4 = g(t, h, y, 4, _modelInputs, k3);

        //original
        // k1 = h * f(x0, y);
        // k2 = h * f(x0 + 0.5 * h, y + 0.5 * k1);
        // k3 = h * f(x0 + 0.5 * h, y + 0.5 * k2);
        // k4 = h * f(x0 + h, y + k3);

        // Update next value of y
        cout << i <<": ";
        for (unsigned int ii = 0; ii < y.size(); ii++)
        {
            //cout << y[ii] << " | ";

            y[ii] += (1.0/6.0)*(k1[ii] + 2*k2[ii] + 2*k3[ii] + k4[ii]);

            cout << y[ii] << " | ";
        }
        cout << endl;

        //y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

        // Update next value of t
        t = t + h;
    }

    return y;
}


// Main method
int main()
{
    // Inputs for the reactor model
    double QH  = -451.509;
    double u   = 19.52;  // %F/V (taxa de diluicao)
    double CA0 = 5.0;
    double T0  = 403.15 ;
    vector<double> reactorInputs = { QH, CA0, T0, u };
    //vector<double> reactorInputs = { 0.0, 0.0, 0.0, 0.0 };

    // Inital state reator model
    vector<double> yDot0 = { 1.0, 0.9, 400 }; // CA0, CB0, T0
    //vector<double> yDot0 = { 1.268, 0.9664, 410.9213 }; // CA0, CB0, T0
    // time domain inputs
    double t0 = 0, tf = 5, h = 0.05;

    // double t0 = 0; double tf = 1;
    // vector<double> yDot0 = {0, 1, 0};
    // double h = 0.1;

    vector<double> state = rungeKutta(t0,  tf, h, yDot0, reactorInputs);
    cout << "The value of state vector at t is : ";
    for (auto &v:state) { std::cout << v << " "; };
    cout << endl; 

    return 0;
}