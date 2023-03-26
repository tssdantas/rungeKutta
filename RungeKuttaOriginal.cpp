#include <iostream>
#include <cmath>

//fonte: https://www.geeksforgeeks.org/runge-kutta-4th-order-method-solve-differential-equation/
 
using namespace std;

// Function to define the differential equation dy/dx
double f(double x, double y)
{
    return x*y; // Example function: y' = xy
}

// Function to implement the Runge-Kutta 4th Order Method
double rungeKutta(double x0, double y0, double x, double h)
{
    // Count number of iterations using step size or
    // step height h
    int n = (int)((x - x0) / h);

    double k1, k2, k3, k4;

    // Iterate for number of iterations
    double y = y0;
    for (int i = 1; i <= n; i++) {
        // Apply Runge-Kutta Formulas to find
        // next value of y
        k1 = h * f(x0, y);
        k2 = h * f(x0 + 0.5 * h, y + 0.5 * k1);
        k3 = h * f(x0 + 0.5 * h, y + 0.5 * k2);
        k4 = h * f(x0 + h, y + k3);

        // Update next value of y
        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);

        // Update next value of x
        x0 = x0 + h;
    }

    return y;
}

// Main method
int main()
{
    double x0 = 0, y = 1, x = 1, h = 0.2;
    cout << "The value of y at x is : " << rungeKutta(x0, y, x, h) << endl;

    return 0;
}