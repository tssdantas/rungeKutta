#include <iostream>
#include <cmath>

using namespace std;

// Define the differential equations
// y -> f(f(x,y,z))
double f(double x, double y, double z) { 
    return z;
}
// z -> f(g(x,y,z))
double g(double x, double y, double z) {
    return -2 * z - y + sin(x);
}

// Runge-Kutta 4th order method
void RungeKutta4(double x0, double y0, double z0, double h, int n) {
    double k1, k2, k3, k4;
    double l1, l2, l3, l4;
    double x = x0, y = y0, z = z0;

    // Print initial values
    cout << "x = " << x << ", y = " << y << ", z = " << z << endl;

    // Iteratively apply the RK4 method
    for (int i = 1; i <= n; i++) {
        k1 = h * f(x, y, z);
        l1 = h * g(x, y, z);
        k2 = h * f(x + h / 2, y + k1 / 2, z + l1 / 2);
        l2 = h * g(x + h / 2, y + k1 / 2, z + l1 / 2);
        k3 = h * f(x + h / 2, y + k2 / 2, z + l2 / 2);
        l3 = h * g(x + h / 2, y + k2 / 2, z + l2 / 2);
        k4 = h * f(x + h, y + k3, z + l3);
        l4 = h * g(x + h, y + k3, z + l3);

        y = y + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        z = z + (l1 + 2 * l2 + 2 * l3 + l4) / 6;
        x = x + h;

        // Print the updated values at each iteration
        cout << "x = " << x << ", y = " << y << ", z = " << z << endl;
    }
}

int main() {
    // Define the initial conditions
    double x0 = 0, y0 = 1, z0 = 0;

    // Define the step size and number of iterations
    double h = 0.1;
    int n = 10;

    // Solve the system of differential equations using RK4 method
    RungeKutta4(x0, y0, z0, h, n);

    return 0;
}