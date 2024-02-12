#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>

using namespace std;

const int tran = 20000, n = 1000, nite = 30000;
const double h = 0.05, k2 = 8.0, eta = 1.0;
const double pi = 4.0 * atan(1.0);

vector<double> omega(n), theta(n), dth(n), tho(n), dth_o(n);

/**
 * Generates a random number uniformly distributed between 0.0 and 1.0 using
 * Mersenne Twister pseudo-random number generator.
 *
 * @return The generated random number.
 */
double rand_uniform()
{
    static mt19937 generator;
    static uniform_real_distribution<double> distribution(0.0, 1.0);
    return distribution(generator);
}

/**
 * Calculates the derivatives for the given parameters and updates the dth vector.
 *
 * @param t the time parameter
 * @param theta the vector of angles
 * @param dth the vector to be updated with derivatives
 * @param k1 the constant
 * @param k2 the constant
 * @param omega the vector of angular velocities
 * @param r1 the radius
 * @param shi1 the angle
 * @param r2 the radius
 * @param shi2 the angle
 *
 * @return void
 *
 * @throws None
 */
void derivs(double t, const vector<double> &theta, vector<double> &dth, double k1, double k2, const vector<double> &omega, double r1, double shi1, double r2, double shi2)
{
    double alpha1 = 0.0, beta1 = 0.0;
    for (int i = 0; i < n; ++i)
    {
        double dth1 = k1 * pow(r1, alpha1 + 1) * sin(shi1 - theta[i]);
        double dth2 = k2 * r2 * pow(r1, beta1 + 1) * sin(shi2 - shi1 - theta[i]);
        dth[i] = omega[i] + dth1 + dth2;
    }
}

/**
 * Runge-Kutta method for solving ordinary differential equations.
 *
 * @param y vector of dependent variables
 * @param dydx vector of derivatives of dependent variables
 * @param n number of variables
 * @param x independent variable
 * @param h step size
 * @param yout output vector of dependent variables
 * @param k1 first constant parameter
 * @param k2 second constant parameter
 * @param omega vector of constant parameters
 * @param r1 constant parameter
 * @param shi1 constant parameter
 * @param r2 constant parameter
 * @param shi2 constant parameter
 */
void rk4(vector<double> &y, const vector<double> &dydx, int n, double x, double h, vector<double> &yout, double k1, double k2, const vector<double> &omega, double r1, double shi1, double r2, double shi2)
{
    vector<double> dym(n), dyt(n), yt(n);
    double hh = h * 0.5;
    double h6 = h / 6.0;
    double xh = x + hh;
    for (int i = 0; i < n; ++i)
    {
        yt[i] = y[i] + hh * dydx[i];
    }
    derivs(xh, yt, dyt, k1, k2, omega, r1, shi1, r2, shi2);
    for (int i = 0; i < n; ++i)
    {
        yt[i] = y[i] + hh * dyt[i];
    }
    derivs(xh, yt, dym, k1, k2, omega, r1, shi1, r2, shi2);
    for (int i = 0; i < n; ++i)
    {
        yt[i] = y[i] + h * dym[i];
        dym[i] = dyt[i] + dym[i];
    }
    derivs(x + h, yt, dyt, k1, k2, omega, r1, shi1, r2, shi2);
    for (int i = 0; i < n; ++i)
    {
        yout[i] = y[i] + h6 * (dydx[i] + dyt[i] + 2.0 * dym[i]);
    }
}

int main()
{
    double alpha = 1;
    for (int i = 0; i < n; ++i)
    {
        omega[i] = alpha * tan((i * pi) / static_cast<double>(n) - ((n + 1) * pi) / static_cast<double>(2 * n));
    }
    double A = n * eta;
    for (int i = 0; i < n; ++i)
    {
        theta[i] = -pi + 2 * pi * rand_uniform();
    }
    double k1 = 0;
    ofstream rt_file(""); // name of the file to be opened
    for (int k_loop = 1; k_loop <= 80; ++k_loop)
    {
        double r = 0.0, r_2 = 0.0, t = 0.0;
        for (int it = 1; it <= nite; ++it)
        {
            double rx1 = 0.0, ry1 = 0.0, rx2 = 0.0, ry2 = 0.0;
            for (int i = 0; i < n; ++i)
            {
                rx1 += cos(theta[i]);
                ry1 += sin(theta[i]);
                rx2 += cos(2 * theta[i]);
                ry2 += sin(2 * theta[i]);
            }
            double r1 = sqrt(rx1 * rx1 + ry1 * ry1) / static_cast<double>(n);
            double r2 = sqrt(rx2 * rx2 + ry2 * ry2) / static_cast<double>(n);
            double shi1 = atan2(ry1, rx1);
            double shi2 = atan2(ry2, rx2);
            if (it > tran)
            {
                r += r1;
                r_2 += r2;
            }
            derivs(t, theta, dth, k1, k2, omega, r1, shi1, r2, shi2);
            rk4(theta, dth, n, t, h, tho, k1, k2, omega, r1, shi1, r2, shi2);
            for (int i = 0; i < n; ++i)
            {
                theta[i] = fmod(tho[i], 2 * pi);
            }
            for (int i = 0; i < n; ++i)
            {
                dth_o[i] = dth[i];
            }
            t += h;
        }
        r /= static_cast<double>(nite - tran);
        r_2 /= static_cast<double>(nite - tran);
        rt_file << k1 << " " << r << " " << r_2 << endl;
        cout << k1 << " " << r << " " << r_2 << endl;
        k1 += 0.5;
    }
    return 0;
}
