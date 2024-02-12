#include <iostream>
#include <cmath>
#include <vector>
#include <chrono>

using namespace std;

void rk4(float alpha2, float lambda2, float Dx, float r1, float r2, int ndim, float *omega1, float *omega2, float lambda1, float *theta1in, float *theta2in, float *dydt1, float *dydt2, float t, float h, float *theta1out, float *theta2out, float shi1, float shi2, float rho1, float rho2, float phi1, float phi2, float lambda3);

void derivs(float alpha2, float lambda2, float Dx, float r1, float r2, int ndim, float t, float *omega1, float *omega2, float lambda1, float *theta1in, float *theta2in, float *dydt1, float *dydt2, float shi1, float shi2, float rho1, float rho2, float phi1, float phi2, float lambda3);

/**
 * Simulates a symplectic all-to-all system with given parameters.
 *
 * @param None
 *
 * @return None
 *
 * @throws None
 */
void symplectic_all_to_all()
{
    int ndim = 1000;
    vector<float> omega1(ndim), omega2(ndim);
    float omega_mean1, omega_mean2;
    vector<float> theta1in(ndim), theta2in(ndim);
    float alpha, alpha2;
    vector<float> theta1out(ndim), theta2out(ndim), dydt1(ndim), dydt2(ndim);
    vector<vector<int>> a2(ndim, vector<int>(ndim));
    vector<vector<int>> a1(ndim, vector<int>(ndim));
    float h, t, r1_final, r2_final, shi1, shi2, lambda1_step;
    float lambda1, lambda2, Dx, coup1, coup2, add1, add2, lambda1_max;
    float rx1, rx2, ry1, ry2, r1, r2, x;
    float rho1, rho2, phi1, phi2, lambda3, lambda1_min;
    float rhox1, rhox2, rhoy1, rhoy2, rho1_final, rho2_final;
    float pi;

    h = 0.01;
    int nstep = 20000;
    t = 0.0;
    int itrans = 10000;

    pi = 4.0 * atan(1.0);

    int val = 1;
    float a = 1.0;
    float b = 0.0;
    lambda2 = 8.0;
    lambda3 = 0.0;
    lambda1_step = 0.2;
    lambda1_max = 3.0;
    lambda1_min = 1.0;

    alpha = 1.0;
    for (int i = 0; i < ndim; i++)
    {
        omega1[i] = alpha * tan((i * pi) / float(ndim) - (float(ndim + 1) * pi) / float(2 * ndim));
        omega2[i] = alpha * tan((i * pi) / float(ndim) - (float(ndim + 1) * pi) / float(2 * ndim));
        omega_mean1 += omega1[i];
        omega_mean2 += omega2[i];
    }
    omega_mean1 /= float(ndim);
    omega_mean2 /= float(ndim);

    for (int i = 0; i < ndim; i++)
    {
        theta1in[i] = -pi + 2 * pi * rand();
        theta2in[i] = 2 * pi;
    }

    r1 = 0.0;
    r2 = 0.0;
    r1_final = 0.0;
    r2_final = 0.0;
    rho1 = 0.0;
    rho2 = 0.0;
    rho1_final = 0.0;
    rho2_final = 0.0;

    for (int it = 1; it <= nstep; it++)
    {
        rx1 = 0;
        rx2 = 0;
        ry1 = 0.0;
        ry2 = 0.0;
        rhox1 = 0;
        rhox2 = 0;
        rhoy1 = 0.0;
        rhoy2 = 0.0;

        for (int i = 0; i < ndim; i++)
        {
            rx1 += cos(theta1in[i]);
            ry1 += sin(theta1in[i]);
            rx2 += cos(2.0 * theta1in[i]);
            ry2 += sin(2.0 * theta1in[i]);
            rhox1 += cos(theta2in[i]);
            rhoy1 += sin(theta2in[i]);
            rhox2 += cos(2.0 * theta2in[i]);
            rhoy2 += sin(2.0 * theta2in[i]);
        }

        r1 = sqrt(rx1 * rx1 + ry1 * ry1) / float(ndim);
        r2 = sqrt(rx2 * rx2 + ry2 * ry2) / float(ndim);
        rho1 = sqrt(rhox1 * rhox1 + rhoy1 * rhoy1) / float(ndim);
        rho2 = sqrt(rhox2 * rhox2 + rhoy2 * rhoy2) / float(ndim);
        shi1 = atan(ry1 / rx1);
        shi2 = atan(ry2 / rx2);
        phi1 = atan(rhoy1 / rhox1);
        phi2 = atan(rhoy2 / rhox2);

        rk4(alpha2, lambda2, Dx, r1, r2, ndim, omega1.data(), omega2.data(), lambda1, theta1in.data(), theta2in.data(), dydt1.data(), dydt2.data(), t, h, theta1out.data(), theta2out.data(), shi1, shi2, rho1, rho2, phi1, phi2, lambda3);

        for (int i = 0; i < ndim; i++)
        {
            theta1in[i] = fmod(theta1out[i], 2.0 * pi);
            theta2in[i] = fmod(theta2out[i], 2.0 * pi);
        }

        if (it > itrans)
        {
            r1_final += r1;
            r2_final += r2;
            rho1_final += rho1;
            rho2_final += rho2;
        }
    }

    r1_final /= float(nstep - itrans);
    r2_final /= float(nstep - itrans);
    rho1_final /= float(nstep - itrans);
    rho2_final /= float(nstep - itrans);

    cout << lambda1 << " " << r1_final << " " << r2_final << endl;
}

/**
 * Perform fourth-order Runge-Kutta integration.
 *
 * @param alpha2 coefficient for the second-order term in the equation
 * @param lambda2 coefficient for the first-order term in the equation
 * @param Dx increment in the x-direction
 * @param r1 first parameter
 * @param r2 second parameter
 * @param ndim dimension of the system
 * @param omega1 pointer to the first frequency array
 * @param omega2 pointer to the second frequency array
 * @param lambda1 coefficient for the first-order term in the equation
 * @param theta1in pointer to the initial angle array for the first pendulum
 * @param theta2in pointer to the initial angle array for the second pendulum
 * @param dydt1 pointer to the derivative of the angle array for the first pendulum
 * @param dydt2 pointer to the derivative of the angle array for the second pendulum
 * @param t time
 * @param h step size
 * @param theta1out pointer to the final angle array for the first pendulum
 * @param theta2out pointer to the final angle array for the second pendulum
 * @param shi1 parameter 1
 * @param shi2 parameter 2
 * @param rho1 parameter 3
 * @param rho2 parameter 4
 * @param phi1 parameter 5
 * @param phi2 parameter 6
 * @param lambda3 coefficient for the third-order term in the equation
 *
 * @throws None
 */
void rk4(float alpha2, float lambda2, float Dx, float r1, float r2, int ndim, float *omega1, float *omega2, float lambda1, float *theta1in, float *theta2in, float *dydt1, float *dydt2, float t, float h, float *theta1out, float *theta2out, float shi1, float shi2, float rho1, float rho2, float phi1, float phi2, float lambda3)
{
    int nmax = 1000;
    float h6, hh, th;
    vector<float> dym1(nmax), dyt1(nmax), yt1(nmax);
    vector<float> dym2(nmax), dyt2(nmax), yt2(nmax);
    vector<vector<int>> a2(ndim, vector<int>(ndim));
    vector<vector<int>> a1(ndim, vector<int>(ndim));
    h6 = h / 6.0;
    hh = h * 0.50;
    th = t + hh;

    derivs(alpha2, lambda2, Dx, r1, r2, ndim, t, omega1, omega2, lambda1, theta1in, theta2in, dydt1, dydt2, shi1, shi2, rho1, rho2, phi1, phi2, lambda3);

    for (int i = 0; i < ndim; i++)
    {
        yt1[i] = theta1in[i] + hh * dydt1[i];
        yt2[i] = theta2in[i] + hh * dydt2[i];
    }

    derivs(alpha2, lambda2, Dx, r1, r2, ndim, th, omega1, omega2, lambda1, yt1.data(), yt2.data(), dyt1.data(), dyt2.data(), shi1, shi2, rho1, rho2, phi1, phi2, lambda3);

    for (int i = 0; i < ndim; i++)
    {
        yt1[i] = theta1in[i] + hh * dyt1[i];
        yt2[i] = theta2in[i] + hh * dyt2[i];
    }

    derivs(alpha2, lambda2, Dx, r1, r2, ndim, th, omega1, omega2, lambda1, yt1.data(), yt2.data(), dym1.data(), dym2.data(), shi1, shi2, rho1, rho2, phi1, phi2, lambda3);

    for (int i = 0; i < ndim; i++)
    {
        yt1[i] = theta1in[i] + h * dym1[i];
        dym1[i] = dyt1[i] + dym1[i];
        yt2[i] = theta2in[i] + h * dym2[i];
        dym2[i] = dyt2[i] + dym2[i];
    }

    derivs(alpha2, lambda2, Dx, r1, r2, ndim, t + h, omega1, omega2, lambda1, yt1.data(), yt2.data(), dyt1.data(), dyt2.data(), shi1, shi2, rho1, rho2, phi1, phi2, lambda3);

    for (int i = 0; i < ndim; i++)
    {
        theta1out[i] = theta1in[i] + h6 * (dydt1[i] + dyt1[i] + 2.0 * dym1[i]);
        theta2out[i] = theta2in[i] + h6 * (dydt2[i] + dyt2[i] + 2.0 * dym2[i]);
    }
}

/**
 * Calculate the derivatives for the given parameters.
 *
 * @param alpha2 description of parameter
 * @param lambda2 description of parameter
 * @param Dx description of parameter
 * @param r1 description of parameter
 * @param r2 description of parameter
 * @param ndim description of parameter
 * @param t description of parameter
 * @param omega1 description of parameter
 * @param omega2 description of parameter
 * @param lambda1 description of parameter
 * @param theta1in description of parameter
 * @param theta2in description of parameter
 * @param dydt1 description of parameter
 * @param dydt2 description of parameter
 * @param shi1 description of parameter
 * @param shi2 description of parameter
 * @param rho1 description of parameter
 * @param rho2 description of parameter
 * @param phi1 description of parameter
 * @param phi2 description of parameter
 * @param lambda3 description of parameter
 *
 * @return description of return value
 *
 * @throws ErrorType description of error
 */
void derivs(float alpha2, float lambda2, float Dx, float r1, float r2, int ndim, float t, float *omega1, float *omega2, float lambda1, float *theta1in, float *theta2in, float *dydt1, float *dydt2, float shi1, float shi2, float rho1, float rho2, float phi1, float phi2, float lambda3)
{
    float sigma1 = 0.0;
    float sigma2 = 0.0;
    float sigma3 = 0.0;

    for (int i = 0; i < ndim; i++)
    {
        dydt1[i] = omega1[i] + lambda1 * r1 * sin(shi1 - theta1in[i]) + lambda2 * r1 * r2 * sin(shi2 - shi1 - theta1in[i]) + lambda3 * r1 * r1 * rho1 * sin(shi1 - theta1in[i]);
        dydt2[i] = omega2[i] + sigma1 * rho1 * sin(phi1 - theta2in[i]) + sigma2 * rho1 * rho2 * sin(phi2 - phi1 - theta2in[i]) + sigma3 * rho1 * rho1 * rho1 * sin(phi1 - theta2in[i]);
    }
}

int main()
{

    auto start = chrono::high_resolution_clock::now();
    symplectic_all_to_all();
    auto end = chrono::high_resolution_clock::now();

    chrono::duration<double> duration = end - start;
    cout << "Runtime: " << duration.count() << " seconds" << endl;

    return 0;
}
