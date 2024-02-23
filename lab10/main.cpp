#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

int nx = 150;
int nt = 1000;
double delta = 0.1;
double dt = 0.05;
double xa = 7.5;
double sigma = 0.5;
double tmax = nt * dt;

using namespace std;

double kronecker(double x, double y)
{
    if (x == y)
        return 1.0;
    return 0.0;
}

double adf(double x, double t, double t_max, double xF)
{
    return cos(50 * t / t_max) * kronecker(x, xF);
}

double calculateEnergy(std::vector<double> u, std::vector<double> v)
{
    double temp1 = pow((u[1] - u[0]) / delta, 2);
    double temp2 = pow((u[nx] - u[nx - 1]) / delta, 2);
    double sum = 0;
    for (int i = 1; i <= nx - 1; i++)
    {
        sum += v[i] * v[i] + pow((u[i + 1] - u[i - 1]) / (2 * delta), 2);
    }
    return delta / 4.0 * (temp1 + temp2) + delta / 2.0 * sum;
}
void Verlet(std::vector<double> u0, std::vector<double> u, std::vector<double> v, std::vector<double> vp, std::vector<double> a, double alpha, double beta)
{
    std::ofstream file_u("ualpha" + std::to_string(alpha) + "beta" + std::to_string(beta) + ".txt");
    std::ofstream file_e("ealpha" + std::to_string(alpha) + "beta" + std::to_string(beta) + ".txt");
    double xF = 0;
    if (alpha == 1.0)
    {
        xF = 2.5;
        std::cout << "jestem";
    }
    for (int i = 0; i < nx; i++)
    {
        double x = i * delta;
        if (alpha != 1.0)
        {
            u[i] = exp(-pow((x - xa), 2) / (2.0 * sigma * sigma));
        }
        else
        {
            u[i] = 0;
        }
    }
    u0 = u;
    for (int i = 1; i < nx - 1; i++)
    {
        a[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / (delta * delta) - beta * (u[i] - u0[i]) / dt + alpha * adf(i * delta, 0, tmax, xF);
    }
    for (int step = 1; step <= nt; step++)
    {
        for (int i = 0; i < nx; i++)
        {
            vp[i] = v[i] + dt / 2.0 * a[i];
        }
        u0 = u;
        for (int i = 0; i < nx; i++)
        {
            u[i] = u[i] + dt * vp[i];
        }
        for (int i = 0; i < nx; i++)
        {
            a[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / (delta * delta) - beta * (u[i] - u0[i]) / dt + alpha * adf(i * delta, step*dt, tmax, xF);
            v[i] = vp[i] + dt * a[i] / 2;
        }
        u[0] = 0;
        u[nx - 1] = 0;
        v[0] = 0;
        v[nx - 1] = 0;
        file_e << calculateEnergy(u, v) << std::endl;

        for (int i = 0; i < nx; i++)
        {
            file_u << u[i] << ' ';
        }
        file_u << '\n';
    }
}

int main()
{
    std::vector<double> u0(nx, 0.0);
    std::vector<double> u(nx, 0.0);
    std::vector<double> v(nx, 0.0);
    std::vector<double> vp(nx, 0.0);
    std::vector<double> a(nx, 0.0);
    Verlet(u0, u, v, vp, a, 0, 0);
    for (int i = 0; i < nx; i++)
    {
        u0[i] = u[i] = v[i] = vp[i] = a[i] = 0;
    }

    Verlet(u0, u, v, vp, a, 0, 0.1);
    for (int i = 0; i < nx; i++)
    {
        u0[i] = u[i] = v[i] = vp[i] = a[i] = 0;
    }

    Verlet(u0, u, v, vp, a, 0, 1);
    for (int i = 0; i < nx; i++)
    {
        u0[i] = u[i] = v[i] = vp[i] = a[i] = 0;
    }

    Verlet(u0, u, v, vp, a, 1, 1);
}
