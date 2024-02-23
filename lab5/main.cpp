#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <sstream>
#define _USE_MATH_DEFINES
double delta = 0.2;
int nx = 128;
int ny = 128;
double xmax = delta * nx;
double ymax = delta * ny;
double TOL = 1e-8;

double vb1(double y)
{
    return sin(M_PI * y * delta / ymax);
}

double vb2(double x)
{
    return -1 * sin(2 * M_PI * x * delta / xmax);
}

double vb3(double y)
{
    return sin(M_PI * y * delta / ymax);
}

double vb4(double x)
{
    return sin(2 * M_PI * x * delta / xmax);
}

double calc_S(double **V, int k)
{
    double S = 0;
    for (int i = 0; i < nx; i += k)
    {
        for (int j = 0; j < ny; j += k)
        {
            S += pow((k * delta), 2) / 2.0 *
                 (pow(((V[i + k][j] - V[i][j] + V[i + k][j + k] - V[i][j + k]) / (2.0 * k * delta)), 2) 
                 + pow(((V[i][j + k] - V[i][j] + V[i + k][j + k] - V[i + k][j]) / (2.0 * k * delta)), 2));
        }
    }
    return S;
}

double precision(double s1, double s2)
{
    return fabs((s1 - s2) / s2);
}

void zeros(double **V)
{
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            V[i][j] = 0;
        }
    }
}

void edges(double **V)
{
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            V[0][j] = vb1(j);
            V[i][0] = vb4(i);
            V[nx][j] = vb3(j);
            V[i][ny] = vb2(i);
        }
    }
}

int main()
{
    std::ofstream k16("k16.txt");
    std::ofstream k8("k8.txt");
    std::ofstream k4("k4.txt");
    std::ofstream k2("k2.txt");
    std::ofstream k1("k1.txt");
    std::ofstream files[5];
    std::stringstream iterations;
    int temp = 0;
    int k[] = {16, 8, 4, 2, 1};
    double x[nx + 1];
    double y[ny + 1];
    double **V;
    V = new double *[nx + 1];
    double S = 0;
    double S_prev = 0;
    for (int i = 0; i < nx + 1; i++)
    {
        V[i] = new double[ny + 1];
    }
    zeros(V);
    edges(V);
    int it = 0;
    for (int step : k)
    {   
        iterations.str("");
        iterations<<"it"<<temp<<".txt";
        files[temp].open(iterations.str());
        S = calc_S(V, step);
        while (true)
        {
            it++;
            for (int i = step; i < nx; i += step)
            {
                for (int j = step; j < ny; j += step)
                {
                    V[i][j] = 0.25 * (V[i + step][j] + V[i - step][j] + V[i][j + step] + V[i][j - step]);
                }
            }
            S_prev = S;
            S = calc_S(V, step);
            files[temp] << S << " " << it << std::endl;
            double err = precision(S, S_prev);
            if (err < TOL)
            {
                for (int i = 0; i <= nx; i += step)
                {
                    for (int j = 0; j <= ny; j += step)
                    {
                        if (step == 16)
                            k16 << V[i][j] << " ";
                        if (step == 8)
                            k8 << V[i][j] << " ";
                        if (step == 4)
                            k4 << V[i][j] << " ";
                        if (step == 2)
                            k2 << V[i][j] << " ";
                        if (step == 1)
                            k1 << V[i][j] << " ";
                    }
                    if (step == 16)
                        k16 << std::endl;
                    if (step == 8)
                        k8 << std::endl;
                    if (step == 4)
                        k4 << std::endl;
                    if (step == 2)
                        k2 << std::endl;
                    if (step == 1)
                        k1 << std::endl;
                }
                break;
            }
            edges(V);
        }
        if (step != 1)
        {
            int step2 = step / 2;
            for (int i = 0; i < nx; i += step)
            {
                for (int j = 0; j < ny; j += step)
                {
                    V[i + step2][j + step2] = 0.25 * (V[i][j] + V[i + step][j] + V[i][j + step] + V[i + step][j + step]);
                    V[i + step][j + step2] = 0.5 * (V[i + step][j] + V[i + step][j + step]);
                    V[i + step2][j + step] = 0.5 * (V[i][j + step] + V[i + step][j + step]);
                    V[i + step2][j] = 0.5 * (V[i][j] + V[i + step][j]);
                    V[i][j + step2] = 0.5 * (V[i][j] + V[i][j + step]);
                }
            }
        edges(V);
        }
        std::cout << "po breaku\n";
        temp++;
    }
    std::cout<<it;
    k16.close();
    k8.close();
    k4.close();
    k2.close();
    k1.close();
    return 0;
}