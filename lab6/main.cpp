#include "mgmres.h"
#include <cmath>
#include <iostream>
#include <fstream>

int calcJ(int l, int nx)
{
    return l / (nx + 1);
}

int calcI(int l, int nx)
{
    return l - calcJ(l, nx) * (nx + 1);
}

double rho1(int i, int j, double delta, int nx, int ny)
{
    double xmax = delta * nx;
    double ymax = delta * ny;
    double sigma = xmax / 10.0;

    double first = std::pow((delta * i - 0.25 * xmax) / sigma, 2);
    double second = std::pow((delta * j - 0.5 * ymax) / sigma, 2);

    return std::exp(-first - second);
}
double rho2(int i, int j, double delta, int nx, int ny)
{
    double xmax = delta * nx;
    double ymax = delta * ny;
    double sigma = xmax / 10.0;

    double first = std::pow((delta * i - 0.75 * xmax) / sigma, 2);
    double second = std::pow((delta * j - 0.5 * ymax) / sigma, 2);

    return -std::exp(-first - second);
}
double calcEpsilon(int l, int nx, int epsilon1, int epsilon2)
{
    int i = calcI(l, nx);
    if (i > nx / 2)
        return epsilon2;
    return epsilon1;
}

void fillMatrix(int nx, int ny, double delta, double epsilon1, double epsilon2,
                double V_odd, double V_ev, bool useRho, bool printMatrix, std::string filename)
{
    int N = (nx + 1) * (ny + 1);
    double *a = (double *)malloc(sizeof(double) * 5 * N);
    int ja[5 * N];
    int ia[N + 1];
    double b[N], V[N];
    for (int i = 0; i < N * 5; i++)
    {
        a[i] = 0;
        ja[i] = 0;
        if (i < N + 1)
            ia[i] = -1;
        if (i < N)
        {
            b[i] = 0;
            V[i] = 0;
        }
    }

    int k = -1;
    int edge;
    double vb;
    int j;
    int i;
    double delta_sqr = std::pow(delta, 2);
    for (int l = 0; l < N; l++)
    {
        edge = 0;
        vb = 0.0;
        j = calcJ(l, nx);
        i = calcI(l, nx);
        // lewy brzeg
        if (i == 0)
        {
            edge = 1;
            vb = V_odd;
        }
        // górny brzeg
        if (j == ny)
        {
            edge = 1;
            vb = V_ev;
        }
        // prawy brzeg
        if (i == nx)
        {
            edge = 1;
            vb = V_odd;
        }
        // dolny brzeg
        if (j == 0)
        {
            edge = 1;
            vb = V_ev;
        }

        if (edge == 1)
        {
            b[l] = vb;
        }
        else if (useRho)
        {
            b[l] = -(rho1(i, j, delta, nx, ny) + rho2(i, j, delta, nx, ny));
        }
        else
        {
            b[l] = 0;
        }
        ia[l] = -1;
        // lewa skrajna przekatna
        if (l - nx - 1 >= 0 && edge == 0)
        {
            k++;
            if (ia[l] < 0)
            {
                ia[l] = k;
            }
            a[k] = calcEpsilon(l, nx, epsilon1, epsilon2) / delta_sqr;
            ja[k] = l - nx - 1;
        }
        // poddiagonala
        if (l - 1 >= 0 && edge == 0)
        {
            k++;
            if (ia[l] < 0)
            {
                ia[l] = k;
            }
            a[k] = calcEpsilon(l, nx, epsilon1, epsilon2) / delta_sqr;
            ja[k] = l - 1;
        }
        // diagonala
        k++;
        if (ia[l] < 0)
        {
            ia[l] = k;
        }
        if (edge == 0)
        {
            double eps2 = calcEpsilon(l, nx, epsilon1, epsilon2);
            double eps1 = calcEpsilon(l + 1, nx, epsilon1, epsilon2);
            double eps1nx = calcEpsilon(l + nx + 1, nx, epsilon1, epsilon2);
            a[k] = -(2 * eps2 + eps1 + eps1nx) / delta_sqr;
        }
        else
        {
            a[k] = 1;
        }
        ja[k] = l;
        // naddiagonala
        if (l < N && edge == 0)
        {
            k++;
            a[k] = calcEpsilon(l + 1, nx, epsilon1, epsilon2) / delta_sqr,
            ja[k] = l + 1;
        }
        // prawa skrajna przekątna
        if (l < N - nx - 1 && edge == 0)
        {
            k++;
            a[k] = calcEpsilon(l + nx + 1, nx, epsilon1, epsilon2) / delta_sqr;
            ja[k] = l + nx + 1;
        }
    }
    int nz_num = k + 1;
    ia[N] = nz_num;
    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, 500, 500, 1e-8, 1e-8);
    if (printMatrix)
    {
        std::ofstream matrixA("matrixA.txt");
        std::ofstream vectorB("vectorB.txt");
        for (int l = 0; l < 5 * N; l++)
        {
            matrixA << l << " " << a[l] << "\n";
        }
        for (int l = 0; l < N; l++)
        {
            vectorB << l << " " << calcI(l, nx) << " " << calcJ(l, nx) << " " << b[l] << "\n";
        }
    }
    else
    {
        std::ofstream file(filename);
        double endline = 0;
        for (int l = 0; l < N; l++)
        {
            file <<l <<" "<< V[l] << "\n";
        }
    }
    free(a);
}

int main()
{
    int itr_max = 500;
    int mr = 500;
    int tol_abs = 1e-8;
    int tol_rel = 1e-8;
    fillMatrix(4, 4, 0.1, 1, 1, 10, -10, false, true, "");
    fillMatrix(50, 50, 0.1, 1, 1, 10, -10, false, false, "5a.txt");
    fillMatrix(100, 100, 0.1, 1, 1, 10, -10, false, false, "5b.txt");
    fillMatrix(200, 200, 0.1, 1, 1, 10, -10, false, false, "5c.txt");
    fillMatrix(100, 100, 0.1, 1, 1, 0, 0, true, false, "6a.txt");
    fillMatrix(100, 100, 0.1, 1, 2, 0, 0, true, false, "6b.txt");
    fillMatrix(100, 100, 0.1, 1, 10, 0, 0, true, false, "6c.txt");
    return 0;
}
