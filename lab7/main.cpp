#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

double delta = 0.01;
int rho = 1;
int mu = 1;
const int nx = 200;
const int ny = 90;
int i_1 = 50;
int j_1 = 55;
int IT_MAX = 20000;

double qwy(double qwe)
{
    return qwe * (std::pow((ny * delta - j_1 * delta), 3)) / std::pow(delta * ny, 3);
}

bool isEdge(int i, int j)
{
    if (j == ny)
    {
        // b
        return true;
    }
    if (i == nx)
    {
        // c
        return true;
    }
    if (i == 0 && (j >= j_1 && j <= ny))
    {
        // a
        return true;
    }
    if ((i >= i_1 && i <= nx) && j == 0)
    {
        // d
        return true;
    }
    // if (i == i_1 && (j >= 0 && j <= j_1))
    if (i <= i_1 && j <= j_1)
    {
        // e
        return true;
    }
    // if ((i >= 0 && i <= i_1) && j == j_1)
    // {
    //     // f
    //     return true;
    // }
    return false;
}
double calcError(double **psi, double **zeta)
{
    double error = 0;
    int j2 = j_1 + 2;
    for (int i = 1; i < nx; i++)
    {
        error += psi[i + 1][j2] + psi[i - 1][j2] + psi[i][j2 + 1] + psi[i][j2 - 1] - 4 * psi[i][j2] - delta * delta * zeta[i][j2];
    }
    return error;
}
void wbPsi(double **psi, double qwe)
{
    // brzeg a
    for (int j = j_1; j <= ny; j++)
    {
        psi[0][j] = qwe / (2.0 * mu) * (std::pow(j * delta, 3) / 3.0 - std::pow(j * delta, 2) / 2.0 * (j_1 * delta + ny * delta) + j * delta * j_1 * delta * ny * delta);
    }
    // brzeg c
    for (int j = 0; j <= ny; j++)
    {
        psi[nx][j] = qwy(qwe) / (2 * mu) * (std::pow(j * delta, 3) / 3.0 - std::pow(j * delta, 2) / 2.0 * ny * delta) + (qwe * std::pow(delta * j_1, 2) * (-j_1 * delta + 3 * ny * delta)) / (12 * mu);
    }
    for (int i = 1; i < nx; i++)
    {
        // brzeg b
        psi[i][ny] = psi[0][ny];
        if (i >= i_1)
        {
            // brzeg d
            psi[i][0] = psi[0][j_1];
        }
    }
    for (int j = 1; j <= j_1; j++)
    {
        // brzeg e
        psi[i_1][j] = psi[0][j_1];
    }
    for (int i = 1; i <= i_1; i++)
    {

        // brzeg f
        psi[i][j_1] = psi[0][j_1];
    }
}
void wbZeta(double **psi, double **zeta, double qwe)
{
    for (int j = 0; j <= ny; j++)
    {
        if (j >= j_1)
        {
            // brzeg a
            zeta[0][j] = qwe / (2.0 * mu) * (2 * delta * j - delta * j_1 - delta * ny);
        }
        if (j >= 1 && j < j_1)
        {
            // brzeg e
            zeta[i_1][j] = 2 / (delta * delta) * (psi[i_1 + 1][j] - psi[i_1][j]);
        }
        // brzeg c
        zeta[nx][j] = qwy(qwe) / (2.0 * mu) * (2 * j * delta - delta * ny);
    }
    for (int i = 1; i < nx; i++)
    {
        if (i > i_1)
        {
            // brzeg d
            zeta[i][0] = 2.0 / (delta * delta) * (psi[i][1] - psi[i][0]);
        }
        // brzeg b
        zeta[i][ny] = 2.0 / (delta * delta) * (psi[i][ny - 1] - psi[i][ny]);
        if (i <= i_1)
        {
            // brzeg f
            zeta[i][j_1] = 2. / (delta * delta) * (psi[i][j_1 + 1] - psi[i][j_1]);
        }
    }
    zeta[i_1][j_1] = 0.5 * (zeta[i_1 - 1][j_1] + zeta[i_1][j_1 - 1]);
}
void relaxation(double **psi, double **zeta, double **u, double **v, int qwe)
{
    wbPsi(psi, qwe);
    int omega;
    double error;
    for (int it = 1; it <= IT_MAX; it++)
    {
        if (it < 2000)
        {
            omega = 0;
        }
        else
        {
            omega = 1;
        }
        for (int i = 1; i < nx; i++)
        {
            for (int j = 1; j < ny; j++)
            {
                if (!isEdge(i, j))
                {
                    psi[i][j] = 0.25 * (psi[i + 1][j] + psi[i - 1][j] + psi[i][j + 1] + psi[i][j - 1] - delta * delta * zeta[i][j]);
                    zeta[i][j] = 0.25 * (zeta[i + 1][j] + zeta[i - 1][j] + zeta[i][j + 1] + zeta[i][j - 1]) - omega * rho / (16.0 * mu) * ((psi[i][j + 1] - psi[i][j - 1]) * (zeta[i + 1][j] - zeta[i - 1][j]) - (psi[i + 1][j] - psi[i - 1][j]) * (zeta[i][j + 1] - zeta[i][j - 1]));
                    u[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2.0 * delta);
                    v[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2.0 * delta);
                }
            }
        }
        wbZeta(psi, zeta, qwe);
        error = calcError(psi, zeta);
    }
    std::cout<<error<<std::endl;
}
int main()
{
    double **psi = new double *[nx + 1];
    for (int i = 0; i < nx + 1; ++i)
        psi[i] = new double[ny + 1]();

    double **zeta = new double *[nx + 1];
    for (int i = 0; i < nx + 1; ++i)
        zeta[i] = new double[ny + 1]();

    double **u = new double *[nx + 1];
    for (int i = 0; i < nx + 1; ++i)
        u[i] = new double[ny + 1]();

    double **v = new double *[nx + 1];
    for (int i = 0; i < nx + 1; ++i)
        v[i] = new double[ny + 1]();
    std::string filename = "-1000";
    std::ofstream fileZeta("zeta" + filename + ".txt");
    std::ofstream filePsi("psi" + filename + ".txt");
    std::ofstream fileU("u" + filename + ".txt");
    std::ofstream fileV("v" + filename + ".txt");
    relaxation(psi, zeta, u, v, -1000);
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            filePsi << psi[i][j] << " ";
            fileZeta << zeta[i][j] << " ";
            fileU << u[i][j] << " ";
            fileV << v[i][j] << " ";
        }
        filePsi << "\n";
        fileZeta << "\n";
        fileU << "\n";
        fileV << "\n";
    }
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            psi[i][j] = 0;
            zeta[i][j] = 0;
            u[i][j] = 0;
            v[i][j] = 0;
        }
    }
    filename = "4000";
    std::ofstream fileZeta1("zeta" + filename + ".txt");
    std::ofstream filePsi1("psi" + filename + ".txt");
    std::ofstream fileU1("u" + filename + ".txt");
    std::ofstream fileV1("v" + filename + ".txt");
    relaxation(psi, zeta, u, v, 4000);
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            filePsi1 << psi[i][j] << " ";
            fileZeta1 << zeta[i][j] << " ";
            fileU1 << u[i][j] << " ";
            fileV1 << v[i][j] << " ";
        }
        filePsi1 << "\n";
        fileZeta1 << "\n";
        fileU1 << "\n";
        fileV1 << "\n";
    }
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            psi[i][j] = 0;
            zeta[i][j] = 0;
            u[i][j] = 0;
            v[i][j] = 0;
        }
    }
    filename = "-4000";
    std::ofstream fileZeta2("zeta" + filename + ".txt");
    std::ofstream filePsi2("psi" + filename + ".txt");
    std::ofstream fileU2("u" + filename + ".txt");
    std::ofstream fileV2("v" + filename + ".txt");
    relaxation(psi, zeta, u, v, -4000);
    for (int i = 0; i < nx + 1; i++)
    {
        for (int j = 0; j < ny + 1; j++)
        {
            filePsi2 << psi[i][j] << " ";
            fileZeta2 << zeta[i][j] << " ";
            fileU2 << u[i][j] << " ";
            fileV2 << v[i][j] << " ";
        }
        filePsi2 << "\n";
        fileZeta2 << "\n";
        fileU2 << "\n";
        fileV2 << "\n";
    }
    for (int i = 0; i < nx + 1; ++i)
    {
        delete[] psi[i];
        delete[] u[i];
        delete[] v[i];
        delete[] zeta[i];
    }
    delete[] psi;
    delete[] zeta;
    delete[] u;
    delete[] v;
    return 0;
}