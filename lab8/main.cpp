#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

using namespace std;
const int nx = 400;
const int ny = 90;
const int i_1 = 200;
const int i_2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = delta * 10;
const double xa = 0.45;
const double ya = 0.45;

int main()
{
    vector<int> i_data, j_data;
    vector<int> u0, u1;
    vector<double> temp_psi;
    double psi[nx + 1][ny + 1];
    double vx[nx + 1][ny + 1];
    double vy[nx + 1][ny + 1];
    double vmax;

    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            vx[i][j] = 0;
            vy[i][j] = 0;
        }
    }

    ifstream f("psi.dat");
    string line;
    int i = 0;

    while (getline(f, line))
    {
        stringstream ss(line);
        string str;
        while (getline(ss, str, '\t'))
        {
            if (i % 3 == 0)
            {
                i_data.push_back(stoi(str));
            }
            else if (i % 3 == 1)
            {
                j_data.push_back(stoi(str));
            }
            else
            {
                temp_psi.push_back(stod(str));
            }
            i++;
        }
    }
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            psi[i][j] = temp_psi.at(i * ny + j);
            cout<<psi[i][j]<<endl;
        }
    }
    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            vx[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (delta * 2.0);
            vy[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (delta * 2.0);
        }
    }
    // Zastawka
    for (int i = i_1; i <= i_2; i++)
    {
        for (int j = 0; j <= j_1; j++)
        {
            vx[i][j] = 0.0;
            vy[i][j] = 0.0;
        }
    }
    // Gorny i dolny brzeg
    for (int i = 1; i < nx; i++)
    {
        vx[i][0] = 0.0;
        vy[i][ny] = 0.0;
    }

    // Lewy i prawy brzeg
    for (int j = 0; j <= ny; j++)
    {
        vx[0][j] = vx[1][j];
        vx[nx][j] = vx[nx - 1][j];
    }

    for (int i = 1; i < nx; i++)
    {
        for (int j = 1; j < ny; j++)
        {
            // cout<<vy[i][j]<<endl;
        }
    }
    vmax = sqrt(pow(vx[0][0], 2) + pow(vy[0][0], 2));
    double vmax_temp;
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            vmax_temp = sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2));
            if (vmax_temp > vmax)
            {
                vmax = vmax_temp;
            }
        }
    }
    double dt = delta / (4.0 * vmax);
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            u0.push_back(exp(-(pow(i * delta - xa, 2) + pow(j * delta - ya, 2)) / (2 * sigma * sigma)) / (2 * M_PI * sigma * sigma));
        }
    }
    double d;
    for (int d_temp = 0; d_temp <= 1; d_temp++)
    {
        d = d_temp * 0.1;
        cout<<d<<endl;
        for (int it = 1; it <= 100; it++)
        {
            u1 = u0;
            for (int k = 1; k <= 20; k++)
            {
                for (int i = 0; i <= nx; i++)
                {
                    for (int j = 1; j <= ny - 1; j++)
                    {
                        if ((i >= i_1 && i <= i_2) && (j >= 0 && j <= j_1))
                        {
                            continue;
                        }
                        else if (i == 0)
                        {
                            u1.at(i * ny + j) = 
                            (1.0 / (1 + 2 * d * dt / (delta * delta))) 
                            * (u0.at(i * ny + j) - dt / 2.0 * vx[i][j] * 
                            ((u0.at((i + 1) * ny + j) - u0.at(nx * ny + j)) / (delta * 2.0) 
                            + (u1.at((i + 1) * ny + j) - u1.at(nx * ny + j)) / (delta * 2.0)) 
                            - dt * vy[i][j] / 2.0 
                            * ((u0.at(i * ny + j + 1) - u0.at(i * ny + j - 1)) / (delta * 2.0) 
                            + (u1.at(i * ny + j + 1) - u1.at(i * ny + j - 1)) / (delta * 2.0)) 
                            + dt * d / 2.0 
                            * ((u0.at((i + 1) * ny + j) + u0.at(nx * ny + j) + u0.at(i * ny + j + 1) + u0.at(i * ny + j - 1) - 4 * u0.at(i*ny+j))/(delta*delta)
                            + (u1.at((i + 1) * ny + j) + u1.at(nx * ny + j) + u1.at(i * ny + j + 1) + u1.at(i * ny + j - 1))/(delta*delta)));
                        }
                        else if (i == nx)
                        {
                            u1.at(i * ny + j) = 
                            (1.0 / (1 + 2 * d * dt / (delta * delta))) 
                            * (u0.at(i * ny + j) - dt / 2.0 * vx[i][j] * 
                            ((u0.at(0 * ny + j) - u0.at((i - 1) * ny + j)) / (delta * 2.0) 
                            + (u1.at(0 * ny + j) - u1.at((i - 1) * ny + j)) / (delta * 2.0)) 
                            - dt * vy[i][j] / 2.0 
                            * ((u0.at(i * ny + j + 1) - u0.at(i * ny + j - 1)) / (delta * 2.0) 
                            + (u1.at(i * ny + j + 1) - u1.at(i * ny + j - 1)) / (delta * 2.0)) 
                            + dt * d / 2.0 
                            * ((u0.at(0 * ny + j) + u0.at((i - 1) * ny + j) + u0.at(i * ny + j + 1) + u0.at(i * ny + j - 1) - 4 * u0.at(i*ny+j))/(delta*delta)
                            + (u1.at(0 * ny + j) + u1.at((i - 1) * ny + j) + u1.at(i * ny + j + 1) + u1.at(i * ny + j - 1))/(delta*delta)));
                        }
                        else
                        {
                            u1.at(i * ny + j) = (1.0 / (1 + 2 * d * dt / (delta * delta))) 
                            * (u0.at(i * ny + j) - dt / 2.0 * vx[i][j] * 
                            ((u0.at((i + 1) * ny + j) - u0.at((i - 1) * ny + j)) / (delta * 2.0) 
                            + (u1.at((i + 1) * ny + j) - u1.at((i - 1) * ny + j)) / (delta * 2.0)) 
                            - dt * vy[i][j] / 2.0 * ((u0.at(i * ny + j + 1) - u0.at(i * ny + j - 1)) / (delta * 2.0) + (u1.at(i * ny + j + 1) - u1.at(i * ny + j - 1)) / (delta * 2.0)) 
                            + dt * d / 2.0 * ((u0.at((i + 1) * ny + j) + u0.at((i - 1) * ny + j) + u0.at(i * ny + j + 1) + u0.at(i * ny + j - 1) - 4 * u0.at(i*ny+j))/(delta*delta)
                            + (u1.at((i + 1) * ny + j) + u1.at((i - 1) * ny + j) + u1.at(i * ny + j + 1) + u1.at(i * ny + j - 1))/(delta*delta)));
                        }
                    }
                }
            }
            u0 = u1;
            double c = 0;
            double xsr = 0;
            for(int i = 0; i < nx; i++){
                for(int j = 0; j <ny;j++){
                    c+= u0.at(i*ny+j)*delta*delta;
                    xsr+= i*delta * u0.at(i*ny+j)*delta*delta;
                }
            }
            // cout<<c<<" ";
        }
    }
    return 0;
}