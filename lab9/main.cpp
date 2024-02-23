#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

const int nx = 40;
const int ny = 40;
const int N = (nx + 1) * (ny + 1);
int delta = 1;
int dt = 1;
int Ta = 40;
int Tb = 0;
int Tc = 30;
int Td = 0;
double kb = 0.1;
double kd = 0.6;
int IT_MAX = 2000;

int calc_l(int i, int j)
{
    return i + j * (nx + 1);
}

int main()
{
    gsl_matrix *A = gsl_matrix_calloc(N, N);
    gsl_matrix *B = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);
    gsl_vector *T_prev = gsl_vector_calloc(N);
    gsl_permutation *perm = gsl_permutation_alloc(N);

    gsl_vector_set_zero(T);
    gsl_vector_set_zero(c);
    gsl_vector_set_zero(d);
    gsl_matrix_set_zero(A);
    gsl_matrix_set_zero(B);

    std::ofstream t100("t100.txt");
    std::ofstream t200("t200.txt");
    std::ofstream t500("t500.txt");
    std::ofstream t1000("t1000.txt");
    std::ofstream t2000("t2000.txt");

    std::ofstream dt100("dt100.txt");
    std::ofstream dt200("dt200.txt");
    std::ofstream dt500("dt500.txt");
    std::ofstream dt1000("dt1000.txt");
    std::ofstream dt2000("dt2000.txt");

    int signum;
    // wnetrze
    for (int i = 0; i <= nx; i++)
    {
        for (int j = 0; j <= ny; j++)
        {
            int l = calc_l(i, j);
            // std::cout<<i<<" "<<j<<" "<<l<<std::endl;
            // wnetrze

            if (i > 0 && i < nx && j > 0 && j < ny)
            {
                gsl_matrix_set(A, l, l - nx - 1, dt / (2.0 * delta * delta));
                gsl_matrix_set(A, l, l - 1, dt / (2.0 * delta * delta));
                gsl_matrix_set(A, l, l + 1, dt / (2.0 * delta * delta));
                gsl_matrix_set(A, l, l + nx + 1, dt / (2.0 * delta * delta));
                gsl_matrix_set(A, l, l, -2.0 * dt / delta / delta - 1.0);

                gsl_matrix_set(B, l, l - nx - 1, -dt / (2.0 * delta * delta));
                gsl_matrix_set(B, l, l - 1, -dt / (2.0 * delta * delta));
                gsl_matrix_set(B, l, l + 1, -dt / (2.0 * delta * delta));
                gsl_matrix_set(B, l, l + nx + 1, -dt / (2.0 * delta * delta));
                gsl_matrix_set(B, l, l, 2.0 * dt / delta / delta - 1.0);
            }
            // wb dirichleta
            if (i == 0 || i == nx)
            {
                gsl_matrix_set(A, l, l, 1);
                gsl_matrix_set(B, l, l, 1);
                gsl_vector_set(c, l, 0);
            }
            // wb neumanna gorny brzeg
            if (i > 0 && i < nx && j == ny)
            {
                gsl_matrix_set(A, l, l - nx - 1, -1.0 / (kb * delta));
                gsl_matrix_set(A, l, l, 1.0 + 1.0 / (kb * delta));
                gsl_vector_set(c, l, Tb);
                for (int l_temp = 0; l_temp < N; l_temp++)
                {
                    gsl_matrix_set(B, l, l_temp, 0);
                }
            }
            // wb neumanna dolny brzeg
            if (i > 0 && i < nx && j == 0)
            {
                gsl_matrix_set(A, l, l + nx + 1, -1.0 / (kd * delta));
                gsl_matrix_set(A, l, l, 1 + 1.0 / (kd * delta));
                gsl_vector_set(c, l, Td);
                for (int l_temp = 0; l_temp < N; l_temp++)
                {
                    gsl_matrix_set(B, l, l_temp, 0);
                }
            }
            // wektor t
            if (i == 0)
            {
                gsl_vector_set(T, l, Ta);
            }
            if (i == nx)
            {
                gsl_vector_set(T, l, Tc);
            }
        }
    }

    gsl_linalg_LU_decomp(A, perm, &signum);
    for (int i = 0; i <= IT_MAX; i++)
    {
        if (i % 100 == 0)
        {
            for (int l = 0; l < N; l++)
            {
                gsl_vector_set(T_prev, l, gsl_vector_get(T, l));
            }
        }
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0.0, d);
        gsl_blas_daxpy(1.0, d, c);
        gsl_linalg_LU_solve(A, perm, d, T);
        if (i == 100)
        {
            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    int l = calc_l(i, j);
                    t100 << gsl_vector_get(T, l) << std::endl;
                    dt100 << (gsl_vector_get(T, l) - gsl_vector_get(T_prev, l)) / dt << std::endl;
                }
            }
        }
        if (i == 200)
        {
            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    int l = calc_l(i, j);
                    t200 << gsl_vector_get(T, l) << std::endl;
                    dt200 << (gsl_vector_get(T, l) - gsl_vector_get(T_prev, l)) / dt << std::endl;
                }
            }
        }
        if (i == 500)
        {
            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    int l = calc_l(i, j);
                    t500 << gsl_vector_get(T, l) << std::endl;
                    dt500 << (gsl_vector_get(T, l) - gsl_vector_get(T_prev, l)) / dt << std::endl;
                }
            }
        }
        if (i == 1000)
        {
            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    int l = calc_l(i, j);
                    t1000 << gsl_vector_get(T, l) << std::endl;
                    dt1000 << (gsl_vector_get(T, l) - gsl_vector_get(T_prev, l)) / dt << std::endl;
                }
            }
        }
        if (i == 2000)
        {
            for (int i = 0; i <= nx; i++)
            {
                for (int j = 0; j <= ny; j++)
                {
                    int l = calc_l(i, j);
                    t2000 << gsl_vector_get(T, l) << std::endl;
                    dt2000 << (gsl_vector_get(T, l) - gsl_vector_get(T_prev, l)) / dt << std::endl;
                }
            }
        }
    }
    gsl_permutation_free(perm);
    gsl_vector_free(T);
    gsl_vector_free(T_prev);
    gsl_vector_free(c);
    gsl_vector_free(d);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
}