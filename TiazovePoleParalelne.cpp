// TiazovePoleParalelne.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Parallel Code

// 1.   P0 -> loads input data
// 2.   MPL Bcast B, L, H, g values;
// 3.   Each P will get their Xx_i, Sx_i, norm_i
// 4.   Each P will get their part of Matrix A -> N x n_local
// 5.   BSGS Solver -> Parallel Solving
//      2x sa nasobi matica x vector -> Parallelne kazdy P da svoju cast a svoj vector
//      Vypis rezidui iba na P0!! 
// 6.   BGather dame dokopy -> MPI_Gather//MPI_Allgather na vsetkych alebo postupne je jedno, dobre je na Vsetkych naraz
//      dostaneme Alpha Vectors z BSGS solving 
//      => cast 5 min
// 7.   Vypocet Vectora U paralelne ->       
//      => cast Cas 30 min
// 8.   P0 potom zoberie vsetky hodnoty vectora U (dame dokopy) a zapise ich do suboru s B,L

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define _USE_MATH_DEFINES 
#include <cmath>

#define R 6378.0
#define D 300
#define GM 398600.5
#define N 8102
#define EPSILON 0.000000000001

using std::cos;
using std::sin;
using std::sqrt;

int main(int argc, char** argv)
{
    double* B = new double[N] {0.0};
    double* L = new double[N] {0.0};
    double* H = new double[N] {0.0};
    double Brad = 0.0, Lrad = 0.0, u2n2 = 0.0; // H=0.0;
    double temp = 0.0;

    // suradnice bodov X_i
    double* X_x = new double[N] {0.0};
    double* X_y = new double[N] {0.0};
    double* X_z = new double[N] {0.0};
    double xNorm = 0.0;

    // suradnice bodov s_j
    double* s_x = new double[N] {0.0};
    double* s_y = new double[N] {0.0};
    double* s_z = new double[N] {0.0};

    // suradnicce normal v x_i
    double* n_x = new double[N] {0.0};
    double* n_y = new double[N] {0.0};
    double* n_z = new double[N] {0.0};

    // g vektor
    double* g = new double[N] {0.0};

    // r vector
    double r_x = 0.0;
    double r_y = 0.0;
    double r_z = 0.0;
    double rNorm = 0.0;
    double rNorm3 = 0.0;

    // dot product of vector r with normal n[i]
    double Kij = 0.0;

    // Load Data
    FILE* file = nullptr;
    file = fopen("C:/Users/karol/Documents/UNI_SCHOOL/STU/Inzinier/Semester_8/Parallel_Algoritmy/TiazovePole/data/BL-8102.dat", "r");

    if (file == nullptr) {
        printf("Error opening data file");
        return -1;
    }

    for (int i = 0; i < N; i++)
    {
        int result = fscanf(file, "%lf %lf %lf %lf %lf", &B[i], &L[i], &H[i], &g[i], &u2n2);

        //Real Data
        g[i] = -g[i] * 0.00001;
        //g[i] = u2n2;
        //g[i] = - GM / (R * R); 

        H[i] = 0.0;

        Brad = B[i] * M_PI / 180.0;
        Lrad = L[i] * M_PI / 180.0;

        X_x[i] = (R + H[i]) * cos(Brad) * cos(Lrad);
        X_y[i] = (R + H[i]) * cos(Brad) * sin(Lrad);
        X_z[i] = (R + H[i]) * sin(Brad);

        s_x[i] = (R + H[i] - D) * cos(Brad) * cos(Lrad);
        s_y[i] = (R + H[i] - D) * cos(Brad) * sin(Lrad);
        s_z[i] = (R + H[i] - D) * sin(Brad);

        xNorm = sqrt(X_x[i] * X_x[i] + X_y[i] * X_y[i] + X_z[i] * X_z[i]);
        n_x[i] = -X_x[i] / xNorm;
        n_y[i] = -X_y[i] / xNorm;
        n_z[i] = -X_z[i] / xNorm;

    }

    fclose(file);

    //for (int i = 0; i < N; i++) // set constant g values
    //    g[i] = -(GM) / (R * R);

    // vytvorenie matice systemu rovnic
    double* A = new double[N * N] {0.0};
    int ij = -1;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            // compute vector r & its norm
            r_x = X_x[i] - s_x[j];
            r_y = X_y[i] - s_y[j];
            r_z = X_z[i] - s_z[j];

            rNorm = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

            rNorm3 = rNorm * rNorm * rNorm;

            // dot product of vector r and normal n_i
            Kij = r_x * n_x[i] + r_y * n_y[i] + r_z * n_z[i];
            //if (i == j && i < 10)
            //    printf("Kij: %.4lf\n", Kij);

            // compute 
            ij = i * N + j;
            A[ij] = (1.0 / (4.0 * M_PI * rNorm3)) * Kij;
        }
    }

    //########## BCGS linear solver ##########//

    double* sol = new double[N]; // vektor x^0 -> na ukladanie riesenia systemu
    double* r_hat = new double[N]; // vektor \tilde{r} = b - A.x^0;
    double* r = new double[N]; // vektor pre rezidua

    double* p = new double[N]; // pomocny vektor na update riesenia
    double* v = new double[N]; // pomocny vektor na update riesenia
    double* s = new double[N]; // pomocny vektor na update riesenia
    double* t = new double[N]; // pomocny vektor na update riesenia
    // Alpha vector
    double* alphaV = new double[N]; // pomocny vektor na update riesenia

    double beta = 0.0;
    double rhoNew = 1.0;
    double rhoOld = 0.0;
    double alpha = 1.0;
    double omega = 1.0;

    double tempDot = 0.0;
    double tempDot2 = 0.0;
    double sNorm = 0.0;

    int MAX_ITER = 1000;
    double TOL = 1.0E-7; //ESILON
    int iter = 1;

    double rezNorm = 0.0;
    for (int i = 0; i < N; i++) // set all to zero
    {
        sol[i] = 0.0;
        p[i] = 0.0; // = 0
        v[i] = 0.0; // = 0
        s[i] = 0.0;
        t[i] = 0.0;
        alphaV[i] = 0.0;

        r[i] = g[i];
        r_hat[i] = g[i];
        rezNorm += r[i] * r[i];

    }

    printf("||r0||: %.6lf\n", sqrt(rezNorm));
    rezNorm = 0.0;

    do
    {
        rhoOld = rhoNew; // save previous rho_{i-2}
        rhoNew = 0.0; // compute new rho_{i-1}
        for (int i = 0; i < N; i++) // dot(r_hat, r)
            rhoNew += r_hat[i] * r[i];

        if (rhoNew == 0.0)
            return -1;

        //if (iter == 1)
        //{
        //    printf("iter 1 setup\n");
        //    for (int i = 0; i < N; i++)
        //        p[i] = r[i];
        //}
        //else
        //{
        beta = (rhoNew / rhoOld) * (alpha / omega);
        for (int i = 0; i < N; i++)
        {
            // update vector p^(i)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);
        }
        //}

        // compute vector v = A.p
        for (int i = 0; i < N; i++)
        {
            v[i] = 0.0;
            for (int j = 0; j < N; j++)
            {
                ij = i * N + j;
                v[i] += A[ij] * p[j];
            }
        }

        // compute alpha
        tempDot = 0.0;
        for (int i = 0; i < N; i++)
            tempDot += r_hat[i] * v[i];

        alpha = rhoNew / tempDot;

        // compute vektor s
        sNorm = 0.0;
        for (int i = 0; i < N; i++)
        {
            alphaV[i] = alpha;      // saving alpha values into vector
            s[i] = r[i] - alpha * v[i];
            sNorm += s[i] * s[i];
        }

        // solution
        sNorm = sqrt(sNorm);
        if (sNorm < TOL) // check if ||s|| is small enough
        {
            for (int i = 0; i < N; i++) // update solution x
                sol[i] = sol[i] + alpha * p[i];

            printf("BCGS stop: ||s|| is small enough, iter: %d\n", iter);
            break;
        }

        // compute vector t = A.s
        for (int i = 0; i < N; i++)
        {
            t[i] = 0.0;
            for (int j = 0; j < N; j++)
            {
                ij = i * N + j;
                t[i] += A[ij] * s[j];
            }
        }

        // compute omega
        tempDot = 0.0; tempDot2 = 0.0;
        for (int i = 0; i < N; i++)
        {
            tempDot += t[i] * s[i];
            tempDot2 += t[i] * t[i];
        }
        omega = tempDot / tempDot2;

        rezNorm = 0.0;
        for (int i = 0; i < N; i++)
        {
            sol[i] = sol[i] + alpha * p[i] + omega * s[i]; // update solution x
            r[i] = s[i] - omega * t[i]; // compute new residuum vector
            rezNorm += r[i] * r[i]; // compute residuum norm
        }

        rezNorm = sqrt(rezNorm);
        printf("iter: %d    ||r||: %.6lf\n", iter, rezNorm);

        if (rezNorm < TOL)
        {
            printf("BCGS stop iter: ||r|| is small enough\n");
            break;
        }

        iter++;

    } while ((iter < MAX_ITER) && (rezNorm > TOL));

    delete[] r_hat;
    delete[] r;
    delete[] p;
    delete[] v;
    delete[] s;
    delete[] t;

    //########## EXPORT ALPHA DATA ##########//
    file = fopen("../outputAlpha.dat", "w");
    if (file == nullptr)
    {
        printf("data export failed\n");
        return -1;
    }

    printf("Vector alpha export started... ");
    for (int i = 0; i < N; i++)
    {
        fprintf(file, "%.5lf\n", alphaV[i]);
    }

    fclose(file);
    printf("done\n");

    //########## EXPORT VECTOR U ##########//
    file = fopen("../outputVectorU.dat", "w");
    if (file == nullptr)
    {
        printf("data export failed\n");
        return -1;
    }

    printf("Vector U export started... ");
    double ui = 0.0, Gij = 0.0;
    for (int i = 0; i < N; i++)
    {
        ui = 0.0;
        for (int j = 0; j < N; j++) // compute solution u(X_i)
        {
            r_x = X_x[i] - s_x[j];
            r_y = X_y[i] - s_y[j];
            r_z = X_z[i] - s_z[j];

            rNorm = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

            Gij = 1.0 / (4.0 * M_PI * rNorm);
            // u_exact = 62.492
            ui += sol[j] * Gij;
        }

        fprintf(file, "%.5lf\n", ui);
    }

    fclose(file);
    printf("done\n");

    //########## EXPORT DATA ##########//
    file = fopen("../outputData.dat", "w");
    if (file == nullptr)
    {
        printf("data export failed\n");
        return -1;
    }

    printf("solution export started... ");
    ui = 0.0, Gij = 0.0;
    for (int i = 0; i < N; i++)
    {
        ui = 0.0;
        for (int j = 0; j < N; j++) // compute solution u(X_i)
        {
            r_x = X_x[i] - s_x[j];
            r_y = X_y[i] - s_y[j];
            r_z = X_z[i] - s_z[j];

            rNorm = sqrt(r_x * r_x + r_y * r_y + r_z * r_z);

            Gij = 1.0 / (4.0 * M_PI * rNorm);

            ui += sol[j] * Gij;
        }

        fprintf(file, "%.5lf\t%.5lf\t%.5lf\n", B[i], L[i], ui);
    }

    fclose(file);
    printf("done\n");

    delete[] X_x;
    delete[] X_y;
    delete[] X_z;
    delete[] s_x;
    delete[] s_y;
    delete[] s_z;
    delete[] n_x;
    delete[] n_y;
    delete[] n_z;
    delete[] g;
    delete[] A;
    delete[] sol;
};