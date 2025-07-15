#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "funcions.h"

const int N = 40000;
const int steps = 100;
const double dt = 1.;
const double Lx = 100.;
const double Ly = 50.;
const int Nx = 100;
const int Ny = 50;

const double vd = 3.;
const double vth = 1.;


void crossProduct(double *A, double *B, double *AcrossB)
{
    AcrossB[0] = A[1]*B[2] - A[2]*B[1];
    AcrossB[1] = A[2]*B[0] - A[0]*B[2];
    AcrossB[2] = A[0]*B[1] - A[1]*B[0];
}

double mod(double a, double b)
{
    return (a < 0) ? fmod((a+(floor(-a/b)+1)*b), b) : (fmod(a,b));
}

void densitat(double **posicions, double **rho)
{
    double dx = Lx / (double) Nx;
    double dy = Ly / (double) Ny;

    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            rho[i][j] = 0.;
        }
    }

    for(int p = 0; p < N; p++)
    {
        int i = floor(posicions[p][0] / dx);
        int j = floor(posicions[p][1] / dy);

        double pos_rel_x = posicions[p][0] - (double) (i*dx);
        double pos_rel_y = posicions[p][1] - (double) (j*dy);

        int i_plus = (int) mod((double) (i+1), (double) Nx);
        int j_plus = (int) mod((double) (j+1), (double) Ny);

        rho[i][j] += (dx - pos_rel_x) * (dy - pos_rel_y);
        rho[i][j_plus] += (dx - pos_rel_x) * pos_rel_y;
        rho[i_plus][j] += pos_rel_x * (dy - pos_rel_y);
        rho[i_plus][j_plus] += pos_rel_x * pos_rel_y;
    }

    for(int i = 0; i < Nx; i++)
    {   
        for( int j = 0; j < Ny; j++)
        {
            rho[i][j]  *= (Lx * Ly)/ N;
            rho[i][j] /= (dx * dx * dy * dy);
            rho[i][j] -= 1.;
        }
    }
}

void potencial(double **rho, double **phi)
{
    double dx = Lx / (double) Nx;
    double dy = Ly / (double) Ny;

    int kmax = 10000;

    for(int k = 0; k < kmax; k++){

        for(int i = 0; i < Nx; i++)
        {   
            int i_min = (int) mod((double) (i-1), (double) Nx);
            int i_plus = (int) mod((double) (i+1), (double) Nx);

            for(int j = 0; j < Ny; j++)
            {
                int j_min = (int) mod((double) (j-1), (double) Ny);
                int j_plus = (int) mod((double) (j+1),(double) Ny);

                phi[i][j] = -0.25*(rho[i][j]*(dx*dy) - phi[i_min][j] - phi[i][j_min] 
                                - phi[i][j_plus] - phi[i_plus][j]);
            }
        }
    }
}

void campoEnodes(double **phi, vector2d **Enodes)
{   
    double dx = Lx / (double) Nx;
    double dy = Ly / (double) Ny;

    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            Enodes[i][j].x = 0.;
            Enodes[i][j].y = 0.;
        }
    }

    for(int j = 0; j < Ny; j++)
    {
        for(int i = 0; i < Nx; i++)
        {
            int i_plus = (int) mod((double) (i+1), (double) Nx);
            int i_min = (int) mod((double) (i-1), (double) Nx);
            Enodes[i][j].x = (phi[i_min][j] - phi[i_plus][j]) /(2.*dx);
        }
    }

    for( int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            int j_plus = (int) mod((double) (j+1), (double) Ny);
            int j_min = (int) mod((double) (j-1), (double) Ny);
            Enodes[i][j].y = (phi[i][j_min] - phi[i][j_plus]) /(2.*dy);
        }
    }
}

void campEpart(vector2d **Enodes, double **posicions, double **Epart)
{   
    double dx = Lx / (double) Nx;
    double dy = Ly / (double) Ny;

    for( int i = 0; i < N; i++)
    {
        Epart[i][0] = 0.;
        Epart[i][1] = 0.;
    }

    for(int p = 0; p < N; p++)
    {
        int i = (int) floor(posicions[p][0] / dx);
        int j = (int) floor(posicions[p][1] / dy);

        double pos_rel_x = posicions[p][0] - (double) (i*dx);
        double pos_rel_y = posicions[p][1] - (double) (j*dy);

        int i_plus = (int) mod((double) (i+1), (double) Nx);
        int j_plus = (int) mod((double) (j+1), (double) Ny);

        double a = (dx - pos_rel_x) * (dy - pos_rel_y);
        double b = (dx - pos_rel_x) * pos_rel_y;
        double c = pos_rel_x * (dy - pos_rel_y);
        double d = pos_rel_x * pos_rel_y;

        Epart[p][0] = Enodes[i][j].x*a + Enodes[i][j_plus].x*b + Enodes[i_plus][j].x*c + Enodes[i_plus][j_plus].x*d;
        Epart[p][1] = Enodes[i][j].y*a + Enodes[i][j_plus].y*b + Enodes[i_plus][j].y*c + Enodes[i_plus][j_plus].y*d;
    }

    for( int i = 0; i < N; i++)
    {
        Epart[i][0] /= (dx*dy);
        Epart[i][1] /= (dx*dy);

    }
}

// avanza medio tiempo la velocidad
void boris(double **velocitats, double **Epart, double dt)
{
    for(int i = 0; i < N; i++)
    {
        velocitats[i][0] -= 0.5 * Epart[i][0] * dt;
        velocitats[i][1] -= 0.5 * Epart[i][1] * dt;
    }
}


// avanza las posiciones y velocidades un tiempo
void avanze(double **posicions, double **velocitats, double **Epart, double **rho, double **phi, vector2d ** Enodes, double dt)
{   
    // calculo de la aceleracion
    aceleracion(posicions, rho, phi, Epart, Enodes);
    // avanze medio tiempo las velocidades
    boris(velocitats, Epart, dt);

    // avanze un tiempo las posiciones
    for(int i = 0; i < N; i++)
    {
        posicions[i][0] += velocitats[i][0] * dt;
        posicions[i][1] += velocitats[i][1] * dt;
        // condiciones de contorno periodicas
        posicions[i][0] = mod(posicions[i][0], Lx);
        posicions[i][1] = mod(posicions[i][1], Ly);
    }

    // calculo de la aceleracion
    aceleracion(posicions, rho, phi, Epart, Enodes);
    // avanze medio tiempo las velocidades
    boris(velocitats, Epart, dt);
}


// calculo del campo electrico en la posicionn de las particulas ( aceleracion )
void aceleracion(double **posicions, double **rho, double **phi, double **Epart, vector2d ** Enodes)
{
    densitat(posicions, rho);
    potencial(rho, phi);
    campoEnodes(phi, Enodes);
    campEpart(Enodes, posicions,Epart);
}

// numero aleatoria entre 0 y 1
double rand01()
{
    return (double)rand() / (double)RAND_MAX ;
}

// distribucion normal
double norm(double media, double sd)
{
    double   r, theta, x, xn;   

    r = sqrt(-2. * log(rand01()));
    theta = 2. * M_PI * rand01();

    // generamos la variable
    x = r * cos(theta);
    xn = (x * sd) + media;

    return(xn);
}