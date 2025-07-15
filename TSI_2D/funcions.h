#ifndef FUNCIONS_H
#define FUNCIONS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

extern const int N;
extern const int steps;
extern const double dt;
extern const double Lx;
extern const double Ly;
extern const int Nx;
extern const int Ny;

extern const double vd;
extern const double vth;

typedef struct Vector2d{
    double x;
    double y;
}vector2d;


void crossProduct(double *A, double *B, double *AcrossB);

double mod(double a, double b);

void densitat(double **posicions, double **rho);

void potencial(double **rho, double **phi);

void campoEnodes( double **phi, vector2d **Enodes);

void campEpart(vector2d **Enodes, double **posicions, double **Epart);

void boris( double **velocitats, double **Epart, double dt);

void avanze(double **posicions, double **velocitats, double **Epart, double **rho, double **phi, vector2d ** Enodes, double dt);

void aceleracion(double **posicions, double **rho, double **phi, double **Epart, vector2d ** Enodes);

double rand01();

double norm(double media, double sd);

#endif