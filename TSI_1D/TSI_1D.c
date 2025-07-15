#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

int N, cells, L, timesteps;
double K0;

double rand01();
double norm(double media, double sd);
void aceleracion(double *dens, double *part, double *Ei, double *Ej, double *phi);
void dominioPeriodico(double *part);
void densitat(double *mesh, double *posicions);
void poisson1d(double *Ei, double *dens, double *phi);
void campoElectrico(double *Ej, double *Ei, double *posicions);
int mod(int x, int n);

int main(){

    N = 40000;
    cells = 800;
    L = 100;
    timesteps = 100;

    double dt = 1.;

    int vb  = 3;  //beam velocity
	int vth = 1;  //beam width
    double A = 0.; //perturbacion


    double part[N*2];

    // posicion de las particulas
    for(int i = 0; i < N; i++) part[i] = rand01()*L;

    // velocidad de las part
    for(int i = N; i < N*2; i++)
    {   
        part[i] = norm(vb, vth) * (1 + A*sin(2*M_PI*part[i-N]/L));

        // la mitad de las particulas con velocidad negativa
        if(i%2 == 0)  part[i] *= - 1;   
    }

    
    FILE *fileP1;  FILE *fileV1; FILE *fileP2;  FILE *fileV2;
    FILE *filePot; FILE *fileK; FILE *filePotK;

	fileP1 = fopen("datos/TSI_1D_pos1.txt", "w");
    fileV1 = fopen("datos/TSI_1D_vel1.txt", "w");
    fileP2 = fopen("datos/TSI_1D_pos2.txt", "w");
    fileV2 = fopen("datos/TSI_1D_vel2.txt", "w");

    filePot = fopen("datos/TSI_1D_potencial.txt", "w");
    filePotK = fopen("datos/TSI_1D_pot_K.txt", "w");
    fileK = fopen("datos/TSI_1D_cinetica.txt", "w");

    //guardar datos inciales en los ficheros
    fprintf(fileP2,"%f",part[0]); fprintf(fileV2,"%f",part[N]);
    fprintf(fileP1,"%f",part[1]); fprintf(fileV1,"%f",part[N+1]);
    for(int n = 2; n < N; n++)
    {
        if(n%2==0) fprintf(fileP2,"\t%f",part[n]);
        if(n%2!=0) fprintf(fileP1,"\t%f",part[n]);

        if(n%2==0) fprintf(fileV2,"\t%f",part[N+n]);
        if(n%2!=0) fprintf(fileV1,"\t%f",part[N+n]);
    }
    fprintf(fileV1,"\n"); fprintf(fileP1,"\n");
    fprintf(fileV2,"\n"); fprintf(fileP2,"\n");

    double dens[cells], Ei[cells];
    double Ej[N];
    double phi_j[N], phi[N];

    //energia cinetica inicial
    K0 = 0.;
    for(int i = N; i < N*2; i++)
    {
        K0 += 0.5 * part[i] * part[i];
    }
    fprintf (fileK, "%f", 1.);

    // loop temporal
    for(int k = 0; k < timesteps; k++)
    {   
        // calculo de la aceleracion (E en la posicion de las particulas)
        aceleracion(dens, part, Ei, Ej, phi);

        // avanze medio tiempo en la velocidad y un tiempo en la posicion
        for(int j = N; j < N*2; j++) part[j] -= (dt/2.)*Ej[j-N];
        for(int j = 0; j < N; j++) part[j] += dt*part[N+j];

        // comprueba que las part estan dentro el dominio (contorno periodico)
        dominioPeriodico(part);

        // calculo de la aceleracion (E en la posicion de las particulas)
        aceleracion(dens, part, Ei, Ej, phi);

        campoElectrico(phi_j, phi, part);

        // avanze medio tiempo en la velocidad 
        for(int j = N; j < N*2; j++) part[j] -=(dt/2.)*Ej[j-N];

        // comprueba que las part estan dentro el domino (contorno periodico)
        dominioPeriodico(part);

        //guardar datos en los ficheros
        fprintf(fileP2,"%f",part[0]); fprintf(fileV2,"%f",part[N]);
        fprintf(fileP1,"%f",part[1]); fprintf(fileV1,"%f",part[N+1]);
        for(int n = 2; n < N; n++)
        {
			if(n%2==0) fprintf(fileP2,"\t%f",part[n]);
            if(n%2!=0) fprintf(fileP1,"\t%f",part[n]);

			if(n%2==0) fprintf(fileV2,"\t%f",part[N+n]);
            if(n%2!=0) fprintf(fileV1,"\t%f",part[N+n]);
		}
        fprintf(fileV1,"\n"); fprintf(fileP1,"\n");
        fprintf(fileV2,"\n"); fprintf(fileP2,"\n");

        double EU = 0.;
        for(int i = 0; i < N; i++) 
        {   
            EU -= phi_j[i];
        };
        fprintf(filePotK, "%f\n", EU/K0);

        EU = 0.;
        for(int i = 0; i < cells; i++) 
        {   
            fprintf(filePot, "%f\t", phi[i]);
            EU += phi[i];
        };
        fprintf(filePot, "%f\n", EU);
        

        // se imprimen los valores de la energia cinetica
        double K = 0.;
        for(int i = N; i < N*2; i++)
        {
            K += 0.5 * part[i] * part[i];
        }
        fprintf (fileK, "\t%f", K/K0);
    }
    fclose(fileP1); fclose(fileV1); fclose(fileP2); fclose(fileV2);
    fclose(filePot); fclose(fileK);
    return 0;
}

//===================================================//
//                 NUMERO ALEATORIO                  //
//---------------------------------------------------//
// Funcion que un numero aletorio en el intervalo    //
// [0,1].                                            //
//===================================================//

double rand01()
{
    return (double)rand() / (double)RAND_MAX ;
}

//===================================================//
//              DISTRIBUCION NORMAL                  //
//---------------------------------------------------//
// Funcion que un valor aleatorio siguiendo una      //
// distribucion normal usando el metodo Box Muller   //
// - media: la media de la distribucion.             //
// - sd: desviacion estandar(en este caso vth)       //
//===================================================//
double norm(double media, double sd)
{
    double   r, theta, x, xn;   

    r = sqrt(-2.0 * log(rand01()));
    theta = 2.0 * M_PI * rand01();

    // generamos la variable
    x = r * cos(theta);
    xn = (x * sd) + media;

    return(xn);
}

//===================================================//
//              CALCULO ACELERACION                  //
//---------------------------------------------------//
// Funcion que dads las posiciones de las particulas //
// calcula el campo electrico en esta, que, cambiado //
// de signo, es la aceleracion de la particula.      //
// - dens: densidad de carga en los nodos            //
// - part: posiciones de las particulas              //
// - Ei: campo electrico en los nodos                //
// - Ej: campo electrico en la posicion de las       //
//      particulas                                   //
// - filePot: archivo donde se guardan los valores   //
//      del potencial                                //
//===================================================//
void aceleracion(double *dens, double *part, double *Ei, double *Ej, double *phi)
{
    densitat(dens, part);
    poisson1d(Ei, dens, phi);
    campoElectrico(Ej, Ei, part);
}

//===================================================//
//                DOMINIO PERIODICO                  //
//---------------------------------------------------//
// Funcion que comprueba que las particulas esten    //
// dentro del dominio y las desplaza cuando es       //
// necesario.                                        //
// - part: posiciones de las particulas              //
//===================================================//
void dominioPeriodico(double *part)
{
    for(int i = 0; i < N; i++)
        {
            if(part[i] > (double) L) part[i] -= (double) L;
            if(part[i] < 0.) part[i] += (double) L;
        }
}
//===================================================//
//            DENSIDAD EN LOS NODOS                  //
//---------------------------------------------------//
// Funcion que calcula la densidad en los modos a    //
// partir de las posiciones de las particulas.       //
// - mesh: array donde se guardaran los valores de   //
//      la densidad.                                 //
// - posicions: array con las posiciones [:N-1] y    //
//      las velocidades [N:N*2-1] de las particulas. //
//===================================================//

void densitat(double *mesh, double *posicions)
{
    double dx = (double) L / cells;

    for(int i = 0; i < cells; i++) mesh[i] = 0.;

    for(int j = 0; j < N; j++)
    {   
        int i = (int) floor(posicions[j]/dx);
        int i_plus = (int) mod((i+1),cells);
        
        if(i < 0 || i > cells)
        {
            //printf("se va\n");
        }else{
            mesh[i] += ((i+1)*dx - posicions[j])/dx;
            mesh[i_plus] += (posicions[j]-i*dx)/dx;
        }
    }

    //normalizamos
    for(int i = 0; i < cells; i++)
    {
        mesh[i] *= (double) L /(N*dx); 
    }
}

//===================================================//
//            CAMPO ELECTRICO NODOS                  //
//---------------------------------------------------//
// Funcion que soluciona la ecuacion de Poisson      //
// usando el metodo de Gauss-Seidel y a partir del   //
// potencial calcula el campo electrico.             //
// INPUT:                                            //
// - Ei: array donde se guardaran los valores  del   //
//      campo electrico  en los nodos.               //
// - dens: array con las densidades en los nodos.    //
//===================================================//

void poisson1d(double *Ei, double *dens, double *phi)
{
    double dx = (double)  L / cells;
    double rho[cells], xold[cells], xnew[cells];

    for (int i = 0; i < cells; i++)
    {
        dens[i] -= 1.;
    }

    for(int n = 0; n < cells; n++) 
    {
        xold[n] = 0.;
        xnew[n] = 0.;
    };
    //int kmax = 100000;

    //for(int k = 0; k < kmax; k++){
    double err = 0;

    while(err==0||err>(pow(10,-6)))//iteracions fins que err < 10^-6
    {
        err = 0;

        for(int j = 0; j < cells; j++)
        {   
            int j_min = (int) mod((j-1),cells);
            int j_plus = (int) mod((j+1),cells);
            xnew[j] = - 0.5*(dens[j]*(dx*dx) - xnew[j_min] - xold[j_plus]);
        };

        //Passem els valors de xnew a xold
        for(int m = 0; m < cells; m++)
        {
            err += sqrt(pow(xold[m]-xnew[m],2));
            xold[m] = xnew[m];
        };
        err /= cells;//noramlitzem error
    }

    for(int n = 0; n < cells; n++) 
    {   
        phi[n] = xold[n];
    };

    
    // encontrar el campo electrico en los nodos a partir del potencial
    for (int j = 0; j < cells; j++)
    {
        int j_plus = (int) mod((j+1),cells);
        int j_min = (int) mod((j-1),cells);

        Ei[j] = (phi[j_min] - phi[j_plus])/(2.*dx);     
    }

}


//===================================================//
//            CAMPO ELECTRICO particulas             //
//---------------------------------------------------//
// Funcion que interpola el campo electrico en la    //
// posicionn de las particulas a partir del de los   //
// nodos.                                            //
// INPUT:                                            //
// - Ej: array donde se guardara el campo electrico  //
//      en la posicion de las part.            //
// - Ei: array con los valores  del campo electrico  //
//      en los nodos                                 //
// - posicions: array con las posiciones [:N-1] y   //
//      las velocidades [N:N*2-1] de las particulas. //
//===================================================//

void campoElectrico(double *Ej, double *Ei, double *posicions)
{
    double dx = (double) L / cells;

    for(int i = 0; i < N; i++) Ej[i] = 0.;

    for(int p = 0; p < N; p++)
    {   
        int i = (int) floor(posicions[p]/dx);
        int j = (int) floor(posicions[p]/dx);

        int i_min = (int) mod((i-1), cells);
        int i_plus = (int) mod((i+1), cells);

        if( i < 0 || i > cells|| j < 0 || j > cells)
        {
            //printf("se va\n");
        }else{
            Ej[p] += ((((i+1)*dx) - posicions[p]) / dx) * Ei[i] +
                        ((posicions[p] - (i*dx)) / dx) * Ei[i_plus];
        }  
    }
}

int mod(int x,int n){
    return (x % n + n) %n;
}
