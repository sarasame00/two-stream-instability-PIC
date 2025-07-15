#include "funcions.h"

int main()
{

    double dx = Lx / (double) Nx;
    double dy = Ly / (double) Ny;

    double margin = dx / 10.;

    double **posicions;
    posicions = (double**) calloc(N, sizeof(double*));
    for(int i = 0; i < N; i++) posicions[i] = calloc(2, sizeof(double)); 

    double **velocitats;
    velocitats = (double**) calloc(N, sizeof(double*));
    for(int i = 0; i < N; i++) velocitats[i] = calloc(2, sizeof(double));

    double **velocitats2;
    velocitats2 = (double**) calloc(N, sizeof(double*));
    for(int i = 0; i < N; i++) velocitats2[i] = calloc(2, sizeof(double));  

    double **rho;
    rho = (double**) calloc(Nx, sizeof(double*));
    for(int i = 0; i < Nx; i++) rho[i] = calloc(Ny, sizeof(double)); 

    double **phi;
    phi = (double**) calloc(Nx, sizeof(double*));
    for(int i = 0; i < Nx; i++) phi[i] = calloc(Ny, sizeof(double)); 

    for(int i = 0; i < Nx; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            phi[i][j] = 0.;
        }
    }

    vector2d **Enodes;
    Enodes = (vector2d**) calloc(Nx, sizeof(vector2d*));
    for(int i = 0; i < Nx; i++) Enodes[i] = calloc(Ny, sizeof(vector2d)); 

    double **Epart;
    Epart = (double**) calloc(N, sizeof(double*));
    for(int i = 0; i < N; i++) Epart[i] = calloc(2, sizeof(double));

    FILE *fileP2;  FILE *fileV2;
    FILE *fileP1;  FILE *fileV1;

	fileP1 = fopen("datos/TSI_2D_50_pos1.txt", "w");
    fileV1 = fopen("datos/TSI_2D_50_vel1.txt", "w");
    fileP2 = fopen("datos/TSI_2D_50_pos2.txt", "w");
    fileV2 = fopen("datos/TSI_2D_50_vel2.txt", "w");

    //FILE *filePot;  FILE *fileDens; 
    FILE *fileE;

    //fileDens = fopen("datos/TSI_2D_50_dens.txt", "w");
	//filePot = fopen("datos/TSI_2D_50_pot.txt", "w");
    fileE = fopen("datos/TSI_2D_50_E.txt", "w");

    //posicion de las particulas
    for(int i = 0; i < N; i++)
    {
        posicions[i][0] = rand01()*(Lx-2.*margin) + margin;
        posicions[i][1] = rand01()*(Ly-2.*margin) + margin;
    }

    //velocidad de las particulas
    for(int i = 0; i < N; i++)
    {   
        velocitats[i][0] = norm(vd, vth);
        velocitats[i][1] = 0.;
        if(i%2 == 0)  velocitats[i][0] *= - 1.;   
    }

    aceleracion(posicions, rho, phi, Epart, Enodes);
    boris(velocitats, Epart, -dt);

    //double EU = 0.;

    //fprintf(filePot, "%f", EU);

    for(int step = 0; step < steps; step++)
    {   
        printf("t = %d\n", step);

       //avanzar un tiempo las posiciones y velocidades
        avanze(posicions, velocitats, Epart, rho, phi, Enodes, dt);

        for(int i = 0; i < N; i++)
        {
            for(int j = 0; j < 2; j++)
            {
                velocitats2[i][j] = velocitats[i][j];
            }
        }
        

        //avanzamos medio tiempo las velocidades para imprimir los datos
        boris(velocitats2, Epart, dt);

        ///guardar datos en los ficheros
        fprintf(fileP2,"%f",posicions[0][0]); fprintf(fileV2,"%f",velocitats2[0][0]);
        fprintf(fileP1,"%f",posicions[1][0]); fprintf(fileV1,"%f",velocitats2[1][0]);
        for(int n = 2; n < N; n++)
        {   
			if(n%2==0) fprintf(fileP2,"\t%f",posicions[n][0]);
            if(n%2!=0) fprintf(fileP1,"\t%f",posicions[n][0]);

			if(n%2==0) fprintf(fileV2,"\t%f",velocitats[n][0]);
            if(n%2!=0) fprintf(fileV1,"\t%f",velocitats[n][0]);
		}
        fprintf(fileV1,"\n"); fprintf(fileP1,"\n");
        fprintf(fileV2,"\n"); fprintf(fileP2,"\n");

        //imprimir valores potencial
        if(step== 30){
            for(int i = 0; i < Nx; i++)
            {
                for(int j = 0; j < Ny; j++)
                {
                    //fprintf(fileDens,"%f\t",rho[i][j]);
                    fprintf(fileE,"%f\t",sqrt(Enodes[i][j].x*Enodes[i][j].x + Enodes[i][j].y*Enodes[i][j].y));
                    //EU += phi[i][j];
                } 
                //fprintf(fileDens,"%f\n",rho[i][Nx-1]);
                fprintf(fileE,"%f\n", sqrt(Enodes[i][Nx-1].x*Enodes[i][Nx-1].x + Enodes[i][Nx-1].y*Enodes[i][Nx-1].y));
            } 
            //fprintf(filePot, "\t%f", EU);
        }


    }

    fclose(fileP1); fclose(fileV1);fclose(fileP2); fclose(fileV2); fclose(fileE);
    //fclose(fileDens), fclose(filePot); 
    return 0;
}

