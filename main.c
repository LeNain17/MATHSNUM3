#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 100
#define MU 0.5
#define TMAX 2.5

void creationMAT_C(double **MAT_C, int NB_LIGNES, int NB_COLONNES){

    for(int i = 0; i<NB_LIGNES; i++){
        for(int j = 0; j<NB_COLONNES; j++){
            MAT_C[i][j] = 0;
        }
    }

    for(int i = 0; i<NB_LIGNES; i++){
        for(int j = 0; j<NB_COLONNES; j++){
            if(i ==j){
                MAT_C[i][j] = 1 - 2*MU;
            }
            if(j == i+1 || i == j+1){
                MAT_C[i][j] = MU;
            }
            if(i == j+1){
                MAT_C[i][j] = MU;
            }
            if(i == NX-1 && j == NX-2){
                MAT_C[i][j] = 2*MU;
            }
        }
    }

}

void creationVEC_SE(double *VEC_SE, int NB_LIGNES){
    for(int i = 0; i<NB_LIGNES; i++){
        VEC_SE[i] = 0;
        if(i == 0){
            VEC_SE[i] = MU;
        }
    }
}

void creationVEC_U(double* VEC_U, double** MAT_C, double* VEC_SE, int NB_LIGNES, int NB_COLONNES, double* VEC_U0, double pasTemps, double* VEC_xmon, int NB_POINTS, double* VEC_tmon, int NB_TEMPS){
    
    double valeur_SommeProduit = 0;

    //INITIALISATION VECTEUR COPIE DE U
    double *VECopie_U = malloc(sizeof(double)*NB_LIGNES);

    double Temps = 0.0;
    printf("TEMPS ini : %f\n", Temps);

    while(Temps<2.5){

        for(int i=1; i<NB_LIGNES; i++){
            valeur_SommeProduit = 0;
            if(Temps == 0){
                for(int j = 0; j<NB_COLONNES; j++){
                    valeur_SommeProduit = valeur_SommeProduit + MAT_C[i][j]*VEC_U0[j];
                }
            }
            else{
                for(int j = 0; j<NB_COLONNES; j++){
                    valeur_SommeProduit = valeur_SommeProduit + MAT_C[i][j]*VECopie_U[j];
                    //printf("TEMPS : %d, MAT_C[%d][%d] = %f, VEC_U[%d] = %f; val SP = %f\n", Temps, i, j, MAT_C[i][j], j, VECopie_U[j], valeur_SommeProduit);
                }
            }

            VEC_U[i] = valeur_SommeProduit + VEC_SE[i];
        }

        VEC_U[0] = 1;
        VEC_U[NB_LIGNES-1] = 0;

        for(int k = 0; k<NB_LIGNES; k++){
            VECopie_U[k] = VEC_U[k];
        }

        FILE* fichierUt = NULL;
        fichierUt = fopen("variaUt.txt", "a");

        for(int ku = 0; ku<NB_LIGNES; ku++){
            for(int kx = 0; kx<NB_POINTS; kx++){
                //fprintf(fichier, "KU : %d div : %f KX : %d, vec_xmon[kx] = %f\n", ku, ((ku*1.0/(NX*1.0))*1.0), kx, VEC_xmon[kx]);
                if((ku*1.0/(NX*1.0))*1.0 == VEC_xmon[kx]){
                    //fprintf(fichier, "p%.2f : %f; ku : %d, kx : %d", VEC_xmon[kx], VEC_U[ku], ku, kx);
                    fprintf(fichierUt, "%f   ", VEC_U[ku]);
                }
            }
        }

        fprintf(fichierUt, "%f\n", Temps);
        fclose(fichierUt);

        FILE* fichierUx = NULL;
        fichierUx = fopen("variaUx.txt", "a");

        for(int tt = 0; tt<NB_TEMPS; tt++){
            if(fabs(VEC_tmon[tt] - Temps) < 1e-6){
                for(int tu = 0; tu<NB_LIGNES; tu++){
                    fprintf(fichierUx, "%.2f %f   ", (tu*1.0)/(NX*1.0),VEC_U[tu]);
                }
                fprintf(fichierUx, "%f\n", Temps);
            }
        }

        fclose(fichierUx);

        Temps = Temps + pasTemps;
    }
}

void AfficherMAT(double **MAT, int NB_LIGNES, int NB_COLONNES){
    printf("MATRICE; \n");
    FILE* fichier = NULL;
    fichier = fopen("MatC.txt", "a");
    for(int i = 0; i<NB_LIGNES; i++){
        for(int j = 0; j<NB_COLONNES; j++){
           fprintf(fichier, "   %f", MAT[i][j]);
        }
        fprintf(fichier, "\n");
    }
    fclose(fichier);
}

void AfficherVEC(double *VEC, int NB_LIGNES){
    printf("VECTEUR; \n");
    FILE* fichier = NULL;
    fichier = fopen("MatU.txt", "a");
    for(int i = 0; i<NB_LIGNES; i++){
        //printf("%f", VEC[i]);
        //printf("\n");
        fprintf(fichier, "%f", VEC[i]);
        fprintf(fichier, "\n");
    }
    fclose(fichier);
}


/*----------------------------------------------------------*/
/*---------------------------MAIN---------------------------*/
/*----------------------------------------------------------*/

int main(){ 

    int NB_LIGNES = NX;
    int NB_COLONNES = NX;

    double pasEspace = 1.0/NX;
    double pasTemps = MU*pasEspace*pasEspace*1.0;

    //suppression des anciennes valeurs du fichier "CoordoF.txt"
    FILE* fichierSUPPRmatC = NULL;
    fichierSUPPRmatC = fopen("MatC.txt", "w");
    fclose(fichierSUPPRmatC);

    //suppression des anciennes valeurs du fichier "CoordoF.txt"
    FILE* fichierSUPPRmatU = NULL;
    fichierSUPPRmatU = fopen("MatU.txt", "w");
    fclose(fichierSUPPRmatU);

    //suppression des anciennes valeurs du fichier "CoordoF.txt"
    FILE* fichierSUPPRvariaUt = NULL;
    fichierSUPPRvariaUt = fopen("variaUt.txt", "w");
    fclose(fichierSUPPRvariaUt);

    //suppression des anciennes valeurs du fichier "CoordoF.txt"
    FILE* fichierSUPPRvariaUx = NULL;
    fichierSUPPRvariaUx = fopen("variaUx.txt", "w");
    fclose(fichierSUPPRvariaUx);

    printf("%f\n", pasTemps);
    printf("pas espace : %f\n", pasEspace);

    //INITIALISATION MATRICE C
    double **MAT_C = malloc(sizeof(*MAT_C) * NB_LIGNES);
        for (int i = 0; i < NB_LIGNES; i++){
            MAT_C[i] = malloc(sizeof(**MAT_C) * NB_COLONNES);
        }

    //INITIALISATION VECTEUR SE
    double *VEC_SE = malloc(sizeof(double)* NB_LIGNES);

    //INITIALISATION VECTEUR U
    double *VEC_U = malloc(sizeof(double)*NB_LIGNES);

    //INITIALISATION VECTEUR U0
    double *VEC_U0 = malloc(sizeof(double)*NB_LIGNES);

    for(int i = 0; i<NB_LIGNES; i++){
        VEC_U0[i] = 0;
    }
    
    int NB_POINTS = 5;
    double *VEC_xmon = malloc(sizeof(double)*NB_POINTS);
    
    double valeurs_xmon[] = {0.1, 0.2, 0.25, 0.5, 0.75};

    for(int w = 0; w < NB_POINTS; w++){
        VEC_xmon[w] = valeurs_xmon[w]*1.0;
        //printf("vec_xmon[%d] = %f", w, VEC_xmon[w]);
    }

    int NB_TEMPS = 5;
    double *VEC_tmon = malloc(sizeof(double)*NB_TEMPS);

    double valeurs_tmon[] = {0.20, 0.40, 0.60, 0.80, 2.00};

    for(int wt = 0; wt < NB_TEMPS; wt++){
        VEC_tmon[wt] = valeurs_tmon[wt]*1.0;
        printf("vect_tmon[%d] = %f  ", wt, VEC_tmon[wt]);
    }

    creationVEC_SE((double *)VEC_SE, NB_LIGNES);
    //AfficherVEC((double*) VEC_SE, NB_LIGNES);

    creationMAT_C((double**) MAT_C, NB_LIGNES, NB_COLONNES);
    AfficherMAT((double**) MAT_C, NB_LIGNES, NB_COLONNES);

    creationVEC_U((double*) VEC_U, (double**)MAT_C, (double*)VEC_SE, NB_LIGNES, NB_COLONNES, (double*)VEC_U0, pasTemps, (double*) VEC_xmon, NB_POINTS, (double*) VEC_tmon, NB_TEMPS);
    AfficherVEC((double*)VEC_U, NB_LIGNES);

    free(VEC_SE);
    free(VEC_U);
    free(VEC_U0);
    free(VEC_xmon);
    return 0;
}
