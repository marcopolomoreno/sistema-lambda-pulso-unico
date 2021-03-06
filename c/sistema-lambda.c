#include <stdio.h>
#include <math.h>

//Desenvolvido por Marco Polo Moreno de Souza em 24/02/2022
//Resolve equações de Bloch no domínio do tempo

double q12 = 0;
double q13 = 0.5*2*3.141592653589793*5e6;
double q33 = 2*3.141592653589793*5e6;

double alpha = 0;
double Tp = 100e-15;
double deltaOmega = 0;
double w1 = -2*3.141592653589793*6835e6*0.5;
double w2 = 2*3.141592653589793*6835e6*0.5;
double phi = 0;
double d = 0;
double Omega = 2*3.141592653589793*1e12;
double t = 0;

double h = 1e-15;

double a11 = 0.5, a22 = 0.5;
double a33 = 0, a12 = 0, b12 = 0;
double a13 = 0, b13 = 0, a23 = 0, b23 = 0;

#define print "%g %g %g %g %g\n", 1e15*t, a11, a22, a33, a11+a22+a33


double bloch(double a11, double a22, double a33, double a12, double b12, 
             double a13, double b13, double a23, double b23, int j)  //sistema de 3 níveis lambda
{
    if (j==1) return 2*Omega*(cos(alpha)*b13-sin(alpha)*a13)                +0.5*q33*a33;
    if (j==2) return 2*Omega*(cos(alpha)*b23-sin(alpha)*a23)                +0.5*q33*a33;
    if (j==3) return 2*Omega*(sin(alpha)*(a13+a23)-cos(alpha)*(b13+b23))        -q33*a33;
    
    if (j==4) return -(w2-w1)*b12+Omega*(cos(alpha)*(b13+b23)-sin(alpha)*(a13+a23))-a12*q12;
    if (j==5) return  (w2-w1)*a12+Omega*(cos(alpha)*(a23-a13)+sin(alpha)*(b23-b13))-b12*q12;
   
    if (j==6) return  -(d-w1)*b13+Omega*cos(alpha)*b12+Omega*sin(alpha)*(a12+a11-a33)-a13*q13;
    if (j==7) return   (d-w1)*a13+Omega*sin(alpha)*b12+Omega*cos(alpha)*(a33-a11-a12)-b13*q13;
    if (j==8) return  -(d-w2)*b23-Omega*cos(alpha)*b12+Omega*sin(alpha)*(a12+a22-a33)-a23*q13;
    if (j==9) return   (d-w2)*a23-Omega*sin(alpha)*b12+Omega*cos(alpha)*(a33-a22-a12)-b23*q13;
}

int main()
{
    FILE *arq;
    arq = fopen("dados.txt", "w");

    double k1[10], k2[10], k3[10], k4[10];
    int k, p;

    fprintf(arq, "tempo rho11 rho22 rho33 soma\n");
    fprintf(arq, "fs\n");

    for (k=0; k<=1000; k++){

        for (p=1; p<=9; p++)
            k1[p] = bloch( a11, a22, a33, a12, b12, a13, b13, a23, b23, p );
        
        for (p=1; p<=9; p++)
            k2[p] = bloch( a11 + 0.5*h*k1[1], a22 + 0.5*h*k1[2], a33 + 0.5*h*k1[3], 
                           a12 + 0.5*h*k1[4], b12 + 0.5*h*k1[5], a13 + 0.5*h*k1[6], 
                           b13 + 0.5*h*k1[7], a23 + 0.5*h*k1[8], b23 + 0.5*h*k1[9], p );

        for (p=1; p<=9; p++)
            k3[p] = bloch( a11 + 0.5*h*k2[1], a22 + 0.5*h*k2[2], a33 + 0.5*h*k2[3], 
                           a12 + 0.5*h*k2[4], b12 + 0.5*h*k2[5], a13 + 0.5*h*k2[6], 
                           b13 + 0.5*h*k2[7], a23 + 0.5*h*k2[8], b23 + 0.5*h*k2[9], p );

        for (p=1; p<=9; p++)
            k4[p] = bloch( a11 + h*k3[1], a22 + h*k3[2], a33 + h*k3[3], 
                           a12 + h*k3[4], b12 + h*k3[5], a13 + h*k3[6], 
                           b13 + h*k3[7], a23 + h*k3[8], b23 + h*k3[9], p );

        t = t + h;
        
        a11 = a11 + h/6.0 * ( k1[1] + 2*k2[1] + 2*k3[1] + k4[1] );
        a22 = a22 + h/6.0 * ( k1[2] + 2*k2[2] + 2*k3[2] + k4[2] );
        a33 = a33 + h/6.0 * ( k1[3] + 2*k2[3] + 2*k3[3] + k4[3] );

        a12 = a12 + h/6.0 * ( k1[4] + 2*k2[4] + 2*k3[4] + k4[4] );
        b12 = b12 + h/6.0 * ( k1[5] + 2*k2[5] + 2*k3[5] + k4[5] );

        a13 = a13 + h/6.0 * ( k1[6] + 2*k2[6] + 2*k3[6] + k4[6] );
        b13 = b13 + h/6.0 * ( k1[7] + 2*k2[7] + 2*k3[7] + k4[7] );

        a12 = a23 + h/6.0 * ( k1[8] + 2*k2[8] + 2*k3[8] + k4[8] );
        b23 = b23 + h/6.0 * ( k1[9] + 2*k2[9] + 2*k3[9] + k4[9] );

        printf(print);
        fprintf(arq, print);
    }

}