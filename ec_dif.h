//------------------ Archivo  ec_dif.h ------------------------------------
#ifndef EC_DIF_H
#define EC_DIF_H
#include <iostream>
#include <vector>
#include <math.h>
#include "cmath"
#include <fstream>


using namespace std;

//Todas los métodos son para ec. dif de segundo orden y sus argumentos son funcion a integrar con tres variables posicion, velocidad , tiempo, luego t0 tiempo inicial, tf tiempo final, dt paso de tiempo, x0 posicion a t0, y v0 velocidad a t0. 

// metodo de integracion, funcion a integrar , posiones iniciales en 3d, velocidad inicial 3d, (tiempo inicial, tiempo final, dt)
void integrador(vector<double> (*)(double,double,double,double,double), vector<double> (*)(vector<double>,vector<double>,double) , vector<double>,vector<double>, vector<double>  );

void integrador_adaptativo(vector<double> (*)(double,double,double,double,double), vector<double> (*)(vector<double>,vector<double>,double) , vector<double>,vector<double>, vector<double>  ); 


//aceleracion,  
vector<double> euler( double, double,double, double, double ); 



//void euler( double (*)(double, double, double), double , double , double , double , double );

void euler_cromer( double (*)(double, double, double), double , double , double , double , double );

void punto_medio( double (*)(double, double,double), double , double , double , double , double );

void verlet( double (*)(double, double, double), double , double , double , double , double );

void leap_frog( double (*)(double, double, double), double , double , double , double , double );

void RK2( double (*)(double, double,double), double , double , double , double , double );

void RK4( double (*)(double, double,double), double , double , double , double ,double ); 

void RK4_adaptativo( double (*)(double, double,double), double , double , double , double ,double ); 

void euler_predictor_corrector( double (*)(double, double, double), double , double , double , double , double ); 

void PEFRL( double (*)(double, double, double), double , double , double , double , double );


//void integrador( void (*)) 


#endif

