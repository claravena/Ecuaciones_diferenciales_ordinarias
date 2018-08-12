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

void euler( double (*funcion)(double, double, double), double , double , double , double , double );

void euler_cromer( double (*funcion)(double, double, double), double , double , double , double , double );

void punto_medio( double (*funcion)(double, double,double), double , double , double , double , double );

void verlet( double (*funcion)(double, double, double), double , double , double , double , double );

void leap_frog( double (*funcion)(double, double, double), double , double , double , double , double );

void RK2( double (*funcion)(double, double,double), double , double , double , double , double );

void RK4( double (*funcion)(double, double,double), double , double , double , double ,double ); 

void RK4_adaptativo( double (*funcion)(double, double,double), double , double , double , double ,double ); 

void euler_predictor_corrector( double (*funcion)(double, double, double), double , double , double , double , double ); 

#endif

