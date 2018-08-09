#include <stdio.h>
#include <string.h>
#include "ec_dif.h"




using namespace std;

double f(double x, double v, double t){
  double M=5.98*pow(10,24);
  double G= 6.67*pow(10,-11);
  double d=6370000;
    
  return -x-v;
}


int main(){
  double t0= 0.0;
  double tf= 10;
  double dt=0.01;
  double x0=5.0;
  double v0=0.0;
  //char a = euler; 
 euler_predictor_corrector(f,t0,tf,dt,x0,v0);
  //cout<<a<<endl;

  return 0;
   
}

//double (*parametro)(double double),
