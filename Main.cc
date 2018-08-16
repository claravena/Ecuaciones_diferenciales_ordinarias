#include <stdio.h>
#include <string.h>
#include "ec_dif.h"



using namespace std;

double f(double x, double v, double t){
  return -x-v;
}



int main(){
  double t0= 0.0;
  double tf= 10;
  double dt=0.01;
  double x0=5.0;
  double v0=0.0;
  //char a = euler; 
  RK4_adaptativo(f,t0,tf,dt,x0,v0);
  //cout<<a<<endl;

  return 0;
   
}

