#include <stdio.h>
#include "ec_dif.h"



using namespace std;

vector<double> f(vector<double> r, vector<double> v, double t){
  vector<double> funcion={-r[0],0,0};     
  return funcion;
}

int main(){
  vector<double> posicion={5.,0,0};
  vector<double> velocidad={0,0,0};
  vector<double> tiempo={0,10,0.01};
  integrador_adaptativo(euler,f,posicion,velocidad,tiempo); 
    //double t0= 0.0;
    //double tf= 10;
    //double dt=0.01;
    //double x0=5.0;
    //double v0=0.0;
    //euler(f,t0,tf,dt,x0,v0);
    //cout<<a<<endl;
  vector<double> vec =f(posicion,velocidad,0.0);
  cout << vec[0] << "\n";
  
  return 0;
   
}

