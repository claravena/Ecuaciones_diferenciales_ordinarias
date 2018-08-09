//------------------ Archivo  ec_dif.h ------------------------------------
#ifndef EC_DIF_H
#define EC_DIF_H
#include <iostream>
#include <vector>
#include <math.h>
#include "cmath"
#include <fstream>


using namespace std;

void euler( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0);

void euler_cromer( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0);

void punto_medio( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0, double v0);

void verlet( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0);

void leap_frog( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0);

void RK2( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0, double v0);

void RK4( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0,double v0); 

void RK4_adaptativo( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0,double v0); 

#endif


void euler( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  
  
  //ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double a= parametro(vec_x[i], vec_v[i], t);
    vec_v.push_back(vec_v[i]+h*a);
    vec_x.push_back(vec_x[i]+h*vec_v[i]);

    //cout<<x<<" "<<vec_y[i]<<endl;
    ofstream archivo("datos_velocidad_euler.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion_euler.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}  

void euler_cromer( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  
  //ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double a= parametro(vec_x[i],vec_v[i], t);
    vec_v.push_back(vec_v[i]+h*a);
    vec_x.push_back(vec_x[i]+h*vec_v[i+1]);

    //cout<<x<<" "<<vec_y[i]<<endl;
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}  


void punto_medio( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};

  //ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double a= parametro(vec_x[i], vec_v[i], t);
    vec_v.push_back(vec_v[i]+h*a);
    vec_x.push_back(vec_x[i]+0.5*h*(vec_v[i+1]+vec_v[i]));

    //cout<<x<<" "<<vec_y[i]<<endl;
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}  

void verlet( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  
  //ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double a_0= parametro(vec_x[i], vec_v[i], t);
    vec_x.push_back(vec_x[i]+vec_v[i]*h+0.5*a_0*h*h);
    double a_1= parametro(vec_x[i+1], vec_v[i+1],(t+h));
    vec_v.push_back(vec_v[i]+0.5*h*(a_0+a_1));
    //cout << vec_v[i]<<" "<<h<<endl; 

    //cout<<x<<" "<<vec_y[i]<<endl;
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}

//no funciona para sistemas con primera derivada. 
void leap_frog( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};

  //inicializamos el primer término con euler.
  //vec_v.push_back(vec_v[0]+h*parametro(vec_x[0], vec_v[0], t));
  //vec_x.push_back(vec_x[0]+h*vec_v[0]);
  
   //ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double v_media=(vec_v[i]+parametro(vec_x[i], vec_v[i], t)*0.5*h);
    //cout <<parametro(vec_x[i], vec_v[i], t)<<" "<<vec_x[i]<<" "<<vec_v[i]<<" "<<t<<endl; 



    vec_x.push_back(vec_x[i]+v_media*h);
    vec_v.push_back(v_media+0.5*h*parametro(vec_x[i+1],vec_v[i+1], (t+h))); 
    //cout<<x<<" "<<vec_y[i]<<endl;
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
} 

void RK2( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double x= x0;
  double v= v0;
  double t=t0;
  ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double k1= h*v;
    double l1= h*parametro(x,v,t);
    double k2= h*v*l1;
    double l2= h*parametro(x+k1,v+l1,t+h);
    x+= 0.5*(k1+k2);
    v+= 0.5*(l1+l2);
    t+= h;
    //cout<<x<<" "<<y<<endl;
    ofstream archivo("datos_posiciones_RK2.dat", ios::app);
    archivo <<t<<" "<< x<<endl;
    archivo.close();
    ofstream archivo1("datos_velocidad_RK2.dat", ios::app);
    archivo1 <<t<<" "<< v<<endl;
    archivo1.close();
  }
}  


//RK4 está perfecto. 
void RK4( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0,double v0){
  double h= dt;
  double x= x0;
  double v= v0; 
  double t=t0;

  for(int i=0; i<((tf-t0)/dt); ++i){
    ofstream archivo("datos_posicion_RK4.dat", ios::app);
    archivo <<t<<" "<< x<<endl;
    archivo.close();
    ofstream archivo1("datos_velocidad_RK4.dat", ios::app);
    archivo1 <<t<<" "<< v<<endl;
    archivo.close();
    double k1= h*v; 
    double l1= h*parametro(x,v,t);
    double k2= h*(v+0.5*l1);
    double l2= h*parametro(x+0.5*k1, v+0.5*l1,t+0.5*h);
    double k3 = h*(v+0.5*l2);
    double l3= h*parametro(x+0.5*k2, v+0.5*l2,t+0.5*h);
    double k4=h*(v+l3);
    double l4= h*parametro(x+k3, v+l3,t+h);
    x+= (k1+2*k2+2*k3+k4)/6.0;
    v+=(l1+2*l2+2*l3+l4)/6;
    t+= h;
    //cout<<t<<" "<<x<<endl;
  }
}

void RK4_adaptativo( double (*parametro)(double, double,double), double t0, double tf, double dt, double x0,double v0){
  double h= dt;
  double x= x0;
  double v= v0; 
  double t=t0;
  vector<double> vec_x={x0};
  vector<double> vec_v={v0};
  
  for(int i=0; i<((tf-t0)/dt); ++i){
    ofstream archivo("datos_posicion_RK4.dat", ios::app);
    archivo <<t<<" "<< x<<endl;
    archivo.close();
    ofstream archivo1("datos_velocidad_RK4.dat", ios::app);
    archivo1 <<t<<" "<< v<<endl;
    archivo.close();
    double k1= h*v; 
    double l1= h*parametro(x,v,t);
    double k2= h*(v+0.5*l1);
    double l2= h*parametro(x+0.5*k1, v+0.5*l1,t+0.5*h);
    double k3 = h*(v+0.5*l2);
    double l3= h*parametro(x+0.5*k2, v+0.5*l2,t+0.5*h);
    double k4=h*(v+l3);
    double l4= h*parametro(x+k3, v+l3,t+h);
    x+= (k1+2*k2+2*k3+k4)/6.0;
    vec_x.push_back(x);
    v+=(l1+2*l2+2*l3+l4)/6;
    vec_v.push_back(v);
    cout << t<<endl;
    //Paso adaptativo. 
    double delta0=0.001;
    double delta_estimado=abs(vec_x[i+1]-vec_x[i]);
    double t_estimado= pow(abs(delta0/delta_estimado),1./5.0)*h;
    double S1=0.9;
    double S2=4.0; 
    if ((S1*t_estimado) > (S2*h)){
      t+= S2*h;
    }
    else if ((S1*t_estimado)< (h/S2)){
      t+= h/S2;

    }
    else {
      t+=S1*t_estimado; 
      
    }
  }
}

void euler_predictor_corrector( double (*parametro)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  vec_v.push_back(v0+h*parametro(x0, v0, t));
  vec_x.push_back(x0+h*vec_v[0]);
  cout << vec_v[0]<< " "<<vec_x[0]<<endl; 
  for(int i=1; i<((tf-t0)/dt); ++i){
    t+= h;
    double a= parametro(vec_x[i], vec_v[i], t);
    double a1= parametro(vec_x[i-1], vec_v[i-1], (t+h)); 
    vec_v.push_back(vec_v[i]+0.5*h*(a+a1));
    vec_x.push_back(vec_x[i]+0.5*h*(vec_v[i]+vec_v[i-1]));
    //cout<<t<<" "<<vec_x[i]<<endl;
    ofstream archivo("datos_velocidad_euler_pc.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion_euler_pc.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
  }
}  



