#include "ec_dif.h"


void integrador(vector<double> (*metodo)(double a, double pos,double vel,double t_1,double dt), vector<double> (*funcion)(vector<double> r_pos, vector<double> v_pos, double tiempos), vector<double> r ,vector<double> v, vector<double> tiempo ){
  double t= tiempo[0];
  double tf=tiempo[1]; 
  double h= tiempo[2]; 
  vector<double> vec_r=r;
  vector<double> vec_v=v;
  
  while(t<=tf){
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[0]<<" "<<vec_v[1]<<" "<<vec_r[2]<<endl;
    cout <<t<<" "<< vec_v[0]<<" "<<vec_v[1]<<" "<<vec_r[2]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_r[0]<<" "<<vec_r[1]<<" "<<vec_r[2]<<endl;
    cout <<t<<" "<< vec_r[0]<<" "<<vec_r[1]<<" "<<vec_r[2]<<endl;
    archivo1.close();
    vector<double> a = funcion(vec_r,vec_v,t);
    vector<double> vec_aux_r=vec_r; //copian el paso anterior de posicion;
    vector<double> vec_aux_v=vec_v; //copian el paso anterio de posicion;
    for(int j=0; j<3; ++j){
      vector<double> vec_estado=metodo(a[j],vec_r[j],vec_v[j],t,h);
      vec_r[j]=vec_estado[0];
      vec_v[j]=vec_estado[1];
    }
    t+= h;
    
  }
}

void integrador_adaptativo(vector<double> (*metodo)(double a, double pos,double vel,double t_1,double dt), vector<double> (*funcion)(vector<double> r_pos, vector<double> v_pos, double tiempos), vector<double> r ,vector<double> v, vector<double> tiempo ){
  double t= tiempo[0];
  double tf=tiempo[1]; 
  double h= tiempo[2]; 
  vector<double> vec_r=r;
  vector<double> vec_v=v;
  double delta_esperado=0.00000001;
  double S1=0.9;
  double S2=4.0; 
 

  while(t<=tf){
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[0]<<" "<<vec_v[1]<<" "<<vec_r[2]<<endl;
    cout <<t<<" "<< vec_v[0]<<" "<<vec_v[1]<<" "<<vec_r[2]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_r[0]<<" "<<vec_r[1]<<" "<<vec_r[2]<<endl;
    cout <<t<<" "<< vec_r[0]<<" "<<vec_r[1]<<" "<<vec_r[2]<<endl;
    archivo1.close();
    vector<double> a = funcion(vec_r,vec_v,t);
    vector<double> vec_aux_r=vec_r; //copian el paso anterior de posicion;
    vector<double> vec_aux_v=vec_v; //copian el paso anterio de posicion;
    for(int j=0; j<3; ++j){
      vector<double> vec_estado=metodo(a[j],vec_r[j],vec_v[j],t,h);
      vec_r[j]=vec_estado[0];
      vec_v[j]=vec_estado[1];
    }
    t+= h;
    double norma_actual=sqrt(vec_r[0]*vec_r[0]+vec_r[1]*vec_r[1]+vec_r[2]*vec_r[2]);
    double norma_anterior=sqrt(vec_aux_r[0]*vec_aux_r[0]+vec_aux_r[1]*vec_aux_r[1]+vec_aux_r[2]*vec_aux_r[2]);
    double delta_estimado=abs(norma_actual-norma_anterior);
    double t_estimado= pow(abs(delta_esperado/delta_estimado),1./5.0)*h;
    if ((S1*t_estimado) > (S2*h)){
      h= S2*tiempo[2];
    }
    else if ((S1*t_estimado)< (h/S2)){
      h= tiempo[2]/2.;
    
    }
    else {
      h=S1*t_estimado; 
    }
  }
}








vector<double> euler( double a, double x, double v,double t,double dt){
  vector<double> vec_estado={x,v,t};
  double h=dt;
  vec_estado[1]=(v+h*a);
  vec_estado[0]=(x+h*v);
  return vec_estado; 
}




void euler( double (*funcion)(double, double, double),double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
 
  for(int i=0; i<((tf-t0)/dt); ++i){
    double a= funcion(vec_x[i], vec_v[i], t);
    vec_v.push_back(vec_v[i]+h*a);
    vec_x.push_back(vec_x[i]+h*vec_v[i]);
    ofstream archivo("datos_velocidad_euler.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion_euler.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}

void euler_cromer( double (*funcion)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  
  //ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double a= funcion(vec_x[i],vec_v[i], t);
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


void punto_medio( double (*funcion)(double, double,double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};

  for(int i=0; i<((tf-t0)/dt); ++i){
    double a= funcion(vec_x[i], vec_v[i], t);
    vec_v.push_back(vec_v[i]+h*a);
    vec_x.push_back(vec_x[i]+0.5*h*(vec_v[i+1]+vec_v[i]));
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}  

void verlet( double (*funcion)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  
   for(int i=0; i<((tf-t0)/dt); ++i){
    double a_0=funcion(vec_x[i], vec_v[i], t);
    vec_x.push_back(vec_x[i]+vec_v[i]*h+0.5*a_0*h*h);
    double a_1= funcion(vec_x[i+1], vec_v[i+1],(t+h));
    vec_v.push_back(vec_v[i]+0.5*h*(a_0+a_1));
 
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
}

void leap_frog( double (*funcion)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  for(int i=0; i<((tf-t0)/dt); ++i){
    double v_media=(vec_v[i]+funcion(vec_x[i], vec_v[i], t)*0.5*h);
    vec_x.push_back(vec_x[i]+v_media*h);
    vec_v.push_back(v_media+0.5*h*funcion(vec_x[i+1],vec_v[i+1], (t+h)));
    ofstream archivo("datos_velocidad.dat", ios::app);
    archivo <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_posicion.dat", ios::app);
    archivo1 <<t<<" "<< vec_x[i]<<endl;
    archivo1.close();
    t+= h;
  }
} 

void RK2( double (*funcion)(double, double,double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double x= x0;
  double v= v0;
  double t=t0;
  ofstream archivo("datos.dat");
  for(int i=0; i<((tf-t0)/dt); ++i){
    double k1= h*v;
    double l1= h*funcion(x,v,t);
    double k2= h*v*l1;
    double l2= h*funcion(x+k1,v+l1,t+h);
    x+= 0.5*(k1+k2);
    v+= 0.5*(l1+l2);
    t+= h;
    ofstream archivo("datos_posiciones_RK2.dat", ios::app);
    archivo <<t<<" "<< x<<endl;
    archivo.close();
    ofstream archivo1("datos_velocidad_RK2.dat", ios::app);
    archivo1 <<t<<" "<< v<<endl;
    archivo1.close();
  }
}  
 
void RK4( double (*funcion)(double, double,double), double t0, double tf, double dt, double x0,double v0){
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
    double l1= h*funcion(x,v,t);
    double k2= h*(v+0.5*l1);
    double l2= h*funcion(x+0.5*k1, v+0.5*l1,t+0.5*h);
    double k3 = h*(v+0.5*l2);
    double l3= h*funcion(x+0.5*k2, v+0.5*l2,t+0.5*h);
    double k4=h*(v+l3);
    double l4= h*funcion(x+k3, v+l3,t+h);
    x+= (k1+2*k2+2*k3+k4)/6.0;
    v+=(l1+2*l2+2*l3+l4)/6;
    t+= h;
  }
}

void RK4_adaptativo( double (*funcion)(double, double,double), double t0, double tf, double dt, double x0,double v0){
  double h= dt;
  double x= x0;
  double v= v0; 
  double t=t0;
  vector<double> vec_x={x0};
  vector<double> vec_v={v0};
  
  for(int i=0; i<10000000; ++i){
    ofstream archivo("datos_posicion_RK4a..dat", ios::app);
    archivo <<t<<" "<< x<<endl;
    archivo.close();
    ofstream archivo1("datos_velocidad_RK4a.dat", ios::app);
    archivo1 <<t<<" "<< v<<endl;
    archivo.close();
    double k1= h*v; 
    double l1= h*funcion(x,v,t);
    double k2= h*(v+0.5*l1);
    double l2= h*funcion(x+0.5*k1, v+0.5*l1,t+0.5*h);
    double k3 = h*(v+0.5*l2);
    double l3= h*funcion(x+0.5*k2, v+0.5*l2,t+0.5*h);
    double k4=h*(v+l3);
    double l4= h*funcion(x+k3, v+l3,t+h);
    x+= (k1+2*k2+2*k3+k4)/6.0;
    vec_x.push_back(x);
    v+=(l1+2*l2+2*l3+l4)/6;
    vec_v.push_back(v);
    t+=h;
    //cout << t<<endl;
    //Paso adaptativo. 
    double delta0=0.00000001;
    double delta_estimado=abs(vec_x[i+1]-vec_x[i]);
    double t_estimado= pow(abs(delta0/delta_estimado),1./5.0)*h;
    double S1=0.9;
    double S2=4.0; 
    if ((S1*t_estimado) > (S2*h)){
      h= S2*dt;
    }
    else if ((S1*t_estimado)< (h/S2)){
      h= dt/2.;

    }
    else {
      h=S1*t_estimado; 
    }
  }
}

void euler_predictor_corrector( double (*funcion)(double, double, double), double t0, double tf, double dt, double x0, double v0){
  double h= dt;
  double t=t0;
  vector<double> vec_v={v0};
  vector<double> vec_x={x0};
  vec_v.push_back(v0+h*funcion(x0, v0, t));
  vec_x.push_back(x0+h*vec_v[0]);
  cout << vec_v[0]<< " "<<vec_x[0]<<endl; 
  for(int i=1; i<((tf-t0)/dt); ++i){
    t+= h;
    double a= funcion(vec_x[i], vec_v[i], t);
    double a1= funcion(vec_x[i-1], vec_v[i-1], (t+h)); 
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

void PEFRL( double (*f)(double, double,double), double t0, double tf, double dt, double x0,double v0){
  double h= dt;
  double x= x0;
  double v= v0; 
  double t=t0+4*h;
  double m=1;
  double eta= 0.1786178958448091;
  double lambda=0.2123418310626054;
  double xi=-0.6626458266981849;
  vector<double> vec_x={x0};
  vector<double> vec_v={v0};
  vec_x.push_back(vec_x[0]+eta*h*vec_v[0]);
  vec_v.push_back(vec_v[0]+(1.0-2.0*lambda)*h*f(vec_x[1],vec_v[1],t+h)/(2.*m));
  //cout<<vec_x[1]<<" "<<vec_v[1]<<endl;
  vec_x.push_back(vec_x[1]+xi*h*vec_v[1]);
  vec_v.push_back(vec_v[1]+lambda*h*f(vec_x[2],vec_v[2],t+2*h)/m);
  //cout<<vec_x[2]<<" "<<vec_v[2]<<endl;
  vec_x.push_back(vec_x[2]+(1-2*(xi+eta))*h*vec_v[2]);
  vec_v.push_back(vec_v[2]+lambda*h*f(vec_x[3],vec_v[3],t+3*h)/m);
  //cout<<vec_x[3]<<" "<<vec_v[3]<<endl;
  vec_x.push_back(vec_x[3]+xi*h*vec_v[3]);
  for(int i=4; i<((tf-t0)/dt); ++i){
    vec_v.push_back(vec_v[i-1]+(1-2.*lambda)*h*f(vec_x[i],vec_v[i],t)/(2.*m));
    cout<<vec_x[i-1]<<" "<<vec_v[i-1]<<endl;
    vec_x.push_back(vec_x[i]+eta*h*vec_v[i]);
    t+= h;
    ofstream archivo("datos_posicion_RK4.dat", ios::app);
    archivo <<t<<" "<< vec_x[i]<<endl;
    archivo.close();
    ofstream archivo1("datos_velocidad_RK4.dat", ios::app);
    archivo1 <<t<<" "<< vec_v[i]<<endl;
    archivo.close();
  }
}





