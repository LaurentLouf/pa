#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include "lib.h"
#include <string>

using namespace std;

//---------------------------------Variables globales

manchoul CI;
louf CF;
WallsAngles wallsAngles;
string Jules;
bool vldt;
bool plot;
bool plotW;
int eval_f;

//---------------------------------Fonctions pour propagation

double Norm(double vx, double vy, double vz){return sqrt(pow(vx,2)+pow(vy,2)+pow(vz,2));}

double Re(double R, double rho, double v, double mu){return rho*v*2*R/mu;} // Reynolds number

double Cm(double R, double w, double v, double Cd){return 2*Cd*R*w/v;}// Magnus/Lift coefficient

double Ft(double vx, double vy, double vz){return -Cd*pi*R*R*rho*Norm(vx,vy,vz)/2;} // Drag Force

double Fm(double vx, double vy, double vz, double w){return  Cm(R,w,Norm(vx,vy,vz),Cd)*pi*R*R*rho*pow(Norm(vx,vy,vz),2)/2;} // Magnus Force

//---------------------------------Fonctions calcul intégrales

double Ints(double a, double b, double t1, double t2){return (1/(pow(a,2)+pow(b,2)))*((-b*exp(a*t2)*cos(b*t2)+a*exp(a*t2)*sin(b*t2))-(-b*exp(a*t1)*cos(b*t1)+a*exp(a*t1)*sin(b*t1)));}

double Intc(double a, double b, double t1, double t2){return (1/(pow(a,2)+pow(b,2)))*((a*exp(a*t2)*cos(b*t2)+b*exp(a*t2)*sin(b*t2))-(a*exp(a*t1)*cos(b*t1)+b*exp(a*t1)*sin(b*t1)));}

//---------------------------------Fonctions Runge-Kutta

double f1(double vx, double vy, double vz, double wx, double wy, double wz){
	if(Norm(wx,wy,wz)<0.1){
		return 1/m*(Ft(vx,vy,vz)*vx);
	}else{
		return 1/m*(Ft(vx,vy,vz)*vx+Fm(vx,vy,vz,Norm(wx,wy,wz))*(wy*vz-wz*vy)/Norm(wy*vz-wz*vy,wz*vx-wx*vz,wx*vy-wy*vx));
	}
}

double f2(double vx, double vy, double vz, double wx, double wy, double wz){
	if(Norm(wx,wy,wz)<0.1){
		return 1/m*(Ft(vx,vy,vz)*vy);
	}else{
		return 1/m*(Ft(vx,vy,vz)*vy+Fm(vx,vy,vz,Norm(wx,wy,wz))*(wz*vx-wx*vz)/Norm(wy*vz-wz*vy,wz*vx-wx*vz,wx*vy-wy*vx));
	}
}

double f3(double vx, double vy, double vz, double wx, double wy, double wz){
	if(Norm(wx,wy,wz)<0.1){
		return 1/m*(-m*g+Ft(vx,vy,vz)*vz);
	}else{
		return 1/m*(-m*g+Ft(vx,vy,vz)*vz+Fm(vx,vy,vz,Norm(wx,wy,wz))*(wx*vy-wy*vx)/Norm(wy*vz-wz*vy,wz*vx-wx*vz,wx*vy-wy*vx));
	}
}

//---------------------------------Fonction propagation dans l'air

void propag(){

if(plot==true){cout<< "\033[1;32m-------------------------------------------------- Propagation dans l'air\033[0m"<<endl;}

double x,y,z,wx,wy,wz,vx,vy,vz;

x=CI.xi;
y=CI.yi;
z=CI.zi;

wx=CI.wxi;
wy=CI.wyi;
wz=CI.wzi;

vx=CI.vxi;
vy=CI.vyi;
vz=CI.vzi;

double K11,K12,K13,K21,K22,K23,K31,K32,K33,K41,K42,K43;
int i=0;

if(plot==true){cout<<"Ité "<<i<<" vx="<<vx<<" vy="<<vy<<" vz="<<vz<<" x="<<x<<" y="<<y<<" z="<<z<<endl;}

if(Norm(wy*vz-wz*vy,wz*vx-wx*vz,wx*vy-wy*vx)<10e-6){if(plotW=true){cout<< "\033[1;31m Warning: norme produit vectoriel(w,v) très petite !\033[0m\n";}}

double epsilon=1e-3;

while(i<N && z>0 && y<16.71){ // stop si la balle=sol ou si la balle=mur

	K11=f1(vx,vy,vz,wx,wy,wz);
	K12=f2(vx,vy,vz,wx,wy,wz);
	K13=f3(vx,vy,vz,wx,wy,wz);

	K21=f1(vx+K11*dt/2,vy+K12*dt/2,vz+K13*dt/2,wx,wy,wz);
	K22=f2(vx+K11*dt/2,vy+K12*dt/2,vz+K13*dt/2,wx,wy,wz);
	K23=f3(vx+K11*dt/2,vy+K12*dt/2,vz+K13*dt/2,wx,wy,wz);

	K31=f1(vx+K21*dt/2,vy+K22*dt/2,vz+K23*dt/2,wx,wy,wz);
	K32=f2(vx+K21*dt/2,vy+K22*dt/2,vz+K23*dt/2,wx,wy,wz);
	K33=f3(vx+K21*dt/2,vy+K22*dt/2,vz+K23*dt/2,wx,wy,wz);

	K41=f1(vx+K31*dt,vy+K32*dt,vz+K33*dt,wx,wy,wz);
	K42=f2(vx+K31*dt,vy+K32*dt,vz+K33*dt,wx,wy,wz);
	K43=f3(vx+K31*dt,vy+K32*dt,vz+K33*dt,wx,wy,wz);

	vx=vx+dt/6*(K11+2*K21+2*K31+K41);
	vy=vy+dt/6*(K12+2*K22+2*K32+K42);
	vz=vz+dt/6*(K13+2*K23+2*K33+K43);

	x=x+dt*vx;
	y=y+dt*vy;
	z=z+dt*vz;

	if(abs(y-(11.89+1.71))<1e-2 && z<0.914+R){vldt=false; if(plotW==true){cout<<"\033[1;31m Warning: Trajectoire non valide: balle ne passe pas le filet \033[0m\n";} break;}

	if(vy<0 && y<(CF.yf)){Jules="cible"; break;}
	if(z<0){z=epsilon; break;}
	if(y>16.71){y=16.71-epsilon; break;}

	//cout<<"Ité "<<i<<" vx="<<vx<<" vy="<<vy<<" vz="<<vz<<" x="<<x<<" y="<<y<<" z="<<z<<endl;

	i++;

}

if(plot==true){cout<<"Ité "<<i<<" vx="<<vx<<" vy="<<vy<<" vz="<<vz<<" x="<<x<<" y="<<y<<" z="<<z<<endl;}

if(z==epsilon){Jules="sol";}
if(y==16.71-epsilon){Jules="mur";}

CI.xi=x;
CI.yi=y;
CI.zi=z;

CI.wxi=wx;
CI.wyi=wy;
CI.wzi=wz;

CI.vxi=vx;
CI.vyi=vy;
CI.vzi=vz;

}

//---------------------------------Fonctions Rebond

void Rebond_sol(){

if(plot==true){cout<< "\033[1;32m-------------------------------------------------- Rebond Sol\033[0m"<<endl;}

if(plot==true){cout<<"Ité "<<0<<" vx="<<CI.vxi<<" vy="<<CI.vyi<<" vz="<<CI.vzi<<" wx="<<CI.wxi<<" wy="<<CI.wyi<<" wz="<<CI.wzi<<endl;}

//----Projection dans le repère TN

int i0,i=0;

double theta,phi;

double vt[50];
double vn1,vn2;
double w[50];

vt[0]=sqrt(pow(CI.vxi,2)+pow(CI.vyi,2));
vn1=-CI.vzi;
vn2=-sCOR*vn1;
w[0]=copysign(1,CI.wxi)*sqrt(pow((CI.wxi),2)+pow((CI.wyi),2));

double M1,M2,I1,I2,Ky,Kt,Cy,Ct,Zy,Zt,Omy,Omt,Omdy,Omdt,Tc,t,dt,epsilon,a,b,c,gamma1,gamma2,alpha,beta,A,B,w0;
bool bastien=false;

epsilon=vt[0]*0.01;

M1=2/3.*m;
M2=1/3.*m;

I1=0.532*I;
I2=0.468*I;

Ky=435.45;
Kt=3.51;

Cy=0.0615;

Omy=sqrt(Ky/M1);
Omt=sqrt(Kt*(1/I1+1/I2));

Ct=2*sqrt(pow(Omt,2)-pow(Omy,2)+pow(Cy/(2*M1),2));

Zy=Cy/(2*M1*Omy);
Zt=Ct/(2*Omt);

Omdy=Omy*sqrt(1-pow(Zy,2));
Omdt=Omt*sqrt(1-pow(Zt,2));

Tc=pi/Omdy;

t=0;
dt=Tc/50;

for(i=0;i<50;i++){

	t=t+dt;

	if(abs(vt[i]-R*w[i])>epsilon){

	//cout<<"Ité "<<i<<" t/Tc="<<t/Tc<<" RAG - vt[i]="<<vt[i]<<" w[i]="<<w[i]<<endl;

	//----Phase Roulement avec Glissement

	a=1/(1/I1+1/I2);
	b=Ct;
	c=Kt;

	gamma1=-copysign(1,vt[0]-R*w[0])*sf*R*(Ky*vn1/Omdy-Cy*Zy*Omy*vn1/Omdy)/(I2*(1/I2+1/I1));
	gamma2=-copysign(1,vt[0]-R*w[0])*sf*R*Cy*vn1/(I2*(1/I2+1/I1));

	alpha=-Zy*Omy;
	beta=Omdy;

	A=(2*a*alpha*beta+b*beta)/(pow((a*(alpha*alpha-beta*beta)+b*alpha+c),2)+pow(2*a*alpha*beta+b*beta,2))*(gamma2+gamma1*(a*alpha*alpha-a*beta*beta+b*alpha+c)/(2*a*alpha*beta+b*beta));
	B=(A*(a*alpha*alpha-a*beta*beta+b*alpha+c)-gamma1)/(2*a*alpha*beta+b*beta);

	vt[i+1]=vt[0]-copysign(1,vt[0]-R*w[0])*(sf/m*(Cy*vn1/Omdy*exp(-Zy*Omy*t)*sin(Omdy*t)+Ky*vn1/(Omdy*(pow(Zy*Omy,2)+pow(Omdy,2)))*(-Zy*Omy*exp(-Zy*Omy*t)*sin(Omdy*t)-Omdy*exp(-Zy*Omy*t)*cos(Omdy*t)+Omdy)));

	w[i+1]=w[0]+(Ct/I2*(A*exp(-Zy*Omy*t)*sin(Omdy*t)+B*exp(-Zy*Omy*t)*cos(Omdy*t)-B)+Kt/I2*(A*Ints(-Zy*Omy,Omdy,0,t)+B*Intc(-Zy*Omy,Omdy,0,t))+gamma1*Ints(-Zy*Omy,Omdy,0,t)+gamma2*Intc(-Zy*Omy,Omdy,0,t))/8;

	}else{

	if(bastien==false){/*cout<<"Ité "<<i<<"--------------------------SOL: TRANSITION RAG->RSG-----------------------"<<endl;*/ bastien=true; i0=i;}

	//cout<<"Ité "<<i<<" t/Tc="<<t/Tc<<" RSG - vt[i]="<<vt[i]<<" w[i]="<<w[i]<<endl;

	//----Phase Roulement sans glissement

	a=1;
	b=Ct;
	c=Kt;

	gamma1=copysign(1,vt[0]-R*w[0])*sf*R/I2*(Ky*vn1/Omdy-Cy*Zy*Omy*vn1/Omdy)/(1/I2+1/I1);
	gamma2=copysign(1,vt[0]-R*w[0])*sf*R/I2*Cy*vn1/(1/I2+1/I1);

	alpha=-Zy*Omy;
	beta=Omdy;

	if(i0==0){
		w0=sqrt((m*(vt[0]*vt[0]+vn1*vn1)-I*w[0]*w[0]-I1*w[0]*w[0]-m*vn1)/(I2+m*R*R));
	}else{
		w0=w[i0];
	}

	w[i+1]=w[i0]+(Ct/(I2+m*R*R)*(-Zt*Omt*(w[i0]-w0)/Omdt*Ints(-Zt*Omt,Omdt,i0*dt,t)+(w[i0]-w0)*Intc(-Zt*Omt,Omdt,i0*dt,t)+alpha*A*Ints(alpha,beta,i0*dt,t)+beta*A*Intc(alpha,beta,i0*dt,t)+alpha*B*Intc(alpha,beta,i0*dt,t)-beta*B*Ints(alpha,beta,i0*dt,t))+Kt/(I2+m*R*R)*((w[i0]-w0)/Omdt*Ints(-Zt*Omt,Omdt,i0*dt,t)+A*Ints(alpha,beta,i0*dt,t)+B*Intc(alpha,beta,i0*dt,t)))/8;

	vt[i+1]=R*w[i+1];

	}

}

if(abs(CI.vxi)<1e-5){
	theta=pi/2;
}else{
	theta=abs(atan((CI.vyi)/(CI.vxi)));
}

if(CI.vxi<=0 && CI.vyi>0){
	theta=pi-theta;
}else if(CI.vxi>=0 && CI.vyi<0){
	theta=-theta;
}else if(CI.vxi<=0 && CI.vyi<0){
	theta=pi+theta;
}

if(abs(CI.wxi)<1e-5){
	phi=pi/2;
}else{
	phi=abs(atan((CI.wyi)/(CI.wxi)));
}

if(CI.wxi<=0 && CI.wyi>0){
	phi=pi-phi;
}else if(CI.wxi>=0 && CI.wyi<0){
	phi=-phi;
}else if(CI.wxi<=0 && CI.wyi<0){
	phi=pi+phi;
}

CI.vxi=vt[49]*cos(theta);
CI.vyi=vt[49]*sin(theta);
CI.vzi=-vn2;

CI.wxi=w[49]*cos(phi);
CI.wyi=w[49]*sin(phi);
CI.wzi=0;

if(plot==true){cout<<"Ité "<<i<<" vx="<<CI.vxi<<" vy="<<CI.vyi<<" vz="<<CI.vzi<<" wx="<<CI.wxi<<" wy="<<CI.wyi<<" wz="<<CI.wzi<<endl;}

}


void Rebond_mur(double alpha1, double beta1){

if(plot==true){cout<< "\033[1;32m-------------------------------------------------- Rebond Mur\033[0m"<<endl;}

if(plot==true){cout<<"Ité "<<0<<" vx="<<CI.vxi<<" vy="<<CI.vyi<<" vz="<<CI.vzi<<" wx="<<CI.wxi<<" wy="<<CI.wyi<<" wz="<<CI.wzi<<endl;}

//----Projection dans le repère TN

double theta,phi,vx1,vy1,vz1,wx1,wy1,wz1;

vx1=copysign(1,CI.vxi)*Norm(cos(alpha1)*(CI.vxi),sin(alpha1)*(CI.vyi),0);
vy1=copysign(1,CI.vyi)*Norm(-sin(alpha1)*cos(beta1)*(CI.vxi),cos(alpha1)*cos(beta1)*(CI.vyi),sin(beta1)*(CI.vzi));
vz1=copysign(1,CI.vzi)*Norm(sin(alpha1)*sin(beta1)*(CI.vxi),-cos(alpha1)*sin(beta1)*(CI.vyi),cos(beta1)*(CI.vzi));

wx1=copysign(1,CI.wxi)*Norm(cos(alpha1)*(CI.wxi),sin(alpha1)*(CI.wyi),0);
wy1=copysign(1,CI.wyi)*Norm(-sin(alpha1)*cos(beta1)*(CI.wxi),cos(alpha1)*cos(beta1)*(CI.wyi),sin(beta1)*(CI.wzi));
wz1=copysign(1,CI.wzi)*Norm(sin(alpha1)*sin(beta1)*(CI.wxi),-cos(alpha1)*sin(beta1)*(CI.wyi),cos(beta1)*(CI.wzi));

int i0,i=0;

double vt[50];
double vn1,vn2;
double w[50];

vt[0]=Norm(vx1,0,vz1);
vn1=vy1;
vn2=-mCOR*vn1;
w[0]=copysign(1,CI.wxi)*Norm(wx1,0,wz1);

double M1,M2,I1,I2,Ky,Kt,Cy,Ct,Zy,Zt,Omy,Omt,Omdy,Omdt,Tc,t,dt,epsilon,a,b,c,gamma1,gamma2,alpha,beta,A,B,w0;
bool bastien=false;

epsilon=vt[0]*0.01;

M1=2/3.*m;
M2=1/3.*m;

I2=0.6*I;
I1=0.4*I;

Ky=435.45;
Kt=3.51;

Kt=Kt/10;

Cy=0.0615;

Omy=sqrt(Ky/M1);
Omt=sqrt(Kt*(1/I1+1/I2));

Ct=2*sqrt(pow(Omt,2)-pow(Omy,2)+pow(Cy/(2*M1),2));

Zy=Cy/(2*M1*Omy);
Zt=Ct/(2*Omt);

Omdy=Omy*sqrt(1-pow(Zy,2));
Omdt=Omt*sqrt(1-pow(Zt,2));

Tc=pi/Omdy;

t=0;
dt=Tc/50;

for(i=0;i<50;i++){

	t=t+dt;

	if(abs(vt[i]-R*w[i]/2)>epsilon){

	//cout<<"Ité "<<i<<" t/Tc="<<t/Tc<<" RAG - vt[i]="<<vt[i]<<" w[i]="<<w[i]<<endl;

	//----Phase Roulement avec Glissement

	a=1/(1/I1+1/I2);
	b=Ct;
	c=Kt;

	gamma1=-copysign(1,vt[0]-R*w[0]/2)*mf*R*(Ky*vn1/Omdy-Cy*Zy*Omy*vn1/Omdy)/(I2*(1/I2+1/I1));
	gamma2=-copysign(1,vt[0]-R*w[0]/2)*mf*R*Cy*vn1/(I2*(1/I2+1/I1));

	alpha=-Zy*Omy;
	beta=Omdy;

	A=(2*a*alpha*beta+b*beta)/(pow((a*(alpha*alpha-beta*beta)+b*alpha+c),2)+pow(2*a*alpha*beta+b*beta,2))*(gamma2+gamma1*(a*alpha*alpha-a*beta*beta+b*alpha+c)/(2*a*alpha*beta+b*beta));
	B=(A*(a*alpha*alpha-a*beta*beta+b*alpha+c)-gamma1)/(2*a*alpha*beta+b*beta);

	vt[i+1]=vt[0]-copysign(1,vt[0]-R*w[0]/2)*(mf/m*(Cy*vn1/Omdy*exp(-Zy*Omy*t)*sin(Omdy*t)+Ky*vn1/(Omdy*(pow(Zy*Omy,2)+pow(Omdy,2)))*(-Zy*Omy*exp(-Zy*Omy*t)*sin(Omdy*t)-Omdy*exp(-Zy*Omy*t)*cos(Omdy*t)+Omdy)));

	w[i+1]=w[0]+(Ct/I2*(A*exp(-Zy*Omy*t)*sin(Omdy*t)+B*exp(-Zy*Omy*t)*cos(Omdy*t)-B)+Kt/I2*(A*Ints(-Zy*Omy,Omdy,0,t)+B*Intc(-Zy*Omy,Omdy,0,t))+gamma1*Ints(-Zy*Omy,Omdy,0,t)+gamma2*Intc(-Zy*Omy,Omdy,0,t))/8-m*g*R*t/I2;

	}else{

	if(bastien==false){/*cout<<"Ité "<<i<<"--------------------------MUR: TRANSITION RAG->RSG-----------------------"<<endl;*/ i0=i; bastien=true;}

	//cout<<"Ité "<<i<<" t/Tc="<<t/Tc<<" RSG - vt[i]="<<vt[i]<<" w[i]="<<w[i]<<endl;

	//----Phase Roulement sans glissement

	a=1;
	b=Ct;
	c=Kt;

	gamma1=-copysign(1,vt[0]-R*w[0]/2)*mf*R/I2*(Ky*vn1/Omdy-Cy*Zy*Omy*vn1/Omdy)/(I2*(1/I2+1/I1));
	gamma2=-copysign(1,vt[0]-R*w[0]/2)*mf*R/I2*Cy*vn1/(I2*(1/I2+1/I1));

	alpha=-Zy*Omy;
	beta=Omdy;

	if(i0==0){
		w0=sqrt((m*(vt[0]*vt[0]+vn1*vn1)-I*w[0]*w[0]-I1*w[0]*w[0]-m*vn1)/(I2+m*R*R));
	}else{
		w0=w[i0];
	}

	w[i+1]=w[i0]+(Ct/(I2+m*R*R)*(-Zt*Omt*(w[i0]-w0)/Omdt*Ints(-Zt*Omt,Omdt,i0*dt,t)+(w[i0]-w0)*Intc(-Zt*Omt,Omdt,i0*dt,t)+alpha*A*Ints(alpha,beta,i0*dt,t)+beta*A*Intc(alpha,beta,i0*dt,t)+alpha*B*Intc(alpha,beta,i0*dt,t)-beta*B*Ints(alpha,beta,i0*dt,t))+Kt/(I2+m*R*R)*((w[i0]-w0)/Omdt*Ints(-Zt*Omt,Omdt,i0*dt,t)+A*Ints(alpha,beta,i0*dt,t)+B*Intc(alpha,beta,i0*dt,t)))/8-m*g*R*t/I2;

	vt[i+1]=R*w[i+1]/2;

	}

}

if(abs(CI.vxi)<1e-5){
	theta=pi/2;
}else{
	theta=abs(atan((CI.vyi)/(CI.vxi)));
}

if(CI.vxi<=0 && CI.vyi>0){
	theta=pi-theta;
}else if(CI.vxi>=0 && CI.vyi<0){
	theta=-theta;
}else if(CI.vxi<=0 && CI.vyi<0){
	theta=pi+theta;
}

if(abs(CI.wxi)<1e-5){
	phi=pi/2;
}else{
	phi=abs(atan((CI.wyi)/(CI.wxi)));
}

if(CI.wxi<0 && CI.wyi>0){
	phi=pi-phi;
}else if(CI.wxi>0 && CI.wyi<0){
	phi=-phi;
}else if(CI.wxi<0 && CI.wyi<0){
	phi=pi+phi;
}

theta=-theta;
phi=-phi;

if(vx1>0){
	if(wx1>0){
		if(wz1>0){
			phi=pi-phi;
		}
	}else{
		if(wz1>0){
			phi=pi-phi;
		}else{
			phi=pi-phi;
			if(bastien==false){
				theta=-theta;
			}
		}
	}
}


if(vx1<0){
	if(wx1>0){
		if(wz1>0){
			theta=-theta;
			phi=pi-phi;
		}else{
			phi=pi-phi;
		}
	}else{
		if(wz1>0){
			if(bastien==true){
				theta=-theta;
			}
		}else{
			theta=-theta;
		}
	}
}


vx1=vt[49]*sin(theta);
vy1=vn2;
vz1=vt[49]*cos(theta);

CI.vxi=cos(alpha1)*vx1-sin(alpha1)*cos(beta1)*vy1+sin(alpha1)*sin(beta1)*vz1;
CI.vyi=sin(alpha1)*vx1+cos(alpha1)*cos(beta1)*vy1-cos(alpha1)*sin(beta1)*vz1;
CI.vzi=sin(beta1)*vy1+cos(beta1)*vz1;

wx1=w[49]*sin(phi);
wy1=0;
wz1=w[49]*cos(phi);

CI.wxi=cos(alpha1)*wx1-sin(alpha1)*cos(beta1)*wy1+sin(alpha1)*sin(beta1)*wz1;
CI.wyi=sin(alpha1)*wx1+cos(alpha1)*cos(beta1)*wy1-cos(alpha1)*sin(beta1)*wz1;
CI.wzi=sin(beta1)*wy1+cos(beta1)*wz1;

if(plot==true){cout<<"Ité "<<i<<" vx="<<CI.vxi<<" vy="<<CI.vyi<<" vz="<<CI.vzi<<" wx="<<CI.wxi<<" wy="<<CI.wyi<<" wz="<<CI.wzi<<endl;}

}

double fitfun(double alpha, double beta){

eval_f++;

propag();
if(vldt==false){return 1000+rand()/1000.;}

if(Jules.compare("sol")==0){
	Rebond_sol();
	if(CI.yi<1.71+11.89){if(plotW==true){cout<<"\033[1;31m Warning: Trajectoire non valide: Premier rebond sur le sol avant filet \033[0m\n";} vldt=false; return 1000+rand()/1000.;}
	propag();
	if(vldt==false){return 1000+rand()/1000.;}
	if(Jules.compare("mur")!=0){if(plotW==true){cout<<"\033[1;31m Warning: Trajectoire non valide: Deuxieme rebond sur le sol \033[0m\n";} vldt=false; return 1000+rand()/1000.;}
	Rebond_mur(alpha,beta);
}else if(Jules.compare("mur")==0){
	Rebond_mur(alpha,beta);
}

propag();
if(vldt==false){return 1000+rand()/1000.;}

if(Jules.compare("cible")==0){return Norm((CI.xi-CF.xf),(CI.yi-CF.yf),(CI.zi-CF.zf));}
Rebond_sol();
propag();
if(vldt==false){return 1000+rand()/1000.;}

return Norm((CI.xi-CF.xf),(CI.yi-CF.yf),(CI.zi-CF.zf));

}

void* Pb_Inv(void* arg){

manchoul tmp1=CI;

double l=10.97;
double ym=1.71+11.89+6.40;

double M,Mnew;
double A1,A2;
double B1,B2;

Mnew=10000;

//--------------------------Recherche linéaire restreinte

double psi1,psi2,psib,psibb,psibbb,psit;

CI=tmp1;

if(CI.xi==l/2){
	psibbb=pi/2;
}else{
	psibbb=atan((ym-CI.yi)/(l/2-CI.xi));
}

if(CI.xi==-l/2){
	psit=pi/2;
}else{
	psit=atan((ym-CI.yi)/(-l/2-CI.xi))+pi;
}

double x3,x4,y3,y4;

if(CF.xf>=0){
	x4=0;
	x3=l/2;
}else{
	x4=-l/2;
	x3=0;
	}

y3=CF.yf;
y4=CF.yf;

psib=atan((ym-y4)/(l/2-x4));

psibb=atan((ym-y3)/(-l/2-x3))+pi;

psi1=(pi-psib-psibbb)/2;

psi2=(psit+psibb-pi)/2;

int N=20;

for(int i=0;i<N;i++){
	for(int j=0;j<N;j++){
		CI=tmp1;
      		vldt=true;
		Jules="joueur";
		A1=-psi1+(psi2+psi1)*(i/(double)(N-1));
		A2=-pi/6+(pi/6+pi/20)*(j/(double)(N-1));
		M=fitfun(A1,A2);
		if(M<Mnew){
			Mnew=M;
			B1=A1;
			B2=A2;
		}
		if(Mnew<1e-1){break;}
	}
}

cout<<"Mfinal="<<Mnew<<endl;

/**alpha=B1;
*beta=B2;*/
wallsAngles.wall1_alpha = B1;
wallsAngles.wall1_beta = B2;

}