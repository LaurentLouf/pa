#ifndef LIB_H
#define LIB_H
 
using namespace std;

//---------------------------------Constantes

const double rho=1.204;  // 20°C
const double R=3.2e-2;   // 3.175<=R<=3.334cm dans la réglementation
const double m=57e-3;    // 56.7<=m<=58.5g dans la réglementation
const double mu=1.85e-5; // Viscosité dynamique de l'air à 20°C
const double Cd=0.55;    // Drag coefficient, à faire varier entre 0.55 et 0.61
const double pi=3.141592653589793;
const double g=9.81;     // Pesanteur
const double dt=1e-3;	 // Time step
const double N=3000;     // Iteration number
const double I=3.1636e-5; // Moment d'inertie

const double sCOR=0.76; // Coefficient de restitution du sol 
const double sf=0.55; // Coefficient de frottement du sol

const double mCOR=0.9; // Coefficient de restitution du mur
const double mf=0.6; // Coefficient de frottement du mur

//---------------------------------Structures globales de fonctionnement

extern bool vldt;
extern bool plot;
extern bool plotW;
extern string Jules;

typedef struct manchoul{

double xi;
double yi;
double zi;
double vxi;
double vyi;
double vzi;
double wxi;
double wyi;
double wzi;

} manchoul;

typedef struct louf{

double xf;
double yf;
double zf;

} louf;

typedef struct{
	double wall1_alpha, wall1_beta;
	double wall2_alpha, wall2_beta;
	double wall3_alpha, wall3_beta;
	double wall4_alpha, wall4_beta;
} WallsAngles;

extern manchoul CI;
extern louf CF;
extern WallsAngles wallsAngles;

//---------------------------------Fonctions pour propagation

double Norm(double, double, double);

double Re(double, double, double, double);

double Cm(double, double, double, double);

double Ft(double, double, double);

double Fm(double, double, double, double);

//---------------------------------Fonctions calcul intégrales

double Intc(double, double, double, double);

double Ints(double, double, double, double);

//---------------------------------Fonctions Runge-Kutta

double f1(double, double, double, double, double, double);

double f2(double, double, double, double, double, double);

double f3(double, double, double, double, double, double);

//---------------------------------Fonction propagation dans l'air
 
void propag();

//---------------------------------Fonctions Rebond

void Rebond_sol();

void Rebond_mur(double, double);

//---------------------------------Probleme inverse

extern int eval_f;

double fitfun(double, double);

void* Pb_Inv(void*);

#endif
