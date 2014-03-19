#include "postDetection.h"

using namespace std;

void getBallPositionFromCamera(int xPixels, int yPixels, double diameterPixels, double* x, double* y, double* z)
{
	double alpha = abs(diameterPixels * CAMERA_PIX2RAD);
	double gamma = abs(CAMERA_PIX2RAD * sqrt(xPixels*xPixels + yPixels*yPixels));
	*z = - BALL_DIAMETER / (2*( tan(alpha/2 + gamma) - tan(gamma)));
	*x = *z * tan(CAMERA_PIX2RAD * (xPixels - CAMERA_X_RES / 2));
	*y = *z * tan(CAMERA_PIX2RAD * (yPixels - CAMERA_Y_RES / 2));
	std::cout << "alpha " << alpha  << "  gamma " << gamma << std::endl;
	return;	
}

void getBallPosition(int xPixels, int yPixels, double diameterPixels, double* x, double* y, double* z)
{
	getBallPositionFromCamera(xPixels, yPixels, diameterPixels, x, y, z);
	*z += CAMERA_HEIGHT ;
	*x += CAMERA_X ;
	*y += CAMERA_Y ;
	return;
}



/*int main(int argc, char const *argv[])
{
	double x,y,z;
	getBallPosition(623, 821,20,&x, &y, &z);
	std::cout << "x : " << x << "  y: " << y << "  z : " << z << std::endl;
	return 0;
}*/