#include "postDetection.h"

using namespace std;

void getBallPositionFromCamera(int xPixels, int yPixels, double diameterPixels, BallState& ballState)
{
	ballState.oldx = ballState.x;
	ballState.oldy = ballState.y;
	ballState.oldz = ballState.z;
	double alpha = abs(diameterPixels * CAMERA_PIX2RAD);
	double gamma = abs(CAMERA_PIX2RAD * sqrt(xPixels*xPixels + yPixels*yPixels));
	ballState.z = - BALL_DIAMETER / (2*( tan(alpha/2 + gamma) - tan(gamma)));
	ballState.x = ballState.z * tan(CAMERA_PIX2RAD * (xPixels - CAMERA_X_RES / 2));
	ballState.y = ballState.z * tan(CAMERA_PIX2RAD * (yPixels - CAMERA_Y_RES / 2));
	//std::cout << "alpha " << alpha  << "  gamma " << gamma << std::endl;
	return;	
}

void getBallPosition(int xPixels, int yPixels, double diameterPixels, BallState& ballState)
{
	getBallPositionFromCamera(xPixels, yPixels, diameterPixels, ballState);
	ballState.z += CAMERA_HEIGHT ;
	ballState.x += CAMERA_X ;
	ballState.y += CAMERA_Y ;
	return;
}

void calculateBallSpeed(BallState& ballState)
{
	ballState.vx = ballState.x - ballState.oldx;
	ballState.vy = ballState.y - ballState.oldy;
	ballState.vz = ballState.z - ballState.oldz;
}

CircleFound getBestCircle(std::vector<CircleFound> const & circlesFound)
{
	// basic : we get the biggest circle
	CircleFound tmpCircle;
	int n = circlesFound.size();
	for(int i =0 ; i < n ; i++)
	{
		if(circlesFound[i].radius > tmpCircle.radius)
			tmpCircle = circlesFound[i];
	}
	if(n == 0)
		tmpCircle.radius = 0;
	return tmpCircle;
}



/*int main(int argc, char const *argv[])
{
	double x,y,z;
	getBallPosition(623, 821,20,&x, &y, &z);
	std::cout << "x : " << x << "  y: " << y << "  z : " << z << std::endl;
	return 0;
}*/