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

void getPlayerPosition(Position& playerPosition)
{
	playerPosition.x = 2.2;
	playerPosition.y = -3.6;
	playerPosition.z = 0;
	return;
}

void setBallTarget(Position& playerPosition, louf& ballTarget)
{
	ballTarget.zf = 0;
	if(DIFFICULTY == 0) // on renvoie la balle devant le joueur
	{
		ballTarget = getLoufStruct(playerPosition);
		ballTarget.yf += 1.5; 
	}
	else if(DIFFICULTY == 1) // on renvoie à gauche si joueur à droite, et inveersement
	{
		ballTarget.yf = 0;
		if(playerPosition.x > 0)
			ballTarget.xf = -3;
		else
			ballTarget.xf = 3;
	}
	else if(DIFFICULTY == 2)
	{
		if(playerPosition.y < 0)
			ballTarget.yf = 4;
		else
			ballTarget.yf = -4;

		if(playerPosition.x < 0)
			ballTarget.xf = 3;
		else
			ballTarget.xf = -3;
	}

}

louf getLoufStruct(Position& position)
{
	louf output;
	output.xf = position.x;
	output.yf = position.y;
	output.zf = position.z;
	return output;
}

void copyStateToMaths(BallState & initial_t1, BallState & initial_t2, manchoul& output)
{
	output.xi = initial_t2.x;
	output.yi = initial_t2.y;
	output.zi = initial_t2.z;
	output.vxi = initial_t2.vx;
	output.vyi = initial_t2.vy;
	output.vzi = initial_t2.vz;
	output.wxi = initial_t1.vx;
	output.wyi = initial_t1.vy;
	output.wzi = initial_t1.vz;
}



/*int main(int argc, char const *argv[])
{
	double x,y,z;
	getBallPosition(623, 821,20,&x, &y, &z);
	std::cout << "x : " << x << "  y: " << y << "  z : " << z << std::endl;
	return 0;
}*/