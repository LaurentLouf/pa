#ifndef POSTDETECTION
#define POSTDETECTION

#include <cmath>
#include <iostream>

#include "imageProcessing.h"

#define BALL_DIAMETER 0.065 // taille officielle entre 6.35 et 6.668 cm
#define CAMERA_HEIGHT 5
#define CAMERA_X 3
#define CAMERA_Y 4
#define CAMERA_ANGLE 0
#define CAMERA_X_RES 1920
#define CAMERA_Y_RES 1080
#define CAMERA_PIX2RAD 0.0006

typedef struct{
	double x,y,z;
	double oldx, oldy, oldz;
	double vx, vy, vz;
}BallState;

void getBallPositionFromCamera(int xPixels, int yPixels, double diameterPixels, double* x, double* y, double* z);
void getBallPosition(int xPixels, int yPixels, double diameterPixels, BallState& ballState);
void calculateBallSpeed(BallState& ballState);
CircleFound getBestCircle(std::vector<CircleFound> const & circlesFound);
//void setWallPosition( structure à def);



#endif // POSTDETECTION