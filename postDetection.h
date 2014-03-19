#ifndef POSTDETECTION
#define POSTDETECTION

#include <cmath>
#include <iostream>

#define BALL_DIAMETER 0.065 // taille officielle entre 6.35 et 6.668 cm
#define CAMERA_HEIGHT 5
#define CAMERA_X 3
#define CAMERA_Y 4
#define CAMERA_ANGLE 0
#define CAMERA_X_RES 1920
#define CAMERA_Y_RES 1080
#define CAMERA_PIX2RAD 0.002

void getBallPositionFromCamera(int xPixels, int yPixels, double diameterPixels, double* x, double* y, double* z);
void getBallPosition(int xPixels, int yPixels, double diameterPixels, double* x, double* y, double* z);


#endif // POSTDETECTION