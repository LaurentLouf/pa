#ifndef IMAGEPROCESSING
#define IMAGEPROCESSING

#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/gpu/gpu.hpp"
#include "opencv2/ocl/ocl.hpp"
#include "opencv2/nonfree/ocl.hpp"

#include <sys/time.h>

#define OCL
#define DISPLAY
#define HOUGHh

using namespace cv;
using namespace std;

typedef struct{
	double filteringTime=0, noiseTime=0, contourTime=0, readTime=0, waitFinal=0, totalTime, finalLoop; 
	struct timeval tv_begin, tv_end, beginTime, endTime;
}Benchmark;

typedef struct{
	int x;
	int y;
	double radius;
}CircleFound;

extern int HUE_CHANNEL ;
extern int SATURATION_CHANNEL ;
extern int LIGHTNESS_CHANNEL ;

extern int H_min ;
extern int H_max ;
extern int S_min ;
extern int S_max ;
extern int V_min ;
extern int V_max ;

extern Mat hlsChannels[3];
extern stringstream stringOutput;

extern const char WINDOW_ORIGIN[];
extern const char WINDOW_THRESHOLD[];
extern const char WINDOW_THRESHOLD_NOISE[];
extern const char WINDOW_THRESHOLD_NOISE_BLUR[];
extern const char WINDOW_CONFIG[];
extern const char WINDOW_TEST[];


extern const char WINDOW_HUE[];
extern const char WINDOW_LIGHT[];
extern const char WINDOW_SATURATION[];

extern const char TRACKBAR_HUE_MIN[];
extern const char TRACKBAR_HUE_MAX[];

extern const char TRACKBAR_SATURATION_MIN[];
extern const char TRACKBAR_SATURATION_MAX[];

extern const char TRACKBAR_VALUE_MIN[];
extern const char TRACKBAR_VALUE_MAX[];


int findBall(VideoCapture& cap, Benchmark& bench, std::vector<CircleFound>& circlesFound); // return 0 if all is ok



#endif // IMAGEPROCESSING