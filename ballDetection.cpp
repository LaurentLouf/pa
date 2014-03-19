#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/gpu/gpu.hpp"

#include <iostream>
#include <ctype.h>
#include <ctime>
#include <sstream>
#include <cmath>

#define DISPLAY
#define HOUGHh

int HUE_CHANNEL = 0 ;
int SATURATION_CHANNEL = 2 ; 
int LIGHTNESS_CHANNEL = 1 ;

using namespace cv;
using namespace std;

int H_min = 80 ;
int H_max = 95 ;
int S_min = 0 ;
int S_max = 255 ;
int V_min = 0 ;
int V_max = 255 ;

static const char WINDOW_ORIGIN[] = "Original Image";
static const char WINDOW_THRESHOLD[] = "Image with threshold applied";
static const char WINDOW_THRESHOLD_NOISE[] = "Image with threshold applied and noise cancelation";
static const char WINDOW_THRESHOLD_NOISE_BLUR[] = "Image with threshold applied and noise cancelation and blur";
static const char WINDOW_CONFIG[] = "Configuration";

static const char WINDOW_HUE[] = "Channel Hue";
static const char WINDOW_LIGHT[] = "Channel Lightness";
static const char WINDOW_SATURATION[] = "Channel Saturation";

static const char TRACKBAR_HUE_MIN[] = "Hue min";
static const char TRACKBAR_HUE_MAX[] = "Hue max";

static const char TRACKBAR_SATURATION_MIN[] = "Saturation min";
static const char TRACKBAR_SATURATION_MAX[] = "Saturation max";

static const char TRACKBAR_VALUE_MIN[] = "Value min";
static const char TRACKBAR_VALUE_MAX[] = "Value max";


Mat hlsChannels[3];


static void onMouse( int event, int x, int y, int /*flags*/, void* channel )
{
    if( event == CV_EVENT_LBUTTONDOWN )
    {
        int chan = *((int*) channel) ;
        int value  = hlsChannels[chan].at<int>((int) x, (int) y) ;
        cout << "X : " << x << ", Y : " << y << " " << chan << " value clicked " << value << endl ;
    }
}


int main( int argc, char** argv )
{
    // Image processing vars
    Mat image, imageFiltered, imageTmp ;
    vector<Vec3f> circles;
    vector<vector<Point> > contours;
    vector<Vec4i> contoursHierarchy;
    vector<vector<Point> >::iterator iteratorContours ;
    Mat erodeElement  = getStructuringElement( MORPH_RECT, Size ( 5, 5) );
    Mat dilateElement = getStructuringElement( MORPH_RECT, Size ( 5, 5) );

    // Circle vars
    Point2f center ;
    float radius ;

    // Video vars
    VideoCapture cap;
    TermCriteria termcrit(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS, 20, 0.03);
    Size subPixWinSize(10,10), winSize(31,31);

    // Display vars
    float begin, end, beginFiltering, endFiltering, filteringTime = 0, beginNoise, endNoise, noiseTime = 0, beginContour, endContour, contourTime = 0, calculationTime = 0, worst = 0, best = 100000000000000     ;
    int i = 0 ;
    stringstream stringOutput;


    // Get the video (filename or device)
    cap.open("MVI_7223.MOV");
    if( !cap.isOpened() )
    {
        cout << "Could not initialize capturing...\n";
        return 0;
    }

    // Create the windows
    namedWindow(WINDOW_ORIGIN) ;
    namedWindow(WINDOW_THRESHOLD);
    namedWindow(WINDOW_THRESHOLD_NOISE);
    namedWindow(WINDOW_THRESHOLD_NOISE_BLUR);  
    namedWindow(WINDOW_CONFIG);  

    createTrackbar(TRACKBAR_HUE_MIN, WINDOW_CONFIG, &H_min, 255) ;  
    createTrackbar(TRACKBAR_HUE_MAX, WINDOW_CONFIG, &H_max, 255) ;  
    createTrackbar(TRACKBAR_SATURATION_MIN, WINDOW_CONFIG, &S_min, 255) ;  
    createTrackbar(TRACKBAR_SATURATION_MAX, WINDOW_CONFIG, &S_max, 255) ;  
    createTrackbar(TRACKBAR_VALUE_MIN, WINDOW_CONFIG, &V_min, 255) ;  
    createTrackbar(TRACKBAR_VALUE_MAX, WINDOW_CONFIG, &V_max, 255) ;  

    moveWindow(WINDOW_ORIGIN, 0, 0) ;
    moveWindow(WINDOW_THRESHOLD, 0, 0);
    moveWindow(WINDOW_THRESHOLD_NOISE, 0, 0);
    moveWindow(WINDOW_THRESHOLD_NOISE_BLUR, 0, 0);
    moveWindow(WINDOW_CONFIG, 0, 0);

    for(;;)
    {
        begin = (float) clock() ;
        i++ ;


        Mat frame;
        cap >> frame;
        if( frame.empty() )
            break;

        frame.copyTo(image);

        // Filtering part : conversion from one color space to another
        beginFiltering = (float) clock() ;
        cvtColor(image, imageFiltered, CV_RGB2HSV);

        // Filtering part : Application of the threshold to keep only the color of the ball
        inRange(imageFiltered, Scalar(H_min, V_min, S_min), Scalar(H_max, V_max, S_max), imageFiltered); 
        endFiltering = (float) clock() ;

        #ifdef DISPLAY
            imshow(WINDOW_ORIGIN, image);
            imshow(WINDOW_THRESHOLD, imageFiltered) ;  
        #endif     


        // Noise cancelation
        beginNoise = (float) clock() ;
        #ifdef GPU
            gpu::GpuMat gpuImage(imageFiltered), gpuCircles ;

            gpu::erode ( gpuImage, gpuImage, erodeElement );
            gpu::erode ( gpuImage, gpuImage, erodeElement );
            gpu::dilate ( gpuImage, gpuImage, dilateElement );
            gpu::dilate ( gpuImage, gpuImage, dilateElement );
        #else
            erode ( imageFiltered, imageFiltered, erodeElement );
            dilate ( imageFiltered, imageFiltered, dilateElement );
            dilate ( imageFiltered, imageFiltered, dilateElement );
            erode ( imageFiltered, imageFiltered, erodeElement );
        #endif
        endNoise = (float) clock() ;
        
        #ifdef DISPLAY
            cout << "New image" << endl ;
            imshow(WINDOW_THRESHOLD_NOISE, imageFiltered) ;     
        #endif  


        // Contour determination
        beginContour = (float) clock() ;
        #ifdef GPU 
            gpu::threshold(gpuImage, gpuImage, 100, 255, THRESH_BINARY) ;
            gpu::HoughCircles(gpuImage, gpuCircles, CV_HOUGH_GRADIENT, 2, 100, 128, 100, 0, 400, 10) ;

            vector<Vec3f>::const_iterator itc = gpuCircles.begin();

           while (itc!=gpuCircles.end()) 
           {
                #ifdef DISPLAY
                    cout << "Circle found : (" << round((*itc)[0]) << ", " << round((*itc)[1])) << " | " << (*itc)[2] << ")" <<endl ;
                #else
                    stringOutput << "Circle found : (" << round((*itc)[0]) << ", " << round((*itc)[1])) << " | " << (*itc)[2] << ")" <<endl ;
                #endif
                 
                ++itc;
           }

        #else
            
            #ifdef HOUGH 
                threshold(imageFiltered, imageFiltered, 100, 255, THRESH_BINARY) ;
                HoughCircles(imageFiltered, circles, CV_HOUGH_GRADIENT, 1, 100, 128, 1000, 0, 400 ) ;

                for( size_t iCircle = 0; iCircle < circles.size(); iCircle++ )
                {

                    #ifdef DISPLAY
                        cout << "Circle found : (" << round(circles[iCircle][0]) << ", " << round(circles[iCircle][1]) << " | " << circles[iCircle][2] << ")" <<endl ;
                    #else
                        stringOutput << "Circle found : (" << round(circles[iCircle][0]) << ", " << round(circles[iCircle][1]) << " | " << circles[iCircle][2] << ")" <<endl ;
                    #endif
               }
            #else   
                threshold(imageFiltered, imageFiltered, 100, 255, THRESH_BINARY) ;
                findContours( imageFiltered, contours, contoursHierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );


                // Find the convex hull object for each contour
                vector<vector<Point> >hull( contours.size() );
                for( int i = 0; i < contours.size(); i++ )
                {  
                    convexHull( Mat(contours[i]), hull[i], false ); 

                    #ifdef DISPLAY
                        // Draw contours + hull results
                        drawContours( imageFiltered, contours, i, Scalar(255,255,255), 1, 8, vector<Vec4i>(), 0, Point() );
                        drawContours( imageFiltered, hull, i, Scalar(255,255,255), 1, 8, vector<Vec4i>(), 0, Point() );
                    #endif
                }



                for ( iteratorContours = contours.begin() ; iteratorContours != contours.end() ; iteratorContours++)
                {   
                    minEnclosingCircle ( *iteratorContours, center, radius ) ;

                    #ifdef DISPLAY
                        cout << "Circle found : (" << round(center.x) << ", " << round(center.y) << " | " << radius << ")" <<endl ;
                        Point centerCircle(cvRound(center.x), cvRound(center.y)) ;
                        circle( imageFiltered, centerCircle, cvRound(radius), Scalar(255,255,255), 2, 8, 0 );
                    #else
                        stringOutput << "Circle found : (" << round(center.x) << ", " << round(center.y) << " | " << radius << ")" <<endl ;
                    #endif
                } 
            #endif
        #endif
        endContour = (float) clock() ;


        #ifdef DISPLAY
            imshow(WINDOW_THRESHOLD_NOISE_BLUR, imageFiltered) ; 
        #endif 


        end = (float) clock() ;
        calculationTime += end - begin ;
        filteringTime += endFiltering - beginFiltering ;
        noiseTime += endNoise - beginNoise ;
        contourTime += endContour - beginContour ;

        if ( (end - begin) > worst )
            worst = end - begin ;
        else if ( (end - begin)  < best )
            best = end - begin ;



        char c = (char)waitKey(5);
        if( c == 27 )
        {
            break;
        }
            
    }

    #ifndef DISPLAY
        cout << stringOutput.str() ;
    #endif

    cout << "Filtering executed in " << filteringTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ")" << endl ;
    cout << "Noise cancelation executed in " << noiseTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ")" << endl ;
    cout << "Contour determination executed in " << contourTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ")" << endl ;
    cout << "Code executed in " << calculationTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ",best : " << best / ((float) CLOCKS_PER_SEC) << "s, worst : " << worst / ((float) CLOCKS_PER_SEC) << "s)" << endl ;

    return 0;
}
