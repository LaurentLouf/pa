#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/gpu/gpu.hpp"
#include "opencv2/ocl/ocl.hpp"
#include "opencv2/nonfree/ocl.hpp"

#include <iostream>
#include <ctype.h>
#include <ctime>
#include <sstream>
#include <cmath>
#include <sys/time.h>

//#define DISPLAY
#define HOUGHh
#define OCL

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

#ifdef DISPLAY
    static const char WINDOW_ORIGIN[] = "Original Image";
    static const char WINDOW_THRESHOLD[] = "Image with threshold applied";
    static const char WINDOW_THRESHOLD_NOISE[] = "Image with threshold applied and noise cancelation";
    static const char WINDOW_THRESHOLD_NOISE_BLUR[] = "Image with threshold applied and noise cancelation and blur";
    static const char WINDOW_CONFIG[] = "Configuration";
    static const char WINDOW_TEST[] = "Test";
#endif

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
    // Benchmark variables
    struct timeval tv_begin, tv_end, beginTime, endTime;
    double filteringTime=0 , noiseTime=0, contourTime=0, readTime=0, totalTime; 

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
    //double begin, end, beginFiltering, endFiltering, filteringTime = 0, beginNoise, endNoise, noiseTime = 0, beginContour, endContour, contourTime = 0, calculationTime = 0, worst = 0, best = 100000000000000     ;
    int i = 0 ;
    stringstream stringOutput;


    // Get the video (filename or device)
    cap.open("MVI_7319.MOV");
    if( !cap.isOpened() )
    {
        cout << "Could not initialize capturing...\n";
        return 0;
    }

    #ifdef DISPLAY
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

        namedWindow(WINDOW_TEST);
        moveWindow(WINDOW_TEST, 0, 0);
    #endif

    #ifdef OCL
        ocl::PlatformsInfo platforms;
        ocl::getOpenCLPlatforms(platforms);
        ocl::DevicesInfo devices;
        ocl::getOpenCLDevices(devices);
        std::cout << "platforms " << platforms.size() << "  devices " << devices.size() << " " << devices[0]->deviceDriverVersion << std::endl;
        ocl::setDevice(devices[0]);

    #endif

    gettimeofday(&beginTime, NULL);

    for(;;)
    {
        //begin = (float) clock() ;
        i++ ;

        gettimeofday(&tv_begin, NULL);
        Mat frame;
        cap >> frame;
        if( frame.empty() )
            break;

        frame.copyTo(image);
        gettimeofday(&tv_end, NULL);
        readTime += (double) (tv_end.tv_sec - tv_begin.tv_sec) + ((double) (tv_end.tv_usec - tv_begin.tv_usec)/1000000);

        // Filtering part : conversion from one color space to another
        //beginFiltering = (float) clock() ;
        gettimeofday(&tv_begin, NULL);
        #ifdef OCL
            ocl::oclMat oclImage(image), oclImageFilteredCvt, oclImageFiltered, oclCircles, thresholdMin, thresholdMax;
            vector<ocl::oclMat> channels;
            ocl::cvtColor(oclImage, oclImageFilteredCvt, CV_RGB2HSV);

            ocl::split(oclImageFilteredCvt, channels);
            ocl::threshold(channels[0], thresholdMin, H_min, 255, THRESH_BINARY);
            ocl::threshold(channels[0], thresholdMax, H_max, 255, THRESH_BINARY_INV);
            ocl::bitwise_and(thresholdMin, thresholdMax, oclImageFiltered); //imageTmp = oclImageFiltered;imshow(WINDOW_TEST, imageTmp);
            
            //inRange(oclImageFiltered, Scalar(H_min, V_min, S_min), Scalar(H_max, V_max, S_max), oclImageFiltered);
            //imageFiltered = oclImageFiltered;
            //inRange seems not to work with ocl ...
            //inRange(imageFiltered, Scalar(H_min, V_min, S_min), Scalar(H_max, V_max, S_max), imageFiltered);
        #else
            cvtColor(image, imageFiltered, CV_RGB2HSV);

            // Filtering part : Application of the threshold to keep only the color of the ball
            inRange(imageFiltered, Scalar(H_min, V_min, S_min), Scalar(H_max, V_max, S_max), imageFiltered);
        #endif 
        //endFiltering = (float) clock() ;
        gettimeofday(&tv_end, NULL);
        filteringTime += (double) (tv_end.tv_sec - tv_begin.tv_sec) + ((double) (tv_end.tv_usec - tv_begin.tv_usec)/1000000);


        #ifdef DISPLAY
            #ifdef OCL
                imageFiltered = oclImageFiltered;
            #endif
            imshow(WINDOW_ORIGIN, image);
            imshow(WINDOW_THRESHOLD, imageFiltered) ;  
        #endif     


        // Noise cancelation
        //beginNoise = (float) clock() ;
        gettimeofday(&tv_begin, NULL);
        #ifdef OCL
            //oclImageFiltered = imageFiltered;

            ocl::erode ( oclImageFiltered, oclImageFiltered, erodeElement );
            ocl::erode ( oclImageFiltered, oclImageFiltered, erodeElement );
            ocl::dilate ( oclImageFiltered, oclImageFiltered, dilateElement );
            ocl::dilate ( oclImageFiltered, oclImageFiltered, dilateElement );
            // Peut être à supprimer par la suite
            imageFiltered = oclImageFiltered;
        #else
            erode ( imageFiltered, imageFiltered, erodeElement );
            dilate ( imageFiltered, imageFiltered, dilateElement );
            dilate ( imageFiltered, imageFiltered, dilateElement );
            erode ( imageFiltered, imageFiltered, erodeElement );
        #endif
        //endNoise = (float) clock() ;
        gettimeofday(&tv_end, NULL);
        noiseTime += (double) (tv_end.tv_sec - tv_begin.tv_sec) + ((double) (tv_end.tv_usec - tv_begin.tv_usec)/1000000);
        
        #ifdef DISPLAY
            cout << "New image" << endl ;
            imshow(WINDOW_THRESHOLD_NOISE, imageFiltered) ;     
        #endif  


        // Contour determination
        //beginContour = (float) clock() ;
        gettimeofday(&tv_begin, NULL);
        #ifdef OCLss
            //source de lenteur ici
            ocl::threshold(oclImageFiltered, oclImageFiltered, 100, 255, THRESH_BINARY) ;
            imageFiltered = oclImageFiltered;
            #ifdef HOUGH
                HoughCircles(imageFiltered, circles, CV_HOUGH_GRADIENT, 1, 100, 128, 1000, 0, 400 ) ;
            #else
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
        //endContour = (float) clock() ;
        gettimeofday(&tv_end, NULL);
        contourTime += (double) (tv_end.tv_sec - tv_begin.tv_sec) + ((double) (tv_end.tv_usec - tv_begin.tv_usec)/1000000);


        #ifdef DISPLAY
            imshow(WINDOW_THRESHOLD_NOISE_BLUR, imageFiltered) ; 
        #endif 


        /*end = (float) clock() ;
        calculationTime += end - begin ;
        filteringTime += endFiltering - beginFiltering ;
        noiseTime += endNoise - beginNoise ;
        contourTime += endContour - beginContour ;

        if ( (end - begin) > worst )
            worst = end - begin ;
        else if ( (end - begin)  < best )
            best = end - begin ;*/



        char c = (char)waitKey(5);
        if( c == 27 )
        {
            break;
        }
            
    }
    gettimeofday(&endTime, NULL);
    totalTime = (double) (endTime.tv_sec - beginTime.tv_sec) + ((double) (endTime.tv_usec - beginTime.tv_usec)/1000000);
    

    #ifndef DISPLAY
        //cout << stringOutput.str() ;
    #endif

    cout << "Total Time : " << totalTime << "s" << endl;
    cout << " Read time : " << readTime << "s" << endl;
    cout << " Filtering time : " << filteringTime << "s" << endl;
    cout << " Noise cancellation time : " << noiseTime << "s" << endl;
    cout << " Contour determination time : " << contourTime << "s" << endl;

    /*cout << "Filtering executed in " << filteringTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ")" << endl ;
    cout << "Noise cancelation executed in " << noiseTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ")" << endl ;
    cout << "Contour determination executed in " << contourTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ")" << endl ;
    cout << "Code executed in " << calculationTime / ((float) i) / ((float) CLOCKS_PER_SEC) << "s on average per frame (frames : " << i << ",best : " << best / ((float) CLOCKS_PER_SEC) << "s, worst : " << worst / ((float) CLOCKS_PER_SEC) << "s)" << endl ;*/

    return 0;
}
