#include "imageProcessing.h"


int findBall(VideoCapture& cap, Benchmark& bench, std::vector<CircleFound>& circlesFound)
{
	// Image processing vars
    Mat imageFiltered, imageTmp ;
    vector<Vec3f> circles;
    vector<vector<Point> > contours;
    vector<Vec4i> contoursHierarchy;
    vector<vector<Point> >::iterator iteratorContours ;
    Mat erodeElement  = getStructuringElement( MORPH_RECT, Size ( 5, 5) );
    Mat dilateElement = getStructuringElement( MORPH_RECT, Size ( 5, 5) );
    circlesFound.clear();

  	// Circle vars
    Point2f center ;
    float radius ;

    gettimeofday(&bench.tv_begin, NULL);
    Mat image;
    bool isStreamOK = cap.read(image);
    if( !isStreamOK )
    {
        return -1;
    }
    //frame.copyTo(image); //inutile et couteux en temps (env 1.5s)
    gettimeofday(&bench.tv_end, NULL);
    bench.readTime += (double) (bench.tv_end.tv_sec - bench.tv_begin.tv_sec) + ((double) (bench.tv_end.tv_usec - bench.tv_begin.tv_usec)/1000000);

    // Filtering part : conversion from one color space to another
    gettimeofday(&bench.tv_begin, NULL);
    #ifdef OCL
        ocl::oclMat oclImage(image), oclImageFilteredCvt, oclImageFiltered, oclCircles, thresholdMin, thresholdMax;
        vector<ocl::oclMat> channels;
        ocl::cvtColor(oclImage, oclImageFilteredCvt, CV_RGB2HSV);

        // 4 lines above equivalent to the function inRange of non-OCL code, but inRange doesn't exist in OCL so I had to "recode" it
        ocl::split(oclImageFilteredCvt, channels);
        ocl::threshold(channels[0], thresholdMin, H_min, 255, THRESH_BINARY);
        ocl::threshold(channels[0], thresholdMax, H_max, 255, THRESH_BINARY_INV);
        ocl::bitwise_and(thresholdMin, thresholdMax, oclImageFiltered); 
    #else
        cvtColor(image, imageFiltered, CV_RGB2HSV);
        // Filtering part : Application of the threshold to keep only the color of the ball
        inRange(imageFiltered, Scalar(H_min, V_min, S_min), Scalar(H_max, V_max, S_max), imageFiltered);
        // Essai pour ne récupérer que le premier plan, à creuser pour détecter le joueur car lent.
        /*//Ptr<BackgroundSubtractor> pMog;
        BackgroundSubtractorMOG pMOG;
        //pMog = BackgroundSubtractorMOG();
        pMOG(image, imageFiltered);*/
    #endif 
    gettimeofday(&bench.tv_end, NULL);
    bench.filteringTime += (double) (bench.tv_end.tv_sec - bench.tv_begin.tv_sec) + ((double) (bench.tv_end.tv_usec - bench.tv_begin.tv_usec)/1000000);

    #ifdef DISPLAY
        #ifdef OCL
            imageFiltered = oclImageFiltered;
        #endif
        imshow(WINDOW_ORIGIN, image);
        imshow(WINDOW_THRESHOLD, imageFiltered) ;  
    #endif     

    // Noise cancelation
    gettimeofday(&bench.tv_begin, NULL);
    #ifdef OCL

        ocl::erode ( oclImageFiltered, oclImageFiltered, erodeElement );
        ocl::erode ( oclImageFiltered, oclImageFiltered, erodeElement );
        ocl::dilate ( oclImageFiltered, oclImageFiltered, dilateElement );
        ocl::dilate ( oclImageFiltered, oclImageFiltered, dilateElement );
        // Peut être à supprimer par la suite (si contourdetermination en OCL efficace)
        imageFiltered = oclImageFiltered;
    #else
        erode ( imageFiltered, imageFiltered, erodeElement );
        dilate ( imageFiltered, imageFiltered, dilateElement );
        dilate ( imageFiltered, imageFiltered, dilateElement );
        erode ( imageFiltered, imageFiltered, erodeElement );
    #endif
    gettimeofday(&bench.tv_end, NULL);
    bench.noiseTime += (double) (bench.tv_end.tv_sec - bench.tv_begin.tv_sec) + ((double) (bench.tv_end.tv_usec - bench.tv_begin.tv_usec)/1000000);
    
    #ifdef DISPLAY
        cout << "New image" << endl ;
        imshow(WINDOW_THRESHOLD_NOISE, imageFiltered) ;     
    #endif  

    // Contour determination
    gettimeofday(&bench.tv_begin, NULL);
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
                CircleFound tmpCircle;
                tmpCircle.x = round(center.x);
                tmpCircle.y = round(center.y);
                tmpCircle.radius = radius;
                circlesFound.push_back(tmpCircle);
                #ifdef DISPLAY
                    cout << "Circle found : (" << round(center.x) << ", " << round(center.y) << " | " << radius << ")" <<endl ;
                    Point centerCircle(cvRound(center.x), cvRound(center.y)) ;
                    circle( imageFiltered, centerCircle, cvRound(radius), Scalar(255,255,255), 2, 8, 0 );
                #else
                    //stringOutput << "Circle found : (" << round(center.x) << ", " << round(center.y) << " | " << radius << ")" <<endl ;
                #endif
            } 

        #endif

    #else
        
        #ifdef HOUGH 
            threshold(imageFiltered, imageFiltered, 100, 255, THRESH_BINARY) ;
            HoughCircles(imageFiltered, circles, CV_HOUGH_GRADIENT, 1, 100, 128, 1000, 0, 400 ) ;

            for( size_t iCircle = 0; iCircle < circles.size(); iCircle++ )
            {
            	CircleFound tmpCircle;
                tmpCircle.x = round(circles[iCircle][0]);
                tmpCircle.y = round(circles[iCircle][1]);
                tmpCircle.radius = circles[iCircle][2];
                circlesFound.push_back(tmpCircle);
                #ifdef DISPLAY
                    cout << "Circle found : (" << round(circles[iCircle][0]) << ", " << round(circles[iCircle][1]) << " | " << circles[iCircle][2] << ")" <<endl ;
                #else
                    //stringOutput << "Circle found : (" << round(circles[iCircle][0]) << ", " << round(circles[iCircle][1]) << " | " << circles[iCircle][2] << ")" <<endl ;
                #endif
           }
        #else   
            threshold(imageFiltered, imageFiltered, 100, 255, THRESH_BINARY) ;
            findContours( imageFiltered, contours, contoursHierarchy, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0) );

            //#ifdef DISPLAY // inutile pour l'algo non ?
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
            //#endif



            for ( iteratorContours = contours.begin() ; iteratorContours != contours.end() ; iteratorContours++)
            {   
                minEnclosingCircle ( *iteratorContours, center, radius ) ;
                CircleFound tmpCircle;
                tmpCircle.x = round(center.x);
                tmpCircle.y = round(center.y);
                tmpCircle.radius = radius;
                circlesFound.push_back(tmpCircle);
                #ifdef DISPLAY
                    cout << "Circle found : (" << round(center.x) << ", " << round(center.y) << " | " << radius << ")" <<endl ;
                    Point centerCircle(cvRound(center.x), cvRound(center.y)) ;
                    circle( imageFiltered, centerCircle, cvRound(radius), Scalar(255,255,255), 2, 8, 0 );
                #else
                    //stringOutput << "Circle found : (" << round(center.x) << ", " << round(center.y) << " | " << radius << ")" <<endl ;
                #endif
            } 
        #endif
    #endif
    gettimeofday(&bench.tv_end, NULL);
    bench.contourTime += (double) (bench.tv_end.tv_sec - bench.tv_begin.tv_sec) + ((double) (bench.tv_end.tv_usec - bench.tv_begin.tv_usec)/1000000);

    #ifdef DISPLAY
        imshow(WINDOW_THRESHOLD_NOISE_BLUR, imageFiltered) ; 
    #endif 

    return 0;
}