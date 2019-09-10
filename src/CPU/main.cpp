//
// File for test function
//
#include <iostream>

#include "PrepareImage.h"
#include "ImageQuaolity.h"
#include "DCT_FILT.h"
#include <ctime>
#include <opencv2/opencv.hpp>


int main(int args, char **argv) {
    /*if(args != 3)
    {
        std::cout << "Enter the image name!" << "\n";
        return -2;
    }*/

    cv::Mat image;
    //image = cv::imread(argv[1], 1);
    image = cv::imread("../imset/Zoneplate.png", 1);

    if(!image.data) {
        return -1;
    }

    std::cout << "Image cols: " << image.cols << "\nImage rows: " << image.rows << "\nImage channels: " << image.channels() << "\n";

    cv::namedWindow("image", cv::WINDOW_AUTOSIZE);
    cv::imshow("image", image);
    cv::waitKey(0);

    double *pixel_ride = new double[image.cols*image.rows];
    int pixel_flag = 0;
    cv::Mat RowClone;
    cv::Mat A(image.size(), image.type());

    double time = 0.0;
    double SMD_image = 0;

    for(int channel = 0; channel < image.channels(); channel++) {
        pixel_flag = 0;
        std::cout << "Channel " << channel << " in process" << "\n";
        for (int i = 0; i < image.rows; i++) {
            RowClone = image.row(i).clone();
            for (int j = channel; j < RowClone.cols * image.channels(); j += image.channels()) {
                *(pixel_ride + pixel_flag) = (double) *(RowClone.data + j);
                pixel_flag++;
            }
        }

        SMD_image = ImProcessing::SMD(pixel_ride, image.cols, image.rows);
        std::cout << "Start SMD in " << channel << " channel = " << SMD_image << "\n";
        ImProcessing::DCT_FILT filter(pixel_ride, image.cols, image.rows, 8, 0.04); // JPEG = 0.012
        unsigned int start_time = clock();
        double *pixel_result = filter.imageFilter();
        unsigned int finish_time = clock();
        std::cout << "Time monothread = " << (float) (finish_time - start_time) / CLOCKS_PER_SEC << " s\n";
        time += (float) (finish_time - start_time) / CLOCKS_PER_SEC;

        SMD_image = ImProcessing::SMD(pixel_result, image.cols, image.rows);
        std::cout << "SMD after precess " << channel << " channel = " << SMD_image << "\n";

        pixel_flag = 0;
        for (int i = 0; i < image.rows; i++)
            for(int j = 0; j < image.cols; j++)
            {
                A.data[A.channels() * (A.cols * i + j) + channel] = *(pixel_result + pixel_flag);
                pixel_flag++;
            }
    }

    cv::namedWindow("filtered", cv::WINDOW_AUTOSIZE);
    cv::imshow("filtered", A);
    cv::waitKey(0);

    std::cout << "All processing time = " << time << " s\n";

    return 0;
}