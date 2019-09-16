//
// Created by rostislav on 10.09.19.
//

#include <cstdlib>
#include "RAWImage.h"

int main(int args, char **argv)
{
    bool file_folder = true; //true == file
    bool cpu_gpu = true; //true == cpu
    bool show_image = true;
    std::string file = "";
    std::string folder = "";
    int block_size = 8;

    for(int i = 1; i < args;) {
        if (std::string(argv[i]) == "-file") {
            file_folder = true;
            file = argv[i + 1];
            i += 2;
        } else if (std::string(argv[i]) == "-folder") {
            file_folder = false;
            folder = argv[i + 1];
            i += 2;
        } else if (std::string(argv[i]) == "-window") {
            if (std::string(argv[i+1]) == "2")
                block_size = 2;
            else if (std::string(argv[i+1]) == "4")
                block_size = 4;
            else if (std::string(argv[i+1]) == "32")
                block_size = 32;
            else if (std::string(argv[i+1]) == "16")
                block_size = 16;
            else
                block_size = 8;
            i += 2;
        } else if (std::string(argv[i]) == "-cpu") {
            cpu_gpu = true;
            i++;
        } else if (std::string(argv[i]) == "-gpu"){
            cpu_gpu = false;
            i++;
        } else{
            i++;
        }
    }

    if(file_folder) {
        using namespace ImProcessing;
        RAWImage image(file, TYPE_BGR);

        image.ApplyDCT(block_size, 0.012, true);

        //float *coeff = image.DCTCoefficients();
        //for (int i = 0; i < 16; i++) printf("%f ", coeff[i]);
        image.show();
    }

    return 0;
}


