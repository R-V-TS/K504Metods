//
// Created by rostislav on 10.09.19.
//

#include <cstdlib>
#include "RAWImage.h"
#include <ctime>
#include <fstream>

int main(int args, char **argv)
{
    bool file_folder = true; //true == file
    bool cpu_gpu = true; //true == cpu
    bool show_image = false;
    bool save_image = false;
    std::string file = "";
    std::string folder = "";
    std::string destination_folder = "";
    std::string method = "DCT";
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
        } else if(std::string(argv[i]) == "-save"){
            save_image = true;
            destination_folder = argv[i+1];
            i += 2;
        } else if(std::string(argv[i]) == "-method"){
            method = argv[i+1];
            i += 2;
        } else if(std::string(argv[i]) == "-show"){
            show_image = true;
            i++;
        }
        else{
            i++;
        }
    }

    unsigned int time_start = clock();
    if(file_folder) {
        using namespace ImProcessing;
        RAWImage image(file, TYPE_BGR);

        if(method == "DCT"){
            image.ApplyDCT(block_size, 0.012, cpu_gpu);
            if(show_image) image.show();
            if(save_image) image.save(destination_folder);
        }

        if(method == "DCT_frequency")
        {
            float *coeff = image.DCTCoefficients(cpu_gpu);
            for (int i = 0; i < 16; i++) printf("%f ", coeff[i]);
            if(save_image) {
                std::ofstream out;
                out.open(destination_folder + "frequency.txt", std::ios_base::app);
                for (int i = 0; i < 16; i++)
                    out << coeff[i] << "/n";
                out.close();
            }
        }

        if(method == "Test")
        {
            image.AddNoise(0, 2);
            float* coff = image.DCTCoefficients(cpu_gpu);
            for (int i = 0; i < 16; i++) printf("%f ", coff[i]);
            printf("\n");
            image.printImageCharacteristics();
            image.ApplyDCT(8,0.012, cpu_gpu);
            float* coeff = image.DCTCoefficients(cpu_gpu);
            for (int i = 0; i < 16; i++) printf("%f ", coeff[i]);
            printf("\n");
            image.printImageCharacteristics();
        }
    }
    unsigned int finish_time = clock();
    printf("\nProgram works %f s\n", (float) (finish_time-time_start)/CLOCKS_PER_SEC);


    return 0;
}


