//
// Created by rostislav on 10.09.19.
//
#include "cstdio"
#include <string>

#ifndef CUDA__ERROR_PRINTER_H
#define CUDA__ERROR_PRINTER_H

namespace ImProcessing{
    void PrintError(std::string text)
    {
        printf("Error: %s\n", text[0]);
    }
}

#endif //CUDA__ERROR_PRINTER_H
