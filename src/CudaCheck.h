#pragma once

#include <cuda_runtime.h>

bool checkCuda()
{
    int count = 0;
    cudaError_t error = cudaGetDeviceCount(&count);

    return error == cudaSuccess;
}
