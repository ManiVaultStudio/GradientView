#pragma once

#include <cuda_runtime.h>
#include <iostream>

bool checkCuda()
{
    typedef cudaError_t(*FnGetDeviceCount)    (int* count);
    HMODULE hCuda = LoadLibrary(TEXT("cudart64_12.dll"));
    if (!hCuda) { std::cout << "ERROR Cannot load Cuda DLL" << std::endl; return false; }// ERROR: cannot load dll, DllMain must have failed because cudart only depends on Kernel dll implicitly. Or cannot find dll in curent directory or in the path.
    FnGetDeviceCount fnGetDeviceCount = (FnGetDeviceCount)GetProcAddress(hCuda, "cudaGetDeviceCount");
    if (!fnGetDeviceCount) { std::cout << "ERROR CudaRT has no entry point for cudaGetDeviceCount" << std::endl; return false; } // ERROR: cudart has no entry point for cudaGetDeviceCount ?!
    int count = 0;
    cudaError_t err = (*fnGetDeviceCount)(&count);
    if (cudaSuccess != err) { std::cout << "ERROR: " << err << " Device enumeration failed " << std::endl; return false; };// ERROR: we don't wanna use CUDA if even device enumeration fails
    if (!count) { std::cout << "ERROR CUDA has no devices" << std::endl; return false; }; // FALLBACK: CUDA has no devices, don't try to use it, fallback to some other BLAS

    std::cout << ">>>>>>>>> CUDA loaded successfully, using accelerated KNN" << std::endl;
    return true;
    //int count = 0;
    //cudaError_t error = cudaGetDeviceCount(&count);

    //return error == cudaSuccess;
}
