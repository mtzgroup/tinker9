// Compile this file by
// nvcc ../tinkerbox_interface/test_tinkerbox_interface.cpp libtinkerbox.so -lcublas -ccbin g++ --compiler-options "-std=c++11" -o test_tinkerbox_interface;
// Run the executable by
// CUDA_VISIBLE_DEVICES=2 LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH ./test_tinkerbox_interface;

#include "tinkerbox.h"

#include <stdio.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>

template <typename T>
void check(T result, char const *const func, const char *const file,
           int const line) {
  if (result) {
    fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
            static_cast<unsigned int>(result), cudaGetErrorName(result), func);
    exit(EXIT_FAILURE);
  }
}
#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

#define CUBLASERR(STATUS) { if (STATUS != CUBLAS_STATUS_SUCCESS) {	\
  if (STATUS == CUBLAS_STATUS_NOT_INITIALIZED) {			\
    printf("CUBLAS error: CUBLAS_STATUS_NOT_INITIALIZED, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_ALLOC_FAILED) {			\
    printf("CUBLAS error: CUBLAS_STATUS_ALLOC_FAILED, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_INVALID_VALUE) {			\
    printf("CUBLAS error: CUBLAS_STATUS_INVALID_VALUE, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_ARCH_MISMATCH) {			\
  printf("CUBLAS error: CUBLAS_STATUS_ARCH_MISMATCH, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_MAPPING_ERROR) {			\
    printf("CUBLAS error: CUBLAS_STATUS_MAPPING_ERROR, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_EXECUTION_FAILED) {			\
    printf("CUBLAS error: CUBLAS_STATUS_EXECUTION_FAILED, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_INTERNAL_ERROR) {			\
    printf("CUBLAS error: CUBLAS_STATUS_INTERNAL_ERROR, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_NOT_SUPPORTED) {			\
    printf("CUBLAS error: CUBLAS_STATUS_NOT_SUPPORTED, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  else if (STATUS == CUBLAS_STATUS_LICENSE_ERROR) {			\
    printf("CUBLAS error: CUBLAS_STATUS_LICENSE_ERROR, file %s, line %d\n", __FILE__, __LINE__); \
  }									\
  exit(1);          \
}}

int main()
{
    cudaSetDevice(0);
    cudaSetDeviceFlags(cudaDeviceMapHost); // Henry 20220111: This is necessary.
    
    int M = 10, N = 10, K = 10;
    double* A, * B, * C;
    A = new double[M * K];
    B = new double[K * N];
    C = new double[M * N];

    for (int i = 0; i < M; i++)
        for (int j = 0; j < K; j++)
            A[i + j * M] = i * j;
    for (int i = 0; i < K; i++)
        for (int j = 0; j < N; j++)
            B[i + j * K] = i + j;

    double* d_A, * d_B, * d_C, * d_A_T, * d_work;
    checkCudaErrors(cudaMalloc((void**)&d_A, M * K * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_B, K * N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_C, M * N * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_A_T, M * K * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_work, M * N * sizeof(double)));

    checkCudaErrors(cudaMemcpy(d_A, A, M * K * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_B, B, K * N * sizeof(double), cudaMemcpyHostToDevice));

    {
        char* xyz = "~/test/tinker-banchmark/water_bulk.xyz";
        initialize_tinker(nullptr, 0, xyz);

        double energy = get_energy_nonpolar_mm_contribution();
        printf("energy = %.10f\n", energy);

        int n = get_n_mm();
        double* gradient = new double[n * 3];
        get_gradients_all_atoms_mm_contribution(gradient);

        for (int i = 0; i < n; i++)
            printf("Atom %d gradient %.10f, %.10f, %.10f \n", i, gradient[i*3+0], gradient[i*3+1], gradient[i*3+2]);
        delete[] gradient;
    }

    cublasHandle_t handle;
    cublasStatus_t stat;
    stat = cublasCreate(&handle);
    CUBLASERR(stat);

    double alpha = 1.0, beta = 0.0;
    stat = cublasDgemm(handle, CUBLAS_OP_T, CUBLAS_OP_N, M, M, M, &alpha, d_A, M, d_B, M, &beta, d_work, M);
    CUBLASERR(stat);

    stat = cublasDgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, M, M, M, &alpha, d_work, M, d_A, M, &beta, d_C, M);
    CUBLASERR(stat);

    checkCudaErrors(cudaMemcpy(C, d_C, M * N * sizeof(double), cudaMemcpyDeviceToHost));

    printf("Henry: cuda result = %.2f\n", C[M * N - 1]);

    stat = cublasDestroy_v2(handle);
    CUBLASERR(stat);

    checkCudaErrors(cudaFree(d_A));
    checkCudaErrors(cudaFree(d_B));
    checkCudaErrors(cudaFree(d_C));
    checkCudaErrors(cudaFree(d_A_T));
    checkCudaErrors(cudaFree(d_work));

    delete[] A;
    delete[] B;
    delete[] C;

    return 0;
}
