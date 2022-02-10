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

    // All the stuff outside this block is to make sure tinker doesn't destroy GPU memory that's already allocated.
    {
        const int32_t n_qm = 4;
        const int32_t qm_indices[8] { 5, 6, 7, 8, 1, 2, 3, 4 };
        const char* xyz = "../tinkerbox_interface/test_system/qmmm.xyz";
        initialize_tinker(qm_indices, n_qm, xyz);
        const int32_t n_mm = get_n_mm();
        const int32_t n_total = n_mm + n_qm;

        {
            int32_t* qm_atom_type = new int32_t[n_qm];
            get_qm_atomic_indices(qm_atom_type);
            for (int i_qm = 0; i_qm < n_qm; i_qm++)
                printf("QM atom %d, atomic type = %d\n", i_qm, qm_atom_type[i_qm]);
            delete[] qm_atom_type;
        }
        
        {
            double* qm_mass = new double[n_qm];
            get_qm_mass(qm_mass);
            for (int i_qm = 0; i_qm < n_qm; i_qm++)
                printf("QM atom %d, mass = %.5f\n", i_qm, qm_mass[i_qm]);
            delete[] qm_mass;
        }
        
        {
            double* mm_mass = new double[n_mm];
            get_mm_mass(mm_mass);
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, mass = %.5f\n", i_mm, mm_mass[i_mm]);
            delete[] mm_mass;
        }
        
        {
            double* qm_coords = new double[n_qm * 3];
            get_qm_xyz(qm_coords);
            for (int i_qm = 0; i_qm < n_qm; i_qm++)
                printf("QM atom %d, coords = (%.10f, %.10f, %.10f)\n", i_qm, qm_coords[i_qm * 3 + 0], qm_coords[i_qm * 3 + 1], qm_coords[i_qm * 3 + 2]);
            delete[] qm_coords;
        }
        
        {
            double* mm_coords = new double[n_mm * 3];
            get_mm_xyz(mm_coords);
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, coords = (%.10f, %.10f, %.10f)\n", i_mm, mm_coords[i_mm * 3 + 0], mm_coords[i_mm * 3 + 1], mm_coords[i_mm * 3 + 2]);
            delete[] mm_coords;
        }

        {
            double* charges = new double[n_mm];
            get_mm_charge(charges);
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, charge = %.10f\n", i_mm, charges[i_mm]);
            delete[] charges;
        }

        {
            double* dipoles = new double[n_mm * 3];
            get_mm_static_point_dipole(dipoles);
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, dipole = (%.10f, %.10f, %.10f)\n", i_mm, dipoles[i_mm * 3 + 0], dipoles[i_mm * 3 + 1], dipoles[i_mm * 3 + 2]);
            delete[] dipoles;
        }

        {
            double* polarizabilities = new double[n_mm];
            get_mm_polarizability(polarizabilities);
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, polarizability = %.10f\n", i_mm, polarizabilities[i_mm]);
            delete[] polarizabilities;
        }

        {
            double energy = get_energy_nonpolar_mm_contribution();
            printf("energy = %.10f\n", energy);
        }   

        // {
        //     double new_xyz[12] = {
        //         0.0, 0.0, 0.0,
        //         1.0, 0.0, 0.0,
        //         0.0, 3.0, 0.0,
        //         1.0, 3.0, 0.0,
        //     }; // Angstrom

        //     const double angstrom_to_bohr = 1.0 / 0.52917721067;
        //     for (int i = 0; i < n_mm * 3; i++)
        //         new_xyz[i] *= angstrom_to_bohr;

        //     set_mm_xyz(new_xyz);

        //     double energy = get_energy_nonpolar_mm_contribution();
        //     printf("new energy = %.10f\n", energy);
        // }

        {
            double* Ed = new double[n_mm * 3];
            double* Ep = new double[n_mm * 3];
            get_electric_field_mm_contribution(Ed, Ep);
            printf("Direct field:\n");
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, Ed = (%.10f, %.10f, %.10f)\n", i_mm, Ed[i_mm * 3 + 0], Ed[i_mm * 3 + 1], Ed[i_mm * 3 + 2]);
            printf("Polarization field:\n");
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, Ep = (%.10f, %.10f, %.10f)\n", i_mm, Ep[i_mm * 3 + 0], Ep[i_mm * 3 + 1], Ep[i_mm * 3 + 2]);
            delete[] Ed;
            delete[] Ep;
        }

        // {
        //     double* Ed = new double[n_mm * 3] {
        //         1.0, 0.0, 0.0,
        //         0.5, 0.0, 0.0,
        //     };
        //     double* Ep = new double[n_mm * 3] {
        //         0.0, 0.0, 0.0,
        //         0.1, 0.0, 0.0,
        //     };
        //     double* mu_d = new double[n_mm * 3];
        //     double* mu_p = new double[n_mm * 3];

        //     printf("Applying direct field:\n");
        //     for (int i_mm = 0; i_mm < n_mm; i_mm++)
        //         printf("MM atom %d, Ed = (%.10f, %.10f, %.10f)\n", i_mm, Ed[i_mm * 3 + 0], Ed[i_mm * 3 + 1], Ed[i_mm * 3 + 2]);
        //     printf("Applying polarization field:\n");
        //     for (int i_mm = 0; i_mm < n_mm; i_mm++)
        //         printf("MM atom %d, Ep = (%.10f, %.10f, %.10f)\n", i_mm, Ep[i_mm * 3 + 0], Ep[i_mm * 3 + 1], Ep[i_mm * 3 + 2]);

        //     evaluate_induced_dipole_from_total_electric_field(Ed, Ep, mu_d, mu_p);

        //     printf("Obtained direct dipole:\n");
        //     for (int i_mm = 0; i_mm < n_mm; i_mm++)
        //         printf("MM atom %d, mu_d = (%.10f, %.10f, %.10f)\n", i_mm, mu_d[i_mm * 3 + 0], mu_d[i_mm * 3 + 1], mu_d[i_mm * 3 + 2]);
        //     printf("Obtained polarization dipole:\n");
        //     for (int i_mm = 0; i_mm < n_mm; i_mm++)
        //         printf("MM atom %d, mu_p = (%.10f, %.10f, %.10f)\n", i_mm, mu_p[i_mm * 3 + 0], mu_p[i_mm * 3 + 1], mu_p[i_mm * 3 + 2]);
                
        //     delete[] Ed;
        //     delete[] Ep;
        //     delete[] mu_d;
        //     delete[] mu_p;
        // }

        {
            double* gradient = new double[n_total * 3];
            get_gradients_all_atoms_mm_contribution(gradient);

            for (int i = 0; i < n_total; i++)
                printf("Atom %d gradient %.10f, %.10f, %.10f\n", i, gradient[i*3+0], gradient[i*3+1], gradient[i*3+2]);
            delete[] gradient;
        }

        {
            const double* torque = new double[n_mm * 3] {
                1.0, 0.0, 0.0,
                1.0, 1.0, 1.0,
                0.0, 0.0, 0.0,
                0.0, 0.0, 0.0,
            };
            double* additional_gradient = new double[n_mm * 3] {
                100.0, 200.0, 300.0,
                100.0, 200.0, 300.0,
                100.0, 200.0, 300.0,
                100.0, 200.0, 300.0,
            };

            append_gradient_from_static_dipole_rotation(torque, additional_gradient);

            for (int i = 0; i < n_mm; i++)
                printf("MM atom %d gradient from rotation = %.10f, %.10f, %.10f\n", i, additional_gradient[i*3+0], additional_gradient[i*3+1], additional_gradient[i*3+2]);
            delete[] torque;
            delete[] additional_gradient;
        }

        {
            double* mu_d = new double[n_mm * 3];
            double* mu_p = new double[n_mm * 3];
            get_mm_induced_dipole(mu_d, mu_p);
            printf("Direct dipole:\n");
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, mu_d = (%.10f, %.10f, %.10f)\n", i_mm, mu_d[i_mm * 3 + 0], mu_d[i_mm * 3 + 1], mu_d[i_mm * 3 + 2]);
            printf("Polarization dipole:\n");
            for (int i_mm = 0; i_mm < n_mm; i_mm++)
                printf("MM atom %d, mu_p = (%.10f, %.10f, %.10f)\n", i_mm, mu_p[i_mm * 3 + 0], mu_p[i_mm * 3 + 1], mu_p[i_mm * 3 + 2]);
            delete[] mu_d;
            delete[] mu_p;
        }
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

    printf("\n\nUnrelated cuda result = %.2f (expected: 2077650.00)\n", C[M * N - 1]);

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
