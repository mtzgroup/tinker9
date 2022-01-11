// Compile this file by "icpc ../tinkerbox_interface/test_tinkerbox_interface.cpp libtinkerbox.so -std=c++11 -o test_tinkerbox_interface;"
// Run the executable by "LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH ./test_tinkerbox_interface;"

#include "tinkerbox.h"

#include <stdio.h>

int main()
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

	return 0;
}
