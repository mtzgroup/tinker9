// Compile this file by "g++ ../tinkerbox_interface/test_tinkerbox_interface.cpp libtinkerbox.so -o test_tinkerbox_interface;"
// Run the executable by "LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH ./test_tinkerbox_interface;"

#include "tinkerbox.h"

int main()
{
	char* xyz = "~/test/tinker-banchmark/water_bulk.xyz";
	initialize_tinker(nullptr, 0, xyz);
	return 0;
}
