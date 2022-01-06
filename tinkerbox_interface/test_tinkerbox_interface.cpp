// Compile this file by "g++ ../tinkerbox_interface/test_tinkerbox_interface.cpp libtinkerbox.so -o test_tinkerbox_interface;"
// Run the executable by "LD_LIBRARY_PATH=$PWD:$LD_LIBRARY_PATH ./test_tinkerbox_interface;"

extern "C" void test_print();

int main()
{
	test_print();
	return 0;
}
