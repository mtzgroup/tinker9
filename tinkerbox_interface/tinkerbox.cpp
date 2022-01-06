#include <stdio.h>

#include "main_tinker9.h"

#include "tinkerbox.h"

extern "C" void test_print()
{
    printf("Henry: Test test\n");

    using namespace tinker;
    x_testgrad(0, NULL);
}
