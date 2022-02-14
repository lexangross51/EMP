#include <iostream>
#include "testing_module.h"

int main()
{
    testing_module test;

    test.set_functions();
    test.run_tests();

    return 0;
}