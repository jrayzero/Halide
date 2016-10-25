#include "Halide.h"
#include <stdio.h>

using namespace Halide;
using namespace Halide::Internal;

int main(int argc, char **argv) {
    Func f("f");
    Var x("x"), y("y");

    Func in("in");
    in(x, y) = x + y;
    Image<int> input = in.realize(100, 100);

    f(x, y) = input(x, y);
    Image<int> f_img(Int(32), 100, 100);
    f.compile_to_coli({f_img}, "fusion_coli", {}, "fusion_coli");

    printf("Success!\n");
    return 0;
}
