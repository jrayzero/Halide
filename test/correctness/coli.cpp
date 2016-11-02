#include "Halide.h"
#include <stdio.h>

using namespace Halide;
using namespace Halide::Internal;

int main(int argc, char **argv) {
    if (1) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 100);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        Image<int> f_img(Int(32), 100, 100);
        f.compile_to_coli({f_img}, "fusion_coli", {}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = 13;
        in.compute_root();

        f(x, y) = cast(Float(32), in(x, y) >> 2);
        Image<int> f_img(Int(32), 100, 100);
        f.compile_to_coli({f_img}, "fusion_coli", {}, "fusion_coli");
    }

    if (0) {
        const int N = 100;

        Image<int> A(N, N), B(N, N);
        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                A(x, y) = rand() & 0xffffffff;
                B(x, y) = rand() & 0xffffffff;
            }
        }

        Func C("C");
        Var x("x"), y("y");

        RDom r(0, N);
        C(x,y) += A(x,r) * B(r,y);

        //C.realize(N, N);
        Image<int> C_img(Int(32), N, N);
        C.compile_to_coli({C_img}, "fusion_coli", {}, "fusion_coli");
    }

    printf("Success!\n");
    return 0;
}
