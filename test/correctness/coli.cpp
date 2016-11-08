#include "Halide.h"
#include <stdio.h>

#define CV_DESCALE(x,n) (((x) + (1 << ((n)-1))) >> (n))

using namespace Halide;
using namespace Halide::Internal;

Expr mixf(Expr x, Expr y, Expr a) {
    return x * (1.0f-a) + y * a;
}

int main(int argc, char **argv) {
    if (0) { // Not working
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(max(x, y), y) >> 2);
        f.compile_to_coli("fusion_coli.cpp", {}, "fusion_coli");
    }

    if (0) { // With reduction
        Func f("f");
        Var x("x"), y("y"), z("z");

        Func input("input");
        input(x, y, z) = cast<uint8_t>(13);
        input.compute_root();

        RDom r(0, 10);
        f(x, y) = cast<uint8_t>(0);
        f(x, y) += input(x, y, r);
        f.compile_to_coli("test_reduction_operator.cpp", {}, "test_reduction_operator");
    }

    if (0) {
        const int N = 100;

        ImageParam A(Int(8), 2, "A");
        ImageParam B(Int(8), 2, "B");

        Func C("C");
        Var x("x"), y("y");

        RDom r(0, N);

        C(x, y) = cast<int8_t>(0);
        C(x, y) += A(x, r) * B(r, y);

        C.compile_to_coli("fusion_coli.cpp", {A, B}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        f.compile_to_coli("fusion_coli.cpp", {}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        f.compile_to_coli("fusion_coli.cpp", {}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = 13;
        in.compute_root();
        in.parallel(x);
        in.reorder(y, x);

        f(x, y) = cast(Float(32), in(x, y) >> 2);
        Image<int> f_img(Int(32), 100, 50);
        f.compile_to_coli("fusion_coli.cpp", {}, "fusion_coli");
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

        C.compile_to_coli("fusion_coli.cpp", {}, "fusion_coli");
    }

    // CVT Color benchmark
    if (0) {
        ImageParam in{UInt(8), 3, "input"};

        Func RGB2Gray{"RGB2Gray"};
        Var x,y,c;

        const Expr yuv_shift = cast<uint32_t>(14);
        const Expr R2Y = cast<uint32_t>(4899);
        const Expr G2Y = cast<uint32_t>(9617);
        const Expr B2Y = cast<uint32_t>(1868);


        RGB2Gray(x, y) = cast<uint8_t>(CV_DESCALE( (in(x, y, 2) * B2Y
                                  + in(x, y, 1) * G2Y
                                  + in(x, y, 0) * R2Y),
                                  yuv_shift));

        RGB2Gray.parallel(y);//.vectorize(x, 8);

        RGB2Gray.compile_to_coli("cvtcolor.cpp", {in}, "cvtcolor");
    }

    if (0) {
        //ImageParam in{UInt(8), 3, "input"};
        Image<uint8_t> in(1024, 512, 3);

        Func RGB2Gray{"RGB2Gray"};
        Var x,y,c;

        RGB2Gray(x, y) = cast<uint8_t>(in(x, y, 2) + in(x, y, 1) + in(x, y, 0));

        RGB2Gray.parallel(y);//.vectorize(x, 8);

        RGB2Gray.realize(1024, 512);

        //RGB2Gray.compile_to_coli("cvtcolor.cpp", {in}, "cvtcolor");
    }

    // Filter 2D nordom benchmark
    if (0) {
        const int RADIUS = 3;

        ImageParam in{Float(32), 2, "input"};
        // kernel is 2*radius x 2*radius
        ImageParam kernel{Float(32), 2, "kernel"};

        Func filter2D_nordom{"filter2D_nordom"};
        Var x("x"), y("y");

        Expr e = 0.0f;

        for (int i=-RADIUS; i<RADIUS; i++) {
            for (int j=-RADIUS; j<RADIUS; j++)  {
                e += in(x+RADIUS+i, y+RADIUS+j) * kernel(RADIUS+i, RADIUS+j);
            }
        }

        filter2D_nordom(x, y) = e;

        filter2D_nordom.parallel(y);//.vectorize(x, 8);

        filter2D_nordom.compile_to_coli("filter2Dnordom.cpp", {in, kernel}, "filter2Dnordom");
    }

    // Gaussian 3x3 benchmark
    if (1) {
        const int kernelX_length = 7;

        ImageParam in{Float(32), 2, "input"};
        ImageParam kernelX{Float(32), 1, "kernelx"};
        ImageParam kernelY{Float(32), 1, "kernely"};

        Func gaussian("gaussian");
        Func gaussian_x("gaussian_x");
        Var x("x"), y("y");

        Expr e,f;
        e = 0.0f;
        for (int i=0; i<kernelX_length; i++) {
            e += in(x+i,y) * kernelX(i);
        }
        gaussian_x(x, y) = e;

        f = 0.0f;
        for (int i=0; i<kernelX_length; i++) {
            f += gaussian_x(x, y+i) * kernelY(i);
        }

        gaussian(x, y) = f;

        gaussian_x.compute_root();
        gaussian.compute_root();

        gaussian.compile_to_coli("gaussian3x3.cpp", {in, kernelX, kernelY}, "gaussian3x3");
    }

    printf("Success!\n");
    return 0;
}
