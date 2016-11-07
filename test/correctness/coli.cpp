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
        f.compile_to_coli("fusion_coli", {}, "fusion_coli");
    }

    if (0) { // With reduction
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(50, 100);

        RDom r(0, 100);
        f(x) += input(x, r);
        f.compile_to_coli("/Users/psuriana/COLi/Halide/fusion_coli.cpp", {}, "fusion_coli");
    }

    if (1) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        f.compile_to_coli("/Users/psuriana/COLi/Halide/fusion_coli.cpp", {}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        f.compile_to_coli("fusion_coli", {}, "fusion_coli");
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
        f.compile_to_coli("fusion_coli", {}, "fusion_coli");
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

        C.compile_to_coli("fusion_coli", {}, "fusion_coli");
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

        RGB2Gray.compile_to_coli("cvtcolor", {in}, "cvtcolor");
    }

    // Resize benchmark
    if (0) {
        ImageParam in(UInt(8), 2, "input");
        Param<int> resampled_rows;
        Param<int> resampled_cols;

        Func resampled("resampled");
        Var x("x"), y("y");

        // Translating this algorithm as close as possible
        Expr o_h = in.height();
        Expr o_w = in.width();
        Expr n_h = resampled_rows;
        Expr n_w = resampled_cols;

        Expr n_r = y;
        Expr n_c = x;

        Expr o_r = (n_r + 0.5f) * (o_h) / (n_h) - 0.5f;
        Expr o_c = (n_c + 0.5f) * (o_w) / (n_w) - 0.5f;

        Expr r = o_r - floor(o_r);
        Expr c = o_c - floor(o_c);

        Expr coord_00_r = clamp( cast<int>(floor(o_r)), 0, o_h - 1 );
        Expr coord_00_c = clamp( cast<int>(floor(o_c)), 0, o_w - 1 );

        Expr coord_01_r = coord_00_r;
        Expr coord_01_c = clamp( coord_00_c + 1, 0, o_w - 1 );

        Expr coord_10_r = clamp( coord_00_r + 1, 0, o_h - 1 );
        Expr coord_10_c = coord_00_c;

        Expr coord_11_r = clamp( coord_00_r + 1, 0, o_h - 1 );
        Expr coord_11_c = clamp( coord_00_c + 1, 0, o_w - 1 );

        Expr A00 = in(coord_00_r, coord_00_c);
        Expr A10 = in(coord_10_r, coord_10_c);
        Expr A01 = in(coord_01_r, coord_01_c);
        Expr A11 = in(coord_11_r, coord_11_c);

        resampled(x, y) = mixf( mixf(A00, A10, r), mixf(A01, A11, r), c);

        //TODO(psuriana): vectorize is not currently working
        resampled.parallel(y);//.vectorize(x, 8);

        resampled.compile_to_coli("fusion_coli", {in, resampled_rows, resampled_cols}, "fusion_coli");
    }

    // Wrap affine benchmark
    if (0) {
        ImageParam in{UInt(8), 2, "input"};
        Param<float> a00;
        Param<float> a01;
        Param<float> a10;
        Param<float> a11;
        Param<float> b00;
        Param<float> b10;

        Func affine{"affine"};
        Var x("x"), y("y");

        // Translating this algorithm as close as possible
        Expr o_r = a11 * y + a10 * x + b00;
        Expr o_c = a01 * y + a00 * x + b10;

        Expr r = o_r - floor(o_r);
        Expr c = o_c - floor(o_c);

        Expr coord_00_r = cast<int>(floor(o_r));
        Expr coord_00_c = cast<int>(floor(o_c));
        Expr coord_01_r = coord_00_r;
        Expr coord_01_c = coord_00_c + 1;
        Expr coord_10_r = coord_00_r + 1;
        Expr coord_10_c = coord_00_c;
        Expr coord_11_r = coord_00_r + 1;
        Expr coord_11_c = coord_00_c + 1;

        Expr src_rows = in.height();
        Expr src_cols = in.width();

        coord_00_r = clamp(coord_00_r, 0, src_rows);
        coord_00_c = clamp(coord_00_c, 0, src_cols);
        coord_01_r = clamp(coord_01_r, 0, src_rows);
        coord_01_c = clamp(coord_01_c, 0, src_cols);
        coord_10_r = clamp(coord_10_r, 0, src_rows);
        coord_10_c = clamp(coord_10_c, 0, src_cols);
        coord_11_r = clamp(coord_11_r, 0, src_rows);
        coord_11_c = clamp(coord_11_c, 0, src_cols);

        Expr A00 = in(coord_00_r, coord_00_c);
        Expr A10 = in(coord_10_r, coord_10_c);
        Expr A01 = in(coord_01_r, coord_01_c);
        Expr A11 = in(coord_11_r, coord_11_c);

        affine(x, y) = mixf( mixf(A00, A10, r), mixf(A01, A11, r), c);

        //TODO(psuriana): vectorize is not currently working
        affine.parallel(y);//.vectorize(x, 8);

        affine.compile_to_coli("fusion_coli", {in, a00, a01, a10, a11, b00, b10}, "fusion_coli");
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

        filter2D_nordom.compile_to_coli("filter2Dnordom", {in, kernel}, "filter2Dnordom");
    }

    // Gaussian 3x3 benchmark
    if (0) {
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

        //gaussian.realize(WIDTH, HEIGHT);

        gaussian.compile_to_coli("gaussian3x3", {in, kernelX, kernelY}, "gaussian3x3");
    }

    printf("Success!\n");
    return 0;
}
