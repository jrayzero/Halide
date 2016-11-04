#include "Halide.h"
#include <stdio.h>

#define CV_DESCALE(x,n) (((x) + (1 << ((n)-1))) >> (n))

using namespace Halide;
using namespace Halide::Internal;

Expr mixf(Expr x, Expr y, Expr a) {
    return x * (1.0f-a) + y * a;
}

int main(int argc, char **argv) {
    if (1) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(max(x, y), y) >> 2);
        Image<int> f_img(Int(32), 100, 50);
        f.compile_to_coli({f_img}, "fusion_coli", {}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        Image<int> f_img(Int(32), 100, 50);
        f.compile_to_coli({f_img}, "fusion_coli", {}, "fusion_coli");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Image<int> input = in.realize(50, 100);

        RDom r(0, 100);
        f(x) += input(x, r);
        Image<int> f_img(Int(32), 100);
        f.compile_to_coli({f_img}, "fusion_coli", {}, "fusion_coli");
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

    // CVT Color benchmark
    if (0) {
        const int N = 100;

        Image<int> in(N, N, 3);
        for (int c = 0; c < 3; c++) {
            for (int y = 0; y < N; y++) {
                for (int x = 0; x < N; x++) {
                    in(x, y, c) = rand() & 0xffffffff;
                }
            }
        }

        Func RGB2Gray{"RGB2Gray"};
        Var x("x"), y("y"), c("c");

        const Expr yuv_shift = cast<uint32_t>(14);
        const Expr R2Y = cast<uint32_t>(4899);
        const Expr G2Y = cast<uint32_t>(9617);
        const Expr B2Y = cast<uint32_t>(1868);


        RGB2Gray(x, y) = cast<uint8_t>(CV_DESCALE( (in(x, y, 2) * B2Y
                                  + in(x, y, 1) * G2Y
                                  + in(x, y, 0) * R2Y),
                                  yuv_shift));

        Image<uint8_t> RGB2Gray_img(UInt(8), N, N);
        RGB2Gray.compile_to_coli({RGB2Gray_img}, "fusion_coli", {}, "fusion_coli");
    }

    // Resize benchmark
    if (0) {
        const int WIDTH = 256;
        const int HEIGHT = 128;

        Image<uint8_t> in(WIDTH, HEIGHT);

        int resampled_rows = 32;
        int resampled_cols = 64;

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

        Image<uint8_t> resampled_img(UInt(8), WIDTH, HEIGHT);
        resampled.compile_to_coli({resampled_img}, "fusion_coli", {}, "fusion_coli");
    }

    if (0) {
        const int WIDTH = 256;
        const int HEIGHT = 128;

        Image<uint8_t> in(WIDTH, HEIGHT);

        float a00 = 0.1;
        float a01 = 0.1;
        float a10 = 0.1;
        float a11 = 0.1;
        float b00 = 0.1;
        float b10 = 0.1;

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

        Image<uint8_t> affine_img(UInt(8), WIDTH, HEIGHT);
        affine.compile_to_coli({affine_img}, "fusion_coli", {}, "fusion_coli");
    }

    printf("Success!\n");
    return 0;
}
