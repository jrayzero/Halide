#include "Halide.h"
#include <stdio.h>

#define CV_DESCALE(x,n) (((x) + (1 << ((n)-1))) >> (n))

#define NOCLAMP 0

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
        Buffer<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(max(x, y), y) >> 2);
        f.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
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
        f.compile_to_tiramisu("test_reduction_operator.cpp", "test_reduction_operator");
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

        C.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Buffer<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        f.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
    }

    if (0) {
        Func f("f");
        Var x("x"), y("y");

        Func in("in");
        in(x, y) = x + y;
        Buffer<int> input = in.realize(100, 50);

        f(x, y) = cast(Float(32), input(x, y) >> 2);
        f.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
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
        Buffer<int> f_img(Int(32), 100, 50);
        f.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
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
        Buffer<int> f_img(Int(32), 100, 50);
        f.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
    }

    if (0) {
        const int N = 100;

        Buffer<int> A(N, N), B(N, N);
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

        C.compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
    }

    // Fusion benchmark
    if (0) {
        //ImageParam in{UInt(8), 3, "input"};
        Buffer<uint8_t> in(2124, 3540, 3);

        Var x("x"), y("y"), c("c");
        Func f("f"), g("g");

        f(x, y, c) = cast<uint8_t>(255 - in(x, y, c));
        g(x, y, c) = cast<uint8_t>(100 + in(x, y, c));

        Buffer<uint8_t> f_im(2124, 3540, 3);
        Buffer<uint8_t> g_im(2124, 3540, 3);
        Pipeline({f, g}).realize({f_im, g_im});

        //Pipeline({f, g}).compile_to_tiramisu("fusion_tiramisu.cpp", "fusion_tiramisu");
    }

    // Filter 2D
    if (0) {
        ImageParam in(UInt(8), 3, "input");
        ImageParam kernel(Float(32), 2, "kernel");

        /*Buffer<uint8_t> in(2124, 3540, 3);
        Buffer<float> kernel(3, 3);
        kernel(0,0) = 0; kernel(0,1) = 1.0f/5; kernel(0,2) = 0;
        kernel(1,0) = 1.0f/5; kernel(1,1) = 1.0f/5; kernel(1,2) = 1.0f/5;
        kernel(2,0) = 0; kernel(2,1) = 1; kernel(2,2) = 0;*/

        Func filter2D("filter2D");
        Var x("x"), y("y"), c("c");

        Expr e = 0.0f;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < 3; i++) {
                e += cast<float>(in(x + i, y + j, c)) * kernel(i, j);
            }
        }

        filter2D(x, y, c) = cast<uint8_t>(e);

        filter2D.parallel(y);//.vectorize(x, 8);

        //filter2D.realize(2116, 3532, 3);

        filter2D.compile_to_tiramisu("filter2D_tiramisu.cpp", "filter2D_tiramisu");
    }

    // Gaussian 5x5 benchmark
    if (0) {
        ImageParam in(UInt(8), 3, "input");
        ImageParam kernelX(Float(32), 1, "kernelx");
        ImageParam kernelY(Float(32), 1, "kernely");

        /*Buffer<uint8_t> in(2124, 3540, 3);
        Buffer<float> kernelX(5);
        Buffer<float> kernelY(5);*/

        Func gaussian("gaussian"), gaussian_x("gaussian_x");
        Var x("x"), y("y"), c("c");

        Expr e = 0.0f;
        for (int i = 0; i < 5; ++i) {
            e += in(x + i, y, c) * kernelX(i);
        }
        gaussian_x(x, y, c) = e;

        Expr f = 0.0f;
        for (int j = 0; j < 5; ++j) {
            f += gaussian_x(x, y + j, c) * kernelY(j);
        }
        gaussian(x, y, c) = cast<uint8_t>(f);

        /*Var x_inner, y_inner, x_outer, y_outer, tile_index;
        gaussian.tile(x, y, x_outer, y_outer, x_inner, y_inner, 4, 4)
                .fuse(x_outer, y_outer, tile_index).compute_root().parallel(tile_index);
        gaussian_x.compute_at(gaussian, y_inner);*/

        gaussian_x.compute_root();

        //gaussian.realize(2116, 3532, 3);
        gaussian.compile_to_tiramisu("gaussian_tiramisu.cpp", "gaussian_tiramisu");
    }

    // Blur x-y
    if (0) {
        ImageParam input(UInt(8), 3);
        //Buffer<uint8_t> input(2124, 3540, 3);
        Func blur_x("blur_x"), blur_y("blur_y");
        Var x("x"), y("y"), c("c"), xi("xi"), yi("yi");

        // The algorithm
        blur_x(x, y, c) = (input(x, y, c) + input(x+1, y, c) + input(x+2, y, c))/3;
        blur_y(x, y, c) = (blur_x(x, y, c) + blur_x(x, y+1, c) + blur_x(x, y+2, c))/3;

        // How to schedule it
        blur_y.split(y, y, yi, 8).parallel(y).vectorize(x, 8).parallel(c);
        blur_x.store_at(blur_y, y).compute_at(blur_y, yi).vectorize(x, 8).parallel(c);

        //blur_y.realize(2116, 3532, 3);
        blur_y.compile_to_tiramisu("blurxy_tiramisu.cpp", "blurxy_tiramisu");
    }

    // RGB to YUV420
    if (1) {
        ImageParam rgb(UInt(8), 3);
        //Buffer<uint8_t> rgb(2124, 3540, 3);

        Var x("x"), y("y"), z("z");
        Func y_part("y_part"), u_part("u_part"), v_part("v_part");

        y_part(x, y) = cast<uint8_t>(((66 * cast<int>(rgb(x, y, 0)) + 129 * cast<int>(rgb(x, y, 1)) +  25 * cast<int>(rgb(x, y, 2)) + 128) % 256) +  16);
        u_part(x, y) = cast<uint8_t>((( -38 * cast<int>(rgb(2*x, 2*y, 0)) -  cast<int>(74 * rgb(2*x, 2*y, 1)) + 112 * cast<int>(rgb(2*x, 2*y, 2) + 128)) % 256) + 128);
        v_part(x, y) = cast<uint8_t>((( 112 * cast<int>(rgb(2*x, 2*y, 0)) -  cast<int>(94 * rgb(2*x, 2*y, 1)) -  18 * cast<int>(rgb(2*x, 2*y, 2) + 128)) % 256) + 128);

        y_part.parallel(y).vectorize(x, 8);
        u_part.parallel(y).vectorize(x, 8);
        v_part.parallel(y).vectorize(x, 8);

        /*const int size = 100;
        Buffer<uint8_t> y_im(size, size), u_im(size/2, size/2), v_im(size/2, size/2);
        Pipeline({y_part, u_part, v_part}).realize({y_im, u_im, v_im});*/

        Pipeline({y_part, u_part, v_part}).compile_to_tiramisu("rgbyuv420_tiramisu.cpp", "rgbyuv420_tiramisu");
    }

    // CVT Color benchmark
    if (0) {
        ImageParam in{UInt(8), 3, "input"};
        //Buffer<uint8_t> in(2116, 3538, 3);

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

        RGB2Gray.parallel(y).vectorize(x, 8);

        //RGB2Gray.realize(2116, 3538);

        RGB2Gray.compile_to_tiramisu("cvtcolor.cpp", "cvtcolor");
    }

    // Rec-filter
    if (0) {
        //ImageParam in(UInt(8), 3, "input"};
        Buffer<uint8_t> in(2116, 3538, 3);

        Param<float> a0{"a0"};
        Param<float> a1{"a1"};
        Param<float> a2{"a2"};

        Func rec_filter{"rec_filter"};
        Var x("x"), y("y"), c("c");
        RDom r(2, in.width()-3, 0, in.height()-1);

        rec_filter(x, y, c) = in(x, y, c);
        rec_filter(r.x, r.y, c) = cast<uint8_t>(a0*rec_filter(r.x, r.y, c) + a1*rec_filter(r.x-1, r.y, c) + a2*rec_filter(r.x-2, r.y, c));

        rec_filter.parallel(y);//.vectorize(x, 8);
        rec_filter.update(0).parallel(r.y);

        rec_filter.realize(2116, 3538, 3);

        rec_filter.compile_to_tiramisu("recfilter_tiramisu.cpp", "recfilter_tiramisu");
    }

    if (0) {
        ImageParam in(Float(32), 2, "input");
        Param<float> alpha;
        Param<float> beta;

        //Buffer<float> in(10000, 10000);
        //float alpha = 0.3;
        //float beta = 0.4;

        Func heat2d("heat2d");
        Var x("x"), y("y");

        RDom r(1, in.width()-2, 1, in.height()-2);
        heat2d(x, y) = 0.0f;
        heat2d(r.x, r.y) = alpha * in(r.x, r.y) +
                           beta * (in(r.x+1, r.y) + in(r.x-1, r.y) + in(r.x, r.y+1) + in(r.x, r.y-1));

        heat2d.parallel(y);//.vectorize(x, 8);
        heat2d.update().parallel(r.y);//.vectorize(r.x, 8);

        //heat2d.realize(10000, 10000);

        heat2d.compile_to_tiramisu("heat2d_tiramisu.cpp", "heat2d_tiramisu");
    }

    if (0) {
        /*ImageParam in(Float(32), 2, "input");
        Param<float> alpha;
        Param<float> beta;
        Param<int> t;*/

        Buffer<float> in(10000, 10000);
        float alpha = 0.3;
        float beta = 0.4;
        int t = 100;

        Func heat3d{"heat3d"};
        Var x, y, z;
        Var k;

        RDom time{1, t, "t"};

        heat3d(x, y, z, k) = undef<float>();
        heat3d(x, y, z, 0) = alpha * in(x, y, z) +
                             beta * (in(x+1, y, z) + in(x-1, y, z)+
                                     in(x, y+1, z) + in(x, y-1, z) +
                                     in(x, y, z+1) + in(x, y, z-1));


        heat3d(x, y, z, time.x) = alpha * heat3d(x, y, z, time.x-1) +
                                  beta * (heat3d(x+1, y, z, time.x-1) + heat3d(x-1, y, z, time.x-1)+
                                          heat3d(x, y+1, z, time.x-1) + heat3d(x, y-1, z, time.x-1) +
                                          heat3d(x, y, z+1, time.x-1) + heat3d(x, y, z-1, time.x-1));

        heat3d.parallel(z).vectorize(x, 8);

        heat3d.realize(10000, 10000);

        //heat3d.compile_to_tiramisu("heat3d_tiramisu.cpp", "heat3d_tiramisu");
    }

    if (0) {
        ImageParam in(Float(32), 2, "input");
        Param<float> alpha;
        Param<float> beta;

        //Buffer<float> in(10000, 10000);
        //float alpha = 0.3;
        //float beta = 0.4;

        Func divergence2d("divergence2d");
        Var x("x"), y("y");

        RDom r(1, in.width()-2, 1, in.height()-2);
        divergence2d(x, y) = 0.0f;
        divergence2d(r.x, r.y) = alpha * (in(r.x+1, r.y) + in(r.x-1, r.y)) +
                                 beta  * (in(r.x, r.y+1) + in(r.x, r.y-1));

        divergence2d.parallel(y);//.vectorize(x, 8);
        divergence2d.update().parallel(r.y);//.vectorize(r.x, 8);

        //divergence2d.realize(10000, 10000);
        divergence2d.compile_to_tiramisu("divergence2d_tiramisu.cpp", "divergence2d_tiramisu");
    }

    if (0) {
        ImageParam in(Float(32), 2, "input");
        Param<float> alpha;
        Param<float> beta;

        //Buffer<float> in(10000, 10000);
        //float alpha = 0.3;
        //float beta = 0.4;

        Func divergence2d("divergence2d");
        Var x("x"), y("y");

        RDom r(1, in.width()-2, 1, in.height()-2);
        divergence2d(x, y) = 0.0f;
        divergence2d(r.x, r.y) = alpha * (in(r.x+1, r.y) + in(r.x-1, r.y)) +
                                 beta  * (in(r.x, r.y+1) + in(r.x, r.y-1));

        divergence2d.parallel(y);//.vectorize(x, 8);
        divergence2d.update().parallel(r.y);//.vectorize(r.x, 8);

        //divergence2d.realize(10000, 10000);
        divergence2d.compile_to_tiramisu("divergence2d_tiramisu.cpp", "divergence2d_tiramisu");
    }

    if (0) {
        ImageParam in{UInt(8), 2, "input"};
        Param<float> a00;
        Param<float> a01;
        Param<float> a10;
        Param<float> a11;
        Param<float> b00;
        Param<float> b10;

        /*Buffer<uint8_t> in(100, 100);
        float a00 = 0.1;
        float a01 = 0.1;
        float a10 = 0.1;
        float a11 = 0.1;
        float b00 = 0.1;
        float b10 = 0.1;*/

        Func affine("affine");
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

    #ifdef NOCLAMP
        coord_00_r = clamp(coord_00_r, 0, src_rows);
        coord_00_c = clamp(coord_00_c, 0, src_cols);
        coord_01_r = clamp(coord_01_r, 0, src_rows);
        coord_01_c = clamp(coord_01_c, 0, src_cols);
        coord_10_r = clamp(coord_10_r, 0, src_rows);
        coord_10_c = clamp(coord_10_c, 0, src_cols);
        coord_11_r = clamp(coord_11_r, 0, src_rows);
        coord_11_c = clamp(coord_11_c, 0, src_cols);
    #else
        coord_00_r = coord_00_r;
        coord_00_c = coord_00_c;
        coord_01_r = coord_01_r;
        coord_01_c = coord_01_c;
        coord_10_r = coord_10_r;
        coord_10_c = coord_10_c;
        coord_11_r = coord_11_r;
        coord_11_c = coord_11_c;
    #endif

        Expr A00 = in(coord_00_r, coord_00_c);
        Expr A10 = in(coord_10_r, coord_10_c);
        Expr A01 = in(coord_01_r, coord_01_c);
        Expr A11 = in(coord_11_r, coord_11_c);

        affine(x, y) = mixf( mixf(A00, A10, r), mixf(A01, A11, r), c);

        affine.parallel(y).vectorize(x, 8);

        //affine.realize(100, 100);
        affine.compile_to_tiramisu("affine_tiramisu.cpp", "affine_tiramisu");
    }

    if (0) {
        ImageParam in{UInt(8), 2, "input"};
        Param<int> resampled_rows;
        Param<int> resampled_cols;

        /*Buffer<uint8_t> in(100, 100);
        int resampled_rows = 5;
        int resampled_cols = 5;*/

        Func resampled{"resampled"};
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

    #ifdef NOCLAMP
        Expr coord_00_r = cast<int>(floor(o_r));
        Expr coord_00_c = cast<int>(floor(o_c));
    #else
        Expr coord_00_r = clamp( cast<int>(floor(o_r)), 0, o_h - 1 );
        Expr coord_00_c = clamp( cast<int>(floor(o_c)), 0, o_w - 1 );
    #endif

        Expr coord_01_r = coord_00_r;

    #ifdef NOCLAMP
        Expr coord_01_c = coord_00_c + 1;
        Expr coord_10_r = coord_00_r + 1;
    #else
        Expr coord_01_c = clamp( coord_00_c + 1, 0, o_w - 1 );
        Expr coord_10_r = clamp( coord_00_r + 1, 0, o_h - 1 );
    #endif

        Expr coord_10_c = coord_00_c;

    #ifdef NOCLAMP
        Expr coord_11_r = coord_00_r + 1;
        Expr coord_11_c = coord_00_c + 1;
    #else
        Expr coord_11_r = clamp( coord_00_r + 1, 0, o_h - 1 );
        Expr coord_11_c = clamp( coord_00_c + 1, 0, o_w - 1 );
    #endif

        Expr A00 = in(coord_00_r, coord_00_c);
        Expr A10 = in(coord_10_r, coord_10_c);
        Expr A01 = in(coord_01_r, coord_01_c);
        Expr A11 = in(coord_11_r, coord_11_c);

        resampled(x, y) = mixf( mixf(A00, A10, r), mixf(A01, A11, r), c);

        resampled.parallel(y).vectorize(x, 8);

        //resampled.realize(100, 100);
        resampled.compile_to_tiramisu("resize_tiramisu.cpp", "resize_tiramisu");
    }

    printf("Success!\n");
    return 0;
}
