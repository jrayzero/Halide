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
        Var x("x"), y("y"), c("c");

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
        Var x("x"), y("y"), c("c");

        RGB2Gray(x, y) = cast<uint8_t>(in(x, y, 2) + in(x, y, 1) + in(x, y, 0));

        RGB2Gray.parallel(y);//.vectorize(x, 8);

        RGB2Gray.realize(1024, 512);

        //RGB2Gray.compile_to_coli("cvtcolor.cpp", {in}, "cvtcolor");
    }

    if (0) {
        //ImageParam in{UInt(8), 3, "input"};
        Image<uint8_t> in(2124, 3540, 3);

        Var x("x"), y("y"), c("c");
        Func f("f"), g("g");

        f(x, y, c) = cast<uint8_t>(255 - in(x, y, c));
        g(x, y, c) = cast<uint8_t>(100 + in(x, y, c));

        Image<uint8_t> f_im(2124, 3540, 3);
        Image<uint8_t> g_im(2124, 3540, 3);
        Pipeline({f, g}).realize({f_im, g_im});

        //Pipeline({f, g}).compile_to_coli("fusion_coli.cpp", {in}, "fusion_coli");
    }

    if (0) {
        ImageParam in(UInt(8), 3, "input");
        ImageParam kernel(UInt(8), 2, "kernel");

        //Image<uint8_t> in(2124, 3540, 3);
        //Image<uint8_t> kernel(3, 3);

        Func filter2D("filter2D");
        Var x("x"), y("y"), c("c");
        RDom r(0, kernel.width(), 0, kernel.height());

        filter2D(x, y, c) += cast<uint8_t>(in(x + r.x, y + r.y, c));

        filter2D.parallel(y);//.vectorize(x, 8);

        //filter2D.realize(2116, 3532, 3);

        filter2D.compile_to_coli("filter2D_coli.cpp", {in, kernel}, "filter2D_coli");
    }

    if (0) {
        ImageParam in(UInt(8), 3, "input");
        ImageParam kernel(Float(32), 2, "kernel");

        /*Image<uint8_t> in(2124, 3540, 3);
        Image<float> kernel(3, 3);
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

        filter2D.compile_to_coli("filter2D_coli.cpp", {in, kernel}, "filter2D_coli");
    }

    if (0) {
        ImageParam in(UInt(8), 3, "input");
        ImageParam kernelX(Float(32), 1, "kernelx");
        ImageParam kernelY(Float(32), 1, "kernely");

        /*Image<uint8_t> in(2124, 3540, 3);
        Image<float> kernelX(5);
        Image<float> kernelY(5);*/

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
        gaussian.compile_to_coli("gaussian_coli.cpp", {in, kernelX, kernelY}, "gaussian_coli");
    }

    if (0) {
        ImageParam input(UInt(8), 3);
        Func blur_x("blur_x"), blur_y("blur_y");
        Var x("x"), y("y"), c("c"), xi("xi"), yi("yi");

        // The algorithm
        //blur_x(x, y, c) = (input(x, y, c) + input(x+1, y, c) + input(x+2, y, c))/3;
        //blur_y(x, y, c) = (blur_x(x, y, c) + blur_x(x, y+1, c) + blur_x(x, y+2, c))/3;

        blur_x(x, y, c) = (input(x, y, c) + input(x+1, y, c) + input(x+2, y, c))/3;
        blur_y(x, y, c) = (blur_x(x, y, c) + blur_x(x, y+1, c) + blur_x(x, y+2, c))/3;

        blur_y.compute_root().store_root().parallel(y).parallel(c);
        blur_x.compute_root().store_root().parallel(y).parallel(c);

        //blur_y.realize(2116, 3532, 3);
        blur_y.compile_to_coli("blurxy_coli.cpp", {input}, "blurxy_coli");
    }

    if (1) {
        ImageParam rgb(Int(16), 3);
        //Image<int16_t> rgb(2124, 3540, 3);

        Var x("x"), y("y");
        Func y_part("y_part"), u_part("u_part"), v_part("v_part");

        y_part(x, y) = cast<uint8_t>(((66 * rgb(x, y, 0) + 129 * rgb(x, y, 1) +  25 * rgb(x, y, 2) + 128) >> 8) +  16);
        u_part(x, y) = cast<uint8_t>((( -38 * rgb(2*x, 2*y, 0) -  74 * rgb(2*x, 2*y, 1) + 112 * rgb(2*x, 2*y, 2) + 128) >> 8) + 128);
        v_part(x, y) = cast<uint8_t>((( 112 * rgb(2*x, 2*y, 0) -  94 * rgb(2*x, 2*y, 1) -  18 * rgb(2*x, 2*y, 2) + 128) >> 8) + 128);

        //u_part.compute_with(y_part, y);
        //v_part.compute_with(u_part, y);

        /*const int size = 100;
        Image<uint8_t> y_im(size, size), u_im(size/2, size/2), v_im(size/2, size/2);
        Pipeline({y_part, u_part, v_part}).realize({y_im, u_im, v_im});*/

        Pipeline({y_part, u_part, v_part}).compile_to_coli("rgbyuv420.cpp", {rgb}, "rgbyuv420");
    }

    if (0) {
        //ImageParam in{UInt(8), 3, "input"};
        Image<uint8_t> in(2116, 3538, 3);

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

        RGB2Gray.parallel(y);//.vectorize(x, 8);

        RGB2Gray.realize(2116, 3538);

        //RGB2Gray.compile_to_coli("cvtcolor.cpp", {in}, "cvtcolor");
    }

    if (0) {
        ImageParam in{Float(32), 2, "input"};
        //Image<float> in(2116, 3538, 2);

        Param<float> a0{"a0"};
        Param<float> a1{"a1"};
        Param<float> a2{"a2"};

        Func rec_filter{"rec_filter"};
        Var x, y;
        RDom r(2, in.width()-3, 0, in.height()-1);

        rec_filter(x, y) = in(x, y);
        rec_filter(r.x, r.y) = a0*rec_filter(r.x, r.y)
                             + a1*rec_filter(r.x-1, r.y)
                             + a2*rec_filter(r.x-2, r.y);


        rec_filter.parallel(y);//.vectorize(x, 8);
        rec_filter.update(0).parallel(r.y);

        //rec_filter.realize(2116, 3538);

        rec_filter.compile_to_coli("recfilter_coli.cpp", {in, a0, a1, a2}, "recfilter_coli");
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

        gaussian.compile_to_coli("gaussian3x3.cpp", {in, kernelX, kernelY}, "gaussian3x3");
    }

    printf("Success!\n");
    return 0;
}
