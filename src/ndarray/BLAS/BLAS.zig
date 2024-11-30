const std = @import("std");
const ci = @import("../../c.zig");
const Complex = std.math.Complex;

// Level 1 BLAS
pub const asum = @import("asum.zig").asum;
pub fn sasum(n: isize, x: [*]const f32, incx: isize) f32 {
    //return c.cblas_sasum(n, x, incx);
    return asum(f32, n, x, incx);
}
pub fn dasum(n: isize, x: [*]const f64, incx: isize) f64 {
    return asum(f64, n, x, incx);
}
pub fn scasum(n: isize, x: [*]const Complex(f32), incx: isize) f32 {
    return asum(Complex(f32), n, x, incx);
}
pub fn dzasum(n: isize, x: [*]const Complex(f64), incx: isize) f64 {
    return asum(Complex(f64), n, x, incx);
}

pub const axpy = @import("axpy.zig").axpy;
pub fn saxpy(n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
    axpy(f32, n, alpha, x, incx, y, incy);
}
pub fn daxpy(n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
    axpy(f64, n, alpha, x, incx, y, incy);
}
pub fn caxpy(n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize) void {
    axpy(Complex(f32), n, alpha, x, incx, y, incy);
}
pub fn zaxpy(n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize) void {
    axpy(Complex(f64), n, alpha, x, incx, y, incy);
}

pub const copy = @import("copy.zig").copy;
pub fn scopy(n: isize, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
    copy(f32, n, x, incx, y, incy);
}
pub fn dcopy(n: isize, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
    copy(f64, n, x, incx, y, incy);
}
pub fn ccopy(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize) void {
    copy(Complex(f32), n, x, incx, y, incy);
}
pub fn zcopy(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize) void {
    copy(Complex(f64), n, x, incx, y, incy);
}

pub const dot = @import("dot.zig").dot;
pub fn sdot(n: isize, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize) f32 {
    return dot(f32, n, x, incx, y, incy);
}
pub fn ddot(n: isize, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize) f64 {
    return dot(f64, n, x, incx, y, incy);
}

pub const dotc = @import("dotc.zig").dotc;
pub fn cdotc(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize) Complex(f32) {
    return dotc(Complex(f32), n, x, incx, y, incy);
}
pub fn zdotc(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize) Complex(f64) {
    return dotc(Complex(f64), n, x, incx, y, incy);
}

pub const dotc_sub = @import("dotc_sub.zig").dotc_sub;
pub fn cdotc_sub(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, ret: *Complex(f32)) void {
    dotc_sub(Complex(f32), n, x, incx, y, incy, ret);
}
pub fn zdotc_sub(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, ret: *Complex(f64)) void {
    dotc_sub(Complex(f64), n, x, incx, y, incy, ret);
}

pub const dotu = @import("dotu.zig").dotu;
pub fn cdotu(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize) Complex(f32) {
    return dotu(Complex(f32), n, x, incx, y, incy);
}
pub fn zdotu(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize) Complex(f64) {
    return dotu(Complex(f64), n, x, incx, y, incy);
}

pub const dotu_sub = @import("dotu_sub.zig").dotu_sub;
pub fn cdotu_sub(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, ret: *Complex(f32)) void {
    dotu_sub(Complex(f32), n, x, incx, y, incy, ret);
}
pub fn zdotu_sub(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, ret: *Complex(f64)) void {
    dotu_sub(Complex(f64), n, x, incx, y, incy, ret);
}

pub const nrm2 = @import("nrm2.zig").nrm2;
pub fn snrm2(n: isize, x: [*]const f32, incx: isize) f32 {
    return nrm2(f32, n, x, incx);
}
pub fn dnrm2(n: isize, x: [*]const f64, incx: isize) f64 {
    return nrm2(f64, n, x, incx);
}
pub fn scnrm2(n: isize, x: [*]const Complex(f32), incx: isize) f32 {
    return nrm2(Complex(f32), n, x, incx);
}
pub fn dznrm2(n: isize, x: [*]const Complex(f64), incx: isize) f64 {
    return nrm2(Complex(f64), n, x, incx);
}

pub const rot = @import("rot.zig").rot;
pub fn srot(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, c: f32, s: f32) void {
    rot(f32, n, x, incx, y, incy, c, s);
}
pub fn drot(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, c: f64, s: f64) void {
    rot(f64, n, x, incx, y, incy, c, s);
}
pub fn crot(n: isize, x: [*]Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize, c: f32, s: f32) void {
    rot(Complex(f32), n, x, incx, y, incy, c, s);
}
pub fn zrot(n: isize, x: [*]Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize, c: f64, s: f64) void {
    rot(Complex(f64), n, x, incx, y, incy, c, s);
}

test {
    std.testing.refAllDeclsRecursive(@This());
}
