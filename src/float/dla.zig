const std = @import("std");

const CN = 134217729.0;

pub inline fn eadd(x: f64, y: f64, z: *f64, zz: *f64) void {
    z.* = x + y;
    zz.* = if (@abs(x) > @abs(y)) (x - z.*) + y else (y - z.*) + x;
}

pub inline fn esub(x: f64, y: f64, z: *f64, zz: *f64) void {
    z.* = x - y;
    zz.* = if (@abs(x) > @abs(y)) (x - z.*) - y else x - (y + z.*);
}

pub inline fn emulv(x: f64, y: f64, z: *f64, zz: *f64) void {
    if (false) {
        z.* = x * y;
        zz.* = @mulAdd(f64, x, y, -z.*);
    } else {
        var __p: f64 = CN * x;
        const hx: f64 = (x - __p) + __p;
        const tx: f64 = x - hx;
        __p = CN * y;
        const hy: f64 = (y - __p) + __p;
        const ty: f64 = y - hy;

        z.* = x * y;
        zz.* = (((hx * hy - z.*) + hx * ty) + tx * hy) + tx * ty;
    }
}

pub inline fn mul12(x: f64, y: f64, z: *f64, zz: *f64) void {
    if (false) {
        emulv(x, y, z, zz);
    } else {
        var __p: f64 = CN * x;
        const hx: f64 = (x - __p) + __p;
        const tx: f64 = x - hx;
        __p = CN * y;
        const hy: f64 = (y - __p) + __p;
        const ty: f64 = y - hy;
        __p = hx * hy;
        const __q: f64 = hx * ty + tx * hy;
        z.* = __p + __q;
        zz.* = ((__p - z.*) + __q) + tx * ty;
    }
}

pub inline fn mul2(x: f64, xx: f64, y: f64, yy: f64, z: *f64, zz: *f64, c: *f64, cc: *f64) void {
    mul12(x, y, c, cc);
    cc.* = (x * yy + xx * y) + cc.*;
    z.* = c.* + cc.*;
    zz.* = (c.* - z.*) + cc.*;
}

pub inline fn div2(x: f64, xx: f64, y: f64, yy: f64, z: *f64, zz: *f64, c: *f64, cc: *f64, u: *f64, uu: *f64) void {
    c.* = x / y;
    mul12(c.*, y, u, uu);
    cc.* = ((((x - u.*) - uu.*) + xx) - c.* * yy) / y;
    z.* = c.* + cc.*;
    zz.* = (c.* - z.*) + cc.*;
}
