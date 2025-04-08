const std = @import("std");

const CN = 134217729.0;

pub inline fn eadd(x: anytype, y: anytype, z: anytype, zz: anytype) void {
    z.* = x + y;
    zz.* = if (@abs(x) > @abs(y)) (x - z.*) + y else (y - z.*) + x;
}

pub inline fn esub(x: anytype, y: anytype, z: anytype, zz: anytype) void {
    z.* = x - y;
    zz.* = if (@abs(x) > @abs(y)) (x - z.*) - y else x - (y + z.*);
}

pub inline fn emulv(x: anytype, y: anytype, z: anytype, zz: anytype) void {
    if (false) {
        z.* = x * y;
        zz.* = @mulAdd(@TypeOf(x, y), x, y, -z.*);
    } else {
        var __p: @TypeOf(x, y) = CN * x;
        const hx: @TypeOf(x, y) = (x - __p) + __p;
        const tx: @TypeOf(x, y) = x - hx;
        __p = CN * y;
        const hy: @TypeOf(x, y) = (y - __p) + __p;
        const ty: @TypeOf(x, y) = y - hy;

        z.* = x * y;
        zz.* = (((hx * hy - z.*) + hx * ty) + tx * hy) + tx * ty;
    }
}
