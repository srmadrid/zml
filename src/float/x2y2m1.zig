const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");
const mul_split = @import("mul_split.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const scast = types.scast;

pub inline fn x2y2m1(x: anytype, y: anytype) EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y))) {
    comptime if (types.numericType(@TypeOf(x)) != .int and types.numericType(@TypeOf(x)) != .float)
        @compileError("float.x2y2m1: x must be an int or float, got " ++ @typeName(@TypeOf(x)));

    comptime if (types.numericType(@TypeOf(y)) != .int and types.numericType(@TypeOf(y)) != .float)
        @compileError("float.x2y2m1: y must be an int or float, got " ++ @typeName(@TypeOf(y)));

    switch (EnsureFloat(Coerce(@TypeOf(x), @TypeOf(y)))) {
        f16 => return scast(f16, x2y2m1_32(scast(f32, x), scast(f32, y))),
        f32 => {
            // glibc/sysdeps/ieee754/dbl-64/x2y2m1f.c
            return x2y2m1_32(scast(f32, x), scast(f32, y));
        },
        f64 => {
            // glibc/sysdeps/ieee754/dbl-64/x2y2m1.c
            return x2y2m1_64(scast(f64, x), scast(f64, y));
        },
        f80 => return scast(f80, x2y2m1_128(scast(f128, x), scast(f128, y))),
        f128 => {
            // glibc/sysdeps/ieee754/ldbl-128/x2y2m1l.c
            return x2y2m1_128(scast(f128, x), scast(f128, y));
        },
        else => unreachable,
    }
}

fn x2y2m1_32(x: f32, y: f32) f32 {
    const dx: f64 = scast(f64, x);
    const dy: f64 = scast(f64, y);
    return scast(f32, (dx - 1) * (dx + 1) + dy * dy);
}

// Calculate X + Y exactly and store the result in *HI + *LO.  It is
// given that |X| >= |Y| and the values are small enough that no
// overflow occurs.  */
inline fn add_split64(hi: *f64, lo: *f64, x: f64, y: f64) void {
    // Apply Dekker's algorithm.
    hi.* = x + y;
    lo.* = (x - hi.*) + y;
}

// Compare absolute values of floating-point values pointed to by P
// and Q for qsort.
fn compare64(p: f64, q: f64) bool {
    const pd: f64 = float.abs(p);
    const qd: f64 = float.abs(q);
    if (pd < qd) {
        return true;
    } else {
        return false;
    }
}

fn x2y2m1_64(x: f64, y: f64) f64 {
    var vals: [5]f64 = undefined;
    mul_split.mul_split64(&vals[1], &vals[0], x, x);
    mul_split.mul_split64(&vals[3], &vals[2], y, y);
    vals[4] = -1;
    sortSmall(f64, &vals, compare64);
    // Add up the values so that each element of VALS has absolute value
    // at most equal to the last set bit of the next nonzero
    // element.
    var i: u32 = 0;
    while (i <= 3) {
        add_split64(&vals[i + 1], &vals[i], vals[i + 1], vals[i]);
        sortSmall(f64, vals[(i + 1)..], compare64);

        i += 1;
    }
    // Now any error from this addition will be small.
    return vals[4] + vals[3] + vals[2] + vals[1] + vals[0];
}

// Calculate X + Y exactly and store the result in *HI + *LO.  It is
// given that |X| >= |Y| and the values are small enough that no
// overflow occurs.  */
inline fn add_split128(hi: *f128, lo: *f128, x: f128, y: f128) void {
    // Apply Dekker's algorithm.
    hi.* = x + y;
    lo.* = (x - hi.*) + y;
}

// Compare absolute values of floating-point values pointed to by P
// and Q for qsort.
fn compare128(p: f128, q: f128) bool {
    const pd: f128 = float.abs(p);
    const qd: f128 = float.abs(q);
    if (pd < qd) {
        return true;
    } else {
        return false;
    }
}

fn x2y2m1_128(x: f128, y: f128) f128 {
    var vals: [5]f128 = undefined;
    mul_split.mul_split128(&vals[1], &vals[0], x, x);
    mul_split.mul_split128(&vals[3], &vals[2], y, y);
    vals[4] = -1;
    sortSmall(f128, &vals, compare128);
    // Add up the values so that each element of VALS has absolute value
    // at most equal to the last set bit of the next nonzero
    // element.
    var i: u32 = 0;
    while (i <= 3) {
        add_split128(&vals[i + 1], &vals[i], vals[i + 1], vals[i]);
        sortSmall(f128, vals[(i + 1)..], compare128);

        i += 1;
    }
    // Now any error from this addition will be small.
    return vals[4] + vals[3] + vals[2] + vals[1] + vals[0];
}

fn sortSmall(comptime T: type, arr: []T, isLessThan: fn (T, T) bool) void {
    var i: i32 = 1;
    while (i < arr.len) {
        const key: T = arr[@intCast(i)];
        var j: i32 = i - 1;
        while (j >= 0 and isLessThan(key, arr[@intCast(j)])) {
            arr[@intCast(j + 1)] = arr[@intCast(j)];

            j -= 1;
        }

        arr[@intCast(j + 1)] = key;

        i += 1;
    }
}
