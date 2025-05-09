const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const mul_split = @import("mul_split.zig");
const Coerce = types.Coerce;
const EnsureFloat = types.EnsureFloat;
const cast = types.cast;

pub inline fn x2y2m1(y: anytype, x: anytype) EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))) {
    comptime if (!types.isFixedPrecision(@TypeOf(y)) or types.isComplex(@TypeOf(y)))
        @compileError("y must be an int or float");

    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(y))) {
        .int => {
            switch (types.numericType(@TypeOf(x))) {
                .int => return x2y2m1(cast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), y, .{}), cast(EnsureFloat(Coerce(@TypeOf(y), @TypeOf(x))), x, .{})),
                .float => return x2y2m1(cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{}), cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{})),
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(x))) {
                .int => return x2y2m1(cast(Coerce(@TypeOf(y), @TypeOf(x)), y, .{}), cast(Coerce(@TypeOf(y), @TypeOf(x)), x, .{})),
                .float => switch (Coerce(@TypeOf(y), @TypeOf(x))) {
                    f16 => return cast(f16, x2y2m1_32(cast(f32, y, .{}), cast(f32, x, .{})), .{}),
                    f32 => {
                        // glibc/sysdeps/ieee754/dbl-64/x2y2m1f.c
                        return x2y2m1_32(cast(f32, y, .{}), cast(f32, x, .{}));
                    },
                    f64 => {
                        // glibc/sysdeps/ieee754/dbl-64/x2y2m1.c
                        return x2y2m1_64(cast(f64, y, .{}), cast(f64, x, .{}));
                    },
                    f80 => return cast(f80, x2y2m1_128(cast(f128, y, .{}), cast(f128, x, .{})), .{}),
                    f128 => {
                        // glibc/sysdeps/ieee754/ldbl-128/x2y2m1l.c
                        return x2y2m1_128(cast(f128, y, .{}), cast(f128, x, .{}));
                    },
                    else => unreachable,
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

fn x2y2m1_32(x: f32, y: f32) f32 {
    const dx: f64 = cast(f64, x, .{});
    const dy: f64 = cast(f64, y, .{});
    return cast(f32, (dx - 1) * (dx + 1) + dy * dy, .{});
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
    const pd: f64 = math.abs(p);
    const qd: f64 = math.abs(q);
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
    const pd: f128 = math.abs(p);
    const qd: f128 = math.abs(q);
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
