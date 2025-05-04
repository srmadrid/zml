const std = @import("std");
const types = @import("../types.zig");
const cast = types.cast;

/// Calculate X * Y exactly and store the result in *HI + *LO.  It is
/// given that the values are small enough that no overflow occurs and
/// large enough (or zero) that no underflow occurs.
pub inline fn mul_split64(hi: *f64, lo: *f64, x: f64, y: f64) void {
    if (true) {
        // Fast built-in fused multiply-add.
        hi.* = x * y;
        lo.* = @mulAdd(f64, x, y, -hi.*);
    } else {
        // Apply Dekker's algorithm.
        hi.* = x * y;
        const C: f64 = cast(f64, (1 << (std.math.floatMantissaBits(f64) + 1) / 2) + 1, .{});
        var x1: f64 = x * C;
        var y1: f64 = y * C;
        x1 = (x - x1) + x1;
        y1 = (y - y1) + y1;
        const x2: f64 = x - x1;
        const y2: f64 = y - y1;
        lo.* = (((x1 * y1 - hi.*) + x1 * y2) + x2 * y1) + x2 * y2;
    }
}

/// Add a + b exactly, such that *hi + *lo = a + b.
/// Assumes |a| >= |b| and rounding to nearest.
pub inline fn fast_two_sum64(hi: *f64, lo: *f64, a: f64, b: f64) void {
    hi.* = a + b;
    const e: f64 = hi.* - a; // exact
    lo.* = b - e; // exact
    // Now *hi + *lo = a + b exactly.
}

/// Multiplication of two floating-point expansions: *hi + *lo is an
/// approximation of (h1+l1)*(h2+l2), assuming |l1| <= 1/2*ulp(h1)
/// and |l2| <= 1/2*ulp(h2) and rounding to nearest.
pub inline fn mul_expansion64(hi: *f64, lo: *f64, h1: f64, l1: f64, h2: f64, l2: f64) void {
    mul_split64(hi, lo, h1, h2);
    const r: f64 = h1 * l2 + h2 * l1;
    // Now add r to (hi,lo) using fast two-sum, where we know |r| < |hi|.
    var e: f64 = undefined;
    fast_two_sum64(hi, &e, hi.*, r);
    lo.* -= e;
}

/// Calculate X / Y and store the approximate result in *HI + *LO.  It is
/// assumed that Y is not zero, that no overflow nor underflow occurs, and
/// rounding is to nearest.
pub inline fn div_split64(hi: *f64, lo: *f64, x: f64, y: f64) void {
    hi.* = x / y;
    var a: f64 = undefined;
    var b: f64 = undefined;
    mul_split64(&a, &b, hi.*, y);
    // a + b = hi*y, which should be near x.
    a = x - a; // huge cancellation
    a = a - b;
    // Now x ~ hi*y + a thus x/y ~ hi + a/y.
    lo.* = a / y;
}

/// Division of two floating-point expansions: *hi + *lo is an
/// approximation of (h1+l1)/(h2+l2), assuming |l1| <= 1/2*ulp(h1)
/// and |l2| <= 1/2*ulp(h2), h2+l2 is not zero, and rounding to nearest.  */
pub inline fn div_expansion64(hi: *f64, lo: *f64, h1: f64, l1: f64, h2: f64, l2: f64) void {
    div_split64(hi, lo, h1, h2);
    // (h1+l1)/(h2+l2) ~ h1/h2 + (l1*h2 - l2*h1)/h2^2
    const r: f64 = (l1 * h2 - l2 * h1) / (h2 * h2);
    // Now add r to (hi,lo) using fast two-sum, where we know |r| < |hi|.
    var e: f64 = undefined;
    fast_two_sum64(hi, &e, hi.*, r);
    lo.* += e;
    // Renormalize since |lo| might be larger than 0.5 ulp(hi).
    fast_two_sum64(hi, lo, hi.*, lo.*);
}

/// Calculate X * Y exactly and store the result in *HI + *LO.  It is
/// given that the values are small enough that no overflow occurs and
/// large enough (or zero) that no underflow occurs.
pub inline fn mul_split128(hi: *f128, lo: *f128, x: f128, y: f128) void {
    if (true) {
        // Fast built-in fused multiply-add.
        hi.* = x * y;
        lo.* = @mulAdd(f128, x, y, -hi.*);
    } else {
        // Apply Dekker's algorithm.
        hi.* = x * y;
        const C: f128 = cast(f128, (1 << (std.math.floatMantissaBits(f128) + 1) / 2) + 1, .{});
        var x1: f128 = x * C;
        var y1: f128 = y * C;
        x1 = (x - x1) + x1;
        y1 = (y - y1) + y1;
        const x2: f128 = x - x1;
        const y2: f128 = y - y1;
        lo.* = (((x1 * y1 - hi.*) + x1 * y2) + x2 * y1) + x2 * y2;
    }
}
