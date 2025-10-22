const std = @import("std");

const types = @import("../types.zig");
const rational = @import("../rational.zig");
const real = @import("../real.zig");
const complex = @import("../complex.zig");
const Complex = complex.Complex;

pub fn neg(allocator: ?std.mem.Allocator, x: anytype) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (types.numericType(X) != .complex)
        @compileError("rational.neg requires x to be a complex, got " ++ @typeName(X));

    var result: X = undefined;
    if (allocator) |a| {
        result = try x.copy(a);
    } else {
        result = x;
    }

    if (comptime types.Scalar(X) == rational.Rational) {
        result.re.num.positive = !result.re.num.positive;
        result.im.num.positive = !result.im.num.positive;
    } else { // real
        result.re.rational.num.positive = !result.re.rational.num.positive;
        result.im.rational.num.positive = !result.im.rational.num.positive;
    }

    result.flags.owns_data = false;
    return result;
}
