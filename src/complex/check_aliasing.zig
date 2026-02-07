const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

const rational_check_aliasing = @import("../rational/check_aliasing.zig").check_aliasing;
// const real_check_aliasing = @import("../real/check_aliasing.zig").check_aliasing;

pub fn check_aliasing(o: anytype, x: anytype) bool {
    switch (comptime types.numericType(types.Scalar(types.Child(@TypeOf(o))))) {
        .rational => switch (comptime types.numericType(@TypeOf(x))) {
            .bool, .int, .float, .cfloat => return false,
            .integer, .rational => return rational_check_aliasing(&o.re, x) or rational_check_aliasing(&o.im, x),
            .real => @compileError("zml.complex.check_aliasing not implemented yet for complex and real types"),
            .complex => switch (comptime types.numericType(types.Scalar(types.Child(@TypeOf(x))))) {
                .rational => return rational_check_aliasing(&o.re, x.re) or rational_check_aliasing(&o.re, x.im) or
                    rational_check_aliasing(&o.im, x.re) or rational_check_aliasing(&o.im, x.im),
                .real => @compileError("zml.complex.check_aliasing not implemented yet for complex and real types"),
                else => unreachable,
            },
            else => unreachable,
        },
        .real => @compileError("zml.complex.check_aliasing not implemented yet for real output types"),
        else => unreachable,
    }
}
