const std = @import("std");
const types = @import("../types.zig");
const float = @import("../float.zig");

/// Returns a value with the magnitude of `x` and the sign of `y`.
pub fn copysign(x: anytype, y: anytype) @TypeOf(x) {
    comptime if (!types.isFixedPrecision(@TypeOf(x)) or types.isComplex(@TypeOf(x)))
        @compileError("x must be an int or float");

    switch (types.numericType(@TypeOf(x))) {
        .int => {
            switch (types.numericType(@TypeOf(y))) {
                .int => {
                    comptime if (@typeInfo(@TypeOf(x)).int.signedness == .unsigned and @typeInfo(@TypeOf(y)).int.signedness == .signed)
                        @compileError("If x is an unsigned int, y must be also be an unsigned int");

                    switch (@typeInfo(@TypeOf(y)).int.signedness) {
                        .unsigned => return float.abs(x),
                        .signed => {
                            if (y < 0) {
                                return -float.abs(x);
                            } else {
                                return float.abs(x);
                            }
                        },
                    }

                    // To make more efficient
                },
                .float => {
                    comptime if (@typeInfo(@TypeOf(x)).int.signedness == .unsigned)
                        @compileError("If x is an unsigned int, y must be also be an unsigned int");

                    // To implement for int and float
                },
                else => unreachable,
            }
        },
        .float => {
            switch (types.numericType(@TypeOf(y))) {
                .int => {
                    comptime if (@typeInfo(@TypeOf(x)).float.signedness == .unsigned)
                        @compileError("Not implemented for float and int");

                    // To implement for float and int
                },
                .float => {
                    // To implement for different float types
                    return std.math.copysign(x, y);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}
