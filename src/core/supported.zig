const std = @import("std");

// To be removed in favor of a more flexible approach.

pub const SupportedNumericType = enum {
    /// u8, u16, u32, u64, u128, i8, i16, i32, i64, i128
    BuiltinInt,
    /// f16, f32, f64, 128
    BuiltinFloat,
    /// bool
    BuiltinBool,
    /// std.math.Complex
    Complex,
    /// BigInt
    CustomInt,
    /// Fraction
    CustomReal,
    /// Complex
    CustomComplex,
    /// Expression
    CustomExpression,
    Unsupported,
};

pub inline fn whatSupportedNumericType(comptime T: type) SupportedNumericType {
    @setEvalBranchQuota(10000);
    const info = @typeInfo(T);

    return switch (info) {
        .Int => .BuiltinInt,
        .Float => .BuiltinFloat,
        .Bool => .BuiltinBool,
        else => { // Maybe just use the ones from the standard library, and fix so this works
            if (T == std.math.big.int.Managed) {
                return .CustomInt;
            } else if (T == std.math.big.Rational) {
                return .CustomReal;
            } else if (T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or T == std.math.Complex(f80) or T == std.math.Complex(f128)) {
                return .Complex;
            } else if (false) {
                return .CustomComplex;
            } else if (false) {
                return .CustomExpression;
            } else {
                @compileError("Unsupported numeric type: " ++ @typeName(T));
            }
        },
    };
}

/// Adds two elements of any supported type and stores the result in the output
/// variable. `left` and `right` must be of the same type, and `out` must be a
/// pointer of that same type.
pub inline fn _add(out: anytype, left: anytype, right: anytype) void {
    const T = @TypeOf(out.*, left, right);
    const supported = whatSupportedNumericType(T);
    switch (supported) {
        .BuiltinInt, .BuiltinFloat => {
            out.* = left + right;
        },
        .Complex => {
            out.* = T.add(left, right);
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
            out.add(left, right);
        },
        .BuiltinBool => @compileError("Addition function not defined for booleans"),
        .Unsupported => unreachable,
    }
}

/// Subtracts two elements of any supported type and stores the result in the
/// output variable. `left` and `right` must be of the same type, and `out` must
/// be a pointer of that same type.
pub inline fn _sub(out: anytype, left: anytype, right: anytype) void {
    const T = @TypeOf(out.*, left, right);
    const supported = whatSupportedNumericType(T);
    switch (supported) {
        .BuiltinInt, .BuiltinFloat => {
            out.* = left - right;
        },
        .Complex => {
            out.* = T.sub(left, right);
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
            out.sub(left, right);
        },
        .BuiltinBool => @compileError("Subtraction function not defined for booleans"),
        .Unsupported => unreachable,
    }
}

/// Multiplies two elements of any supported type and stores the result in the
/// output variable. `left` and `right` must be of the same type, and `out` must
/// be a pointer of that same type.
pub inline fn _mul(out: anytype, left: anytype, right: anytype) void {
    const T = @TypeOf(out.*, left, right);
    const supported = whatSupportedNumericType(T);
    switch (supported) {
        .BuiltinInt, .BuiltinFloat => {
            out.* = left * right;
        },
        .Complex => {
            out.* = T.mul(left, right);
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
            out.mul(left, right);
        },
        .BuiltinBool => @compileError("Multiplication function not defined for booleans"),
        .Unsupported => unreachable,
    }
}

/// Divides two elements of any supported type and stores the result in the
/// output variable. `left` and `right` must be of the same type, and `out` must
/// be a pointer of that same type.
pub inline fn _div(out: anytype, left: anytype, right: anytype) void {
    const T = @TypeOf(out.*, left, right);
    const supported = whatSupportedNumericType(T);
    switch (supported) {
        .BuiltinInt, .BuiltinFloat => {
            out.* = left / right;
        },
        .Complex => {
            out.* = T.div(left, right);
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
            out.div(left, right);
        },
        .BuiltinBool => @compileError("Division function not defined for booleans"),
        .Unsupported => unreachable,
    }
}

/// Returns the scalar type of a given numeric type. For example, if `T` is
/// `std.math.Complex(f64)`, this function will return `f64`, and if `T` is
/// `f64`, this function will return `f64`.
pub fn scalar(comptime T: type) type {
    const supported = whatSupportedNumericType(T);
    switch (supported) {
        .BuiltinBool, .BuiltinInt, .BuiltinFloat => return T,
        .Complex => return std.meta.FieldType(T, .re), // change with @FieldType in zig 0.14
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => @compileError("Unsupported numeric type"),
        .Unsupported => unreachable,
    }
}
