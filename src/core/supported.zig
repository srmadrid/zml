const std = @import("std");

// To be removed in favor of a more flexible approach.

pub const SupportedNumericType = enum {
    /// u8, u16, u32, u64, u128, i8, i16, i32, i64, i128
    BuiltinInt,
    /// f16, f32, f64, 128
    BuiltinFloat,
    /// bool
    BuiltinBool,
    /// BigInt
    CustomInt,
    /// Fraction
    CustomReal,
    /// cf32, cf64, cf128
    CustomComplexFloat,
    /// Complex
    CustomComplex,
    /// Expression
    CustomExpression,
    Unsupported,
};

pub inline fn whatSupportedNumericType(comptime T: type) SupportedNumericType {
    const info = @typeInfo(T);
    const name = @typeName(T);

    return switch (info) {
        .Int => SupportedNumericType.BuiltinInt,
        .Float => SupportedNumericType.BuiltinFloat,
        .Bool => SupportedNumericType.BuiltinBool,
        else => { // Maybe just use the ones from the standard library, and fix so this works
            if (std.mem.eql(u8, name, "bigint.BigInt")) {
                return SupportedNumericType.CustomInt;
            } else if (std.mem.eql(u8, name, "fraction.Fraction")) {
                return SupportedNumericType.CustomReal;
            } else if (std.mem.eql(u8, name, "types.cf.cf(f16)") or
                std.mem.eql(u8, name, "types.cf.cf(f32)") or
                std.mem.eql(u8, name, "types.cf.cf(f64)") or
                std.mem.eql(u8, name, "types.cf.cf(f128)"))
            {
                return SupportedNumericType.CustomComplexFloat;
            } else if (std.mem.eql(u8, name, "complex.Complex")) {
                return SupportedNumericType.CustomComplex;
            } else if (std.mem.eql(u8, name, "expression.Expression")) {
                return SupportedNumericType.CustomExpression;
            } else {
                @compileError("Unsupported numeric type: " ++ name);
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
        .CustomComplexFloat => {
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
        .CustomComplexFloat => {
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
        .CustomComplexFloat => {
            out.* = T.mult(left, right);
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
            out.mult(left, right);
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
        .CustomComplexFloat => {
            out.* = T.div(left, right);
        },
        .CustomInt, .CustomReal, .CustomComplex, .CustomExpression => {
            out.div(left, right);
        },
        .BuiltinBool => @compileError("Division function not defined for booleans"),
        .Unsupported => unreachable,
    }
}
