const std = @import("std");
const NDArray = @import("../ndarray/ndarray.zig").NDArray;

pub const cf16 = @import("types/cfloat.zig").cf16;
pub const cf32 = @import("types/cfloat.zig").cf32;
pub const cf64 = @import("types/cfloat.zig").cf64;
pub const cf80 = @import("types/cfloat.zig").cf80;
pub const cf128 = @import("types/cfloat.zig").cf128;
pub const comptime_complex = @import("types/cfloat.zig").comptime_complex;
pub const IntegerUnmanaged = @import("types/integer.zig").IntegerUnmanaged;
pub const Integer = @import("types/integer.zig").Integer;
pub const RationalUnmanaged = @import("types/rational.zig").RationalUnmanaged;
pub const Rational = @import("types/rational.zig").Rational;
pub const RealUnmanaged = @import("types/real.zig").RealUnmanaged;
pub const Real = @import("types/real.zig").Real;
pub const ComplexUnmanaged = @import("types/complex.zig").ComplexUnmanaged;
pub const Complex = @import("types/complex.zig").Complex;

//pub const ExpressionUnmanaged = @import("../expression/expression.zig").ExpressionUnmanaged;
//pub const Expression = @import("../expression/expression.zig").Expression;

pub const NumericType = enum {
    /// u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, comptime_int
    int,
    /// f16, f32, f64, f80, f128, comptime_float
    float,
    /// bool
    bool,
    /// cf16, cf32, cf64, cf80, cf128, comptime_complex
    cfloat,
    /// Integer
    integer,
    /// Rational
    rational,
    /// Real
    real,
    /// Complex
    complex,
    /// Expression
    expression,
    unsupported,
};

pub inline fn numericType(comptime T: type) NumericType {
    @setEvalBranchQuota(10000);

    switch (@typeInfo(T)) {
        .int, .comptime_int => return .int,
        .float, .comptime_float => return .float,
        .bool => return .bool,
        else => {
            if (T == cf16 or T == cf32 or T == cf64 or T == cf80 or T == cf128 or T == comptime_complex or T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_complex)) {
                return .cfloat;
            } else if (T == IntegerUnmanaged or T == Integer) {
                return .integer;
            } else if (T == RationalUnmanaged or T == Rational) {
                return .rational;
            } else if (T == RealUnmanaged or T == Real) {
                return .real;
            } else if (T == ComplexUnmanaged or T == Complex) {
                return .complex;
                //} else if (T == ExpressionUnmanaged or T == Expression) {
                //    return .expression;
            } else {
                @compileError("Unsupported numeric type: " ++ @typeName(T));
            }
        },
    }
}

/// Checks if the input type is an instance of an NDArray.
pub fn isNDArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            if (!@hasDecl(T, "empty")) return false;

            const N = NDArray(f64);
            const ninfo = @typeInfo(N).@"struct";
            const t: T = .empty;
            const n: N = .empty;

            inline for (tinfo.fields, ninfo.fields) |field_tinfo, field_ninfo| {
                comptime if (std.mem.eql(u8, field_ninfo.name, "data")) {
                    if (!std.mem.eql(u8, field_tinfo.name, "data")) return false;

                    switch (@typeInfo(@TypeOf(@field(t, "data")))) {
                        .pointer => |p| {
                            if (p.size != .slice) return false;
                        },
                        else => return false,
                    }
                } else if (std.mem.eql(u8, field_ninfo.name, "base")) {
                    if (!std.mem.eql(u8, field_tinfo.name, "base")) return false;

                    switch (@typeInfo(@TypeOf(@field(t, "base")))) {
                        .optional => continue, // Checking if child type is an NDArray would lead to infinite recursion
                        else => return false,
                    }
                } else {
                    if (!std.meta.eql(@field(t, field_tinfo.name), @field(n, field_ninfo.name))) return false;
                };
            }

            inline for (tinfo.decls, ninfo.decls) |decl_tinfo, decl_ninfo| {
                if (!std.meta.eql(decl_tinfo, decl_ninfo)) return false;
            }

            if (tinfo.is_tuple != ninfo.is_tuple) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input is an array.
pub fn isArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .slice) return false;

            return true;
        },
        else => return false,
    }
}

/// Coerces the types of the input variables.
pub fn Coerce(comptime K: type, comptime V: type) type {
    std.debug.print("UNFINISHED: Coerce\n", .{});
    comptime var T1: type = undefined;
    comptime var T2: type = undefined;

    if (isNDArray(K)) {
        const t1: K = .empty;
        T1 = @typeInfo(@TypeOf(@field(t1, "data"))).pointer.child;
    } else if (isArray(K)) {
        T1 = @typeInfo(K).pointer.child;
    }

    const t1numeric = numericType(T1);

    if (isNDArray(V)) {
        const t2: V = .empty;
        T2 = @typeInfo(@TypeOf(@field(t2, "data"))).pointer.child;
    } else if (isArray(V)) {
        T2 = @typeInfo(V).pointer.child;
    }

    const t2numeric = numericType(T2);

    switch (t1numeric) {
        .int => switch (t2numeric) {
            .int => return @TypeOf(K, V),
        },
    }
}

/// Returns the scalar type of a given numeric type. For example, if `T` is
/// `std.math.Complex(f64)`, this function will return `f64`, and if `T` is
/// `f64`, this function will return `f64`.
///
/// Edit so it returns the scalar of an array (NDarray(T) -> T, []T -> T, so if
/// NDarray(Complex(f64) -> Complex(f64)) but keep Complex(T) -> T.
pub fn Numeric(comptime T: type) type {
    if (isNDArray(T)) {
        @compileLog("NDArray detected");
        const t: T = .empty;
        return @typeInfo(@TypeOf(@field(t, "data"))).pointer.child;
    } else if (isArray(T)) {
        const K = @typeInfo(T).pointer.child;
        numericType(K); // Will raise compile error if K is not a supported numeric type
        return K;
    }

    const numeric = numericType(T);

    switch (numeric) {
        .int => return T,
        .float => return T,
        .bool => return T,
        .cfloat => switch (T) {
            cf16 => return f16,
            cf32 => return f32,
            cf64 => return f64,
            cf80 => return f80,
            cf128 => return f128,
            comptime_complex => return comptime_complex,
            std.math.Complex(f16) => return f16,
            std.math.Complex(f32) => return f32,
            std.math.Complex(f64) => return f64,
            std.math.Complex(f80) => return f80,
            std.math.Complex(f128) => return f128,
            std.math.Complex(comptime_complex) => return comptime_complex,
            else => @compileError("???"),
        },
        .integer => return Integer,
        .rational => return Rational,
        .real => return Real,
        .complex => return Rational,
        //.expression => return Expression,
        else => unreachable,
    }
}
