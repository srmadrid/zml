const std = @import("std");
const NDArray = @import("../ndarray/ndarray.zig").NDArray;

pub const cfloat = @import("types/cfloat.zig");
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
    /// bool
    bool,
    /// u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, comptime_int
    int,
    /// f16, f32, f64, f80, f128, comptime_float
    float,
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
        .bool => return .bool,
        .int, .comptime_int => {
            if (T != u8 and T != u16 and T != u32 and T != u64 and T != u128 and T != i8 and T != i16 and T != i32 and T != i64 and T != i128 and T != comptime_int) {
                @compileError("Unsupported integer type: " ++ @typeName(T));
            } else {
                return .int;
            }
        },
        .float, .comptime_float => return .float,
        else => {
            if (T == cf16 or T == cf32 or T == cf64 or T == cf80 or T == cf128 or T == comptime_complex or T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_complex)) {
                return .cfloat;
            } else if (T == IntegerUnmanaged or T == Integer) {
                return .integer;
            } else if (T == RationalUnmanaged or T == Rational) {
                return .rational;
            } else if (T == RealUnmanaged or T == Real) {
                return .real;
            } else if (T == ComplexUnmanaged(IntegerUnmanaged) or T == Complex(IntegerUnmanaged) or T == ComplexUnmanaged(RationalUnmanaged) or T == Complex(RationalUnmanaged) or T == ComplexUnmanaged(RealUnmanaged) or T == Complex(RealUnmanaged)) {
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

            if (tinfo.fields.len != ninfo.fields.len) return false;
            if (tinfo.decls.len != ninfo.decls.len) return false;

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

/// Checks if the input is a slice.
pub fn isSlice(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .slice) return false;

            return true;
        },
        else => return false,
    }
}

/// Coerces the input types to the smallest numeric type that can represent both
/// types. If either type is a slice or an NDArray, the corresponding type of
/// the elements is used for the coercion.
pub fn Coerce(comptime K: type, comptime V: type) type {
    //@compileLog("UNFINISHED: Coerce\n");
    comptime var T1: type = K;
    comptime var T2: type = V;

    if (isNDArray(K)) {
        const t1: K = .empty;
        T1 = @typeInfo(@TypeOf(@field(t1, "data"))).pointer.child;
    } else if (isSlice(K)) {
        T1 = @typeInfo(K).pointer.child;
    }

    const t1numeric = numericType(T1);

    if (isNDArray(V)) {
        const t2: V = .empty;
        T2 = @typeInfo(@TypeOf(@field(t2, "data"))).pointer.child;
    } else if (isSlice(V)) {
        T2 = @typeInfo(V).pointer.child;
    }

    const t2numeric = numericType(T2);

    if (T1 == T2) {
        return T1;
    }

    switch (t1numeric) {
        .bool => switch (t2numeric) {
            .bool => return bool,
            else => return T2,
        },
        .int => switch (t2numeric) {
            .bool => return T1,
            .int => {
                if (T1 == comptime_int or T2 == comptime_int) {
                    return comptime_int;
                } else {
                    comptime var t1info = @typeInfo(T1);
                    comptime var t2info = @typeInfo(T2);

                    if (t1info.int.signedness == .unsigned) {
                        if (t2info.int.signedness == .unsigned) {
                            if (t1info.int.bits > t2info.int.bits) {
                                return T1;
                            } else {
                                return T2;
                            }
                        } else {
                            if (t1info.int.bits > t2info.int.bits) {
                                if (t1info.int.bits == 128) {
                                    @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2) ++ " as " ++ @typeName(T1) ++ " already has the maximum amount of bits.");
                                } else {
                                    t2info.int.bits = t1info.int.bits * 2;
                                    return @Type(t2info);
                                }
                            } else if (t1info.int.bits == t2info.int.bits) {
                                t2info.int.bits *= 2;
                                return @Type(t2info);
                            } else {
                                return T2;
                            }
                        }
                    } else {
                        if (t2info.int.signedness == .unsigned) {
                            if (t2info.int.bits > t1info.int.bits) {
                                if (t2info.int.bits == 128) {
                                    @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2) ++ " as " ++ @typeName(T2) ++ " already has the maximum amount of bits.");
                                } else {
                                    t1info.int.bits = t2info.int.bits * 2;
                                    return @Type(t1info);
                                }
                            } else if (t1info.int.bits == t2info.int.bits) {
                                t1info.int.bits *= 2;
                                return @Type(t1info);
                            } else {
                                return T1;
                            }
                        } else {
                            if (t1info.int.bits > t2info.int.bits) {
                                return T1;
                            } else {
                                return T2;
                            }
                        }
                    }
                }
            },
            .unsupported => unreachable,
            else => return T2,
        },
        .float => {
            switch (t2numeric) {
                .bool => return T1,
                .int => return T1,
                .float => {
                    const t1info = @typeInfo(T1);
                    const t2info = @typeInfo(T2);

                    if (t1info.float.bits > t2info.float.bits) {
                        return T1;
                    } else {
                        return T2;
                    }
                },
                .cfloat => {
                    const t1info = @typeInfo(T1);
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    const t2info = @typeInfo(@TypeOf(@field(t2, "re")));

                    if (t1info.float.bits > t2info.float.bits) {
                        return cfloat.cfloat(T1);
                    } else {
                        return T2;
                    }
                },
                .unsupported => unreachable,
                else => return T2,
            }
        },
        .cfloat => {
            switch (t2numeric) {
                .bool => return T1,
                .int => return T1,
                .float => return T1,
                .cfloat => {
                    const t1: T1 = .{ .re = 0, .im = 0 };
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    return cfloat.cfloat(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                },
                .integer => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
                .rational => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
                .real => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
                .unsupported => unreachable,
                else => return T2,
            }
        },
        .integer => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
            .integer => return T1,
            .unsupported => unreachable,
            else => return T2,
        },
        .rational => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
            .integer => return T1,
            .rational => return T1,
            .complex => {
                if (T2 == ComplexUnmanaged(IntegerUnmanaged)) {
                    return ComplexUnmanaged(RationalUnmanaged);
                } else if (T2 == Complex(IntegerUnmanaged)) {
                    return Complex(RationalUnmanaged);
                } else {
                    return T2;
                }
            },
            .unsupported => unreachable,
            else => return T2,
        },
        .real => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
            .integer => return T1,
            .rational => return T1,
            .real => return T1,
            .complex => {
                if (T2 == ComplexUnmanaged(IntegerUnmanaged)) {
                    return ComplexUnmanaged(RealUnmanaged);
                } else if (T2 == Complex(IntegerUnmanaged)) {
                    return Complex(RealUnmanaged);
                } else if (T2 == ComplexUnmanaged(RationalUnmanaged)) {
                    return ComplexUnmanaged(RealUnmanaged);
                } else if (T2 == Complex(RationalUnmanaged)) {
                    return Complex(RealUnmanaged);
                } else {
                    return T2;
                }
            },
            .unsupported => unreachable,
            else => return T2,
        },
        .complex => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
            .integer => return T1,
            .rational => return T1,
            .real => return T1,
            .complex => {
                const t1: T1 = .empty;
                const t2: T2 = .empty;
                if (@hasField(T1, "allocator") or @hasField(T2, "allocator")) {
                    return Complex(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                } else {
                    return ComplexUnmanaged(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                }
            },
            .expression => return T2,
            else => unreachable,
        },
        .expression => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => @compileError("Cannot coerce " ++ @typeName(T1) ++ " and " ++ @typeName(T2)),
            .integer => return T1,
            .rational => return T1,
            .real => return T1,
            .complex => return T1,
            .expression => return T1,
            else => unreachable,
        },
        .unsupported => unreachable,
    }
}

/// Returns the scalar type of a given numeric type:
/// - Atomic types, such as `u8`, `i16`, `f32`, etc., are returned as-is.
/// - Complex types, such as `cf32`, `cf64`, etc., are returned as their real part.
/// - Real arbitrary precision types, such as `Integer`, `Rational`, etc., are returned as-is.
/// - Complex arbitrary precision types are returned as their real part.
/// - Other types are returned as `unsupported`.
pub fn Numeric(comptime T: type) type {
    if (isNDArray(T)) {
        const t: T = .empty;
        return @typeInfo(@TypeOf(@field(t, "data"))).pointer.child;
    } else if (isSlice(T)) {
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
            else => unreachable,
        },
        .integer => return Integer,
        .rational => return Rational,
        .real => return Real,
        .complex => switch (T) {
            ComplexUnmanaged(IntegerUnmanaged) => return IntegerUnmanaged,
            ComplexUnmanaged(RationalUnmanaged) => return RationalUnmanaged,
            ComplexUnmanaged(RealUnmanaged) => return RealUnmanaged,
            Complex(IntegerUnmanaged) => return IntegerUnmanaged,
            Complex(RationalUnmanaged) => return RationalUnmanaged,
            Complex(RealUnmanaged) => return RealUnmanaged,
            else => unreachable,
        },
        //.expression => return Expression,
        else => unreachable,
    }
}

/// Casts a value of any numeric type to any other numeric type. The optional
/// allocator is needed only if the type to be casted to requires allocation.
pub inline fn cast(comptime T: type, val: anytype) T {
    std.debug.print("UNFINISHED: cast\n", .{});
    _ = val;
    return 0;
}
