const std = @import("std");
const NDArray = @import("../ndarray/ndarray.zig").NDArray;

pub const default_uint = usize;
pub const default_int = isize;
pub const default_float = f64;

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
            if (T != u8 and T != u16 and T != u32 and T != u64 and T != u128 and T != usize and T != i8 and T != i16 and T != i32 and T != i64 and T != i128 and T != isize and T != comptime_int) {
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

/// Checks if the input type is of fixed precision.
pub fn isFixedPrecision(comptime T: type) bool {
    switch (numericType(T)) {
        .int => return true,
        .float => return true,
        .cfloat => return true,
        else => return false,
    }
}

/// Checks if the input type is of arbitrary precision.
pub fn isArbitraryPrecision(comptime T: type) bool {
    switch (numericType(T)) {
        .integer => return true,
        .rational => return true,
        .real => return true,
        .complex => return true,
        .expression => return true,
        else => return false,
    }
}

/// Checks if the input type is complex without distunguishing between
/// arbitrary and fixed precision.
pub fn isComplex(comptime T: type) bool {
    switch (numericType(T)) {
        .cfloat => return true,
        .complex => return true,
        else => return false,
    }
}

/// Coerces the input types to the smallest numeric type that can represent both
/// types. If either type is a slice or an NDArray, the corresponding type of
/// the elements is used for the coercion.
pub fn Coerce(comptime K: type, comptime V: type) type {
    comptime var T1: type = K;
    comptime var T2: type = V;

    comptime if (isNDArray(K)) {
        const t1: K = .empty;
        T1 = @typeInfo(@TypeOf(@field(t1, "data"))).pointer.child;
    } else if (isSlice(K)) {
        T1 = @typeInfo(K).pointer.child;
    };

    const t1numeric = numericType(T1);

    comptime if (isNDArray(V)) {
        const t2: V = .empty;
        T2 = @typeInfo(@TypeOf(@field(t2, "data"))).pointer.child;
    } else if (isSlice(V)) {
        T2 = @typeInfo(V).pointer.child;
    };

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
                        return cfloat.Cfloat(T1);
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
                    return cfloat.Cfloat(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                },
                .integer => return ComplexUnmanaged(RationalUnmanaged),
                .rational => return ComplexUnmanaged(RationalUnmanaged),
                .real => return ComplexUnmanaged(RealUnmanaged),
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
            }
        },
        .integer => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => return ComplexUnmanaged(RationalUnmanaged),
            .integer => return T1,
            .unsupported => unreachable,
            else => return T2,
        },
        .rational => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => return ComplexUnmanaged(RationalUnmanaged),
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
            .cfloat => return ComplexUnmanaged(RealUnmanaged),
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
            .cfloat => {
                if (T1 == ComplexUnmanaged(IntegerUnmanaged)) {
                    return ComplexUnmanaged(RationalUnmanaged);
                } else if (T1 == Complex(IntegerUnmanaged)) {
                    return Complex(RationalUnmanaged);
                } else {
                    return T1;
                }
            },
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
            .cfloat => return T1,
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

/// Checks if if `K` can be coerced to `V` without loss of information. This is
/// a more flexible version of `Coerce`, as it does not require `V` to be to the
/// smallest type that can represent both types. The only requirement is that
/// `V` can represent all values of the first two types.
pub fn canCoerce(comptime K: type, comptime V: type) bool {
    comptime var T1: type = K;
    comptime var T2: type = V;

    comptime if (isNDArray(K)) {
        const t1: K = .empty;
        T1 = @typeInfo(@TypeOf(@field(t1, "data"))).pointer.child;
    } else if (isSlice(K)) {
        T1 = @typeInfo(K).pointer.child;
    };

    const t1numeric = numericType(T1);

    comptime if (isNDArray(V)) {
        const t2: V = .empty;
        T2 = @typeInfo(@TypeOf(@field(t2, "data"))).pointer.child;
    } else if (isSlice(V)) {
        T2 = @typeInfo(V).pointer.child;
    };

    const t2numeric = numericType(T2);

    comptime var T3: type = undefined;

    if (T1 == T2) {
        return true;
    }

    switch (t1numeric) {
        .bool => switch (t2numeric) {
            .bool => T3 = bool,
            else => T3 = T2,
        },
        .int => switch (t2numeric) {
            .bool => T3 = T1,
            .int => {
                if (T1 == comptime_int or T2 == comptime_int) {
                    T3 = comptime_int;
                } else {
                    comptime var t1info = @typeInfo(T1);
                    comptime var t2info = @typeInfo(T2);

                    if (t1info.int.signedness == .unsigned) {
                        if (t2info.int.signedness == .unsigned) {
                            if (t1info.int.bits > t2info.int.bits) {
                                T3 = T1;
                            } else {
                                T3 = T2;
                            }
                        } else {
                            if (t1info.int.bits > t2info.int.bits) {
                                if (t1info.int.bits == 128) {
                                    return false;
                                } else {
                                    t2info.int.bits = t1info.int.bits * 2;
                                    T3 = @Type(t2info);
                                }
                            } else if (t1info.int.bits == t2info.int.bits) {
                                t2info.int.bits *= 2;
                                T3 = @Type(t2info);
                            } else {
                                T3 = T2;
                            }
                        }
                    } else {
                        if (t2info.int.signedness == .unsigned) {
                            if (t2info.int.bits > t1info.int.bits) {
                                if (t2info.int.bits == 128) {
                                    return false;
                                } else {
                                    t1info.int.bits = t2info.int.bits * 2;
                                    T3 = @Type(t1info);
                                }
                            } else if (t1info.int.bits == t2info.int.bits) {
                                t1info.int.bits *= 2;
                                T3 = @Type(t1info);
                            } else {
                                T3 = T1;
                            }
                        } else {
                            if (t1info.int.bits > t2info.int.bits) {
                                T3 = T1;
                            } else {
                                T3 = T2;
                            }
                        }
                    }
                }
            },
            .unsupported => unreachable,
            else => T3 = T2,
        },
        .float => {
            switch (t2numeric) {
                .bool => T3 = T1,
                .int => T3 = T1,
                .float => {
                    const t1info = @typeInfo(T1);
                    const t2info = @typeInfo(T2);

                    if (t1info.float.bits > t2info.float.bits) {
                        T3 = T1;
                    } else {
                        T3 = T2;
                    }
                },
                .cfloat => {
                    const t1info = @typeInfo(T1);
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    const t2info = @typeInfo(@TypeOf(@field(t2, "re")));

                    if (t1info.float.bits > t2info.float.bits) {
                        T3 = cfloat.Cfloat(T1);
                    } else {
                        T3 = T2;
                    }
                },
                .unsupported => unreachable,
                else => T3 = T2,
            }
        },
        .cfloat => {
            switch (t2numeric) {
                .bool => T3 = T1,
                .int => T3 = T1,
                .float => T3 = T1,
                .cfloat => {
                    const t1: T1 = .{ .re = 0, .im = 0 };
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    T3 = cfloat.Cfloat(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                },
                .integer => T3 = ComplexUnmanaged(RationalUnmanaged),
                .rational => T3 = ComplexUnmanaged(RationalUnmanaged),
                .real => T3 = ComplexUnmanaged(RealUnmanaged),
                .complex => {
                    if (T2 == ComplexUnmanaged(IntegerUnmanaged)) {
                        T3 = ComplexUnmanaged(RationalUnmanaged);
                    } else if (T2 == Complex(IntegerUnmanaged)) {
                        T3 = Complex(RationalUnmanaged);
                    } else {
                        T3 = T2;
                    }
                },
                .unsupported => unreachable,
                else => T3 = T2,
            }
        },
        .integer => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = ComplexUnmanaged(RationalUnmanaged),
            .integer => T3 = T1,
            .unsupported => unreachable,
            else => T3 = T2,
        },
        .rational => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = ComplexUnmanaged(RationalUnmanaged),
            .integer => T3 = T1,
            .rational => T3 = T1,
            .complex => {
                if (T2 == ComplexUnmanaged(IntegerUnmanaged)) {
                    T3 = ComplexUnmanaged(RationalUnmanaged);
                } else if (T2 == Complex(IntegerUnmanaged)) {
                    T3 = Complex(RationalUnmanaged);
                } else {
                    T3 = T2;
                }
            },
            .unsupported => unreachable,
            else => T3 = T2,
        },
        .real => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = ComplexUnmanaged(RealUnmanaged),
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => {
                if (T2 == ComplexUnmanaged(IntegerUnmanaged)) {
                    T3 = ComplexUnmanaged(RealUnmanaged);
                } else if (T2 == Complex(IntegerUnmanaged)) {
                    T3 = Complex(RealUnmanaged);
                } else if (T2 == ComplexUnmanaged(RationalUnmanaged)) {
                    T3 = ComplexUnmanaged(RealUnmanaged);
                } else if (T2 == Complex(RationalUnmanaged)) {
                    T3 = Complex(RealUnmanaged);
                } else {
                    T3 = T2;
                }
            },
            .unsupported => unreachable,
            else => T3 = T2,
        },
        .complex => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => {
                if (T1 == ComplexUnmanaged(IntegerUnmanaged)) {
                    T3 = ComplexUnmanaged(RationalUnmanaged);
                } else if (T1 == Complex(IntegerUnmanaged)) {
                    T3 = Complex(RationalUnmanaged);
                } else {
                    T3 = T1;
                }
            },
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => {
                const t1: T1 = .empty;
                const t2: T2 = .empty;
                if (@hasField(T1, "allocator") or @hasField(T2, "allocator")) {
                    T3 = Complex(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                } else {
                    T3 = ComplexUnmanaged(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                }
            },
            .expression => T3 = T2,
            else => unreachable,
        },
        .expression => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = T1,
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => T3 = T1,
            .expression => T3 = T1,
            else => unreachable,
        },
        .unsupported => unreachable,
    }

    return T2 == T3;
}

/// Returns the scalar type of a given numeric type, slice, or NDArray:
/// - Atomic types, such as `u8`, `i16`, `f32`, etc., are returned as-is.
/// - Complex types, such as `cf32`, `cf64`, etc., are returned as their real part.
/// - Real arbitrary precision types, such as `Integer`, `Rational`, etc., are returned as-is.
/// - Complex arbitrary precision types are returned as their real part.
/// - Slices and NDArray types are returned as their element type.
pub fn Scalar(comptime T: type) type {
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

/// Casts a value of any numeric type to any other numeric type. It does not
/// check if the cast is valid, so it is up to the user to ensure that the cast
/// is valid. The optional allocator is needed only if the type to be casted to
/// requires allocation to be initialized.
pub inline fn cast(
    comptime T: type,
    value: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) T {
    const T1 = T;
    const T2 = @TypeOf(value);

    _ = options; // To be used to initialize arbitrary precision types

    if (T1 == T2) {
        return value;
    }

    const t1numeric = numericType(T1);
    const t2numeric = numericType(T2);

    switch (t2numeric) {
        .bool => switch (t1numeric) {
            .bool => return value,
            .int => return @intFromBool(value),
            .float => return @floatFromInt(@intFromBool(value)),
            .cfloat => return .{
                .re = @floatFromInt(@intFromBool(value)),
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from bool to integer"),
            .rational => @compileError("Not implemented yet: casting from bool to rational"),
            .real => @compileError("Not implemented yet: casting from bool to real"),
            .complex => @compileError("Not implemented yet: casting from bool to complex"),
            .expression => @compileError("Not implemented yet: casting from bool to expression"),
            .unsupported => unreachable,
        },
        .int => switch (t1numeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from int to integer"),
            .rational => @compileError("Not implemented yet: casting from int to rational"),
            .real => @compileError("Not implemented yet: casting from int to real"),
            .complex => @compileError("Not implemented yet: casting from int to complex"),
            .expression => @compileError("Not implemented yet: casting from int to expression"),
            .unsupported => unreachable,
        },
        .float => switch (t1numeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .cfloat => return .{
                .re = @floatCast(value),
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from float to integer"),
            .rational => @compileError("Not implemented yet: casting from float to rational"),
            .real => @compileError("Not implemented yet: casting from float to real"),
            .complex => @compileError("Not implemented yet: casting from float to complex"),
            .expression => @compileError("Not implemented yet: casting from float to expression"),
            .unsupported => unreachable,
        },
        .cfloat => switch (t1numeric) {
            .bool => @compileError("Cannot cast cfloat to bool"),
            .int => @compileError("Cannot cast cfloat to int"),
            .float => @compileError("Cannot cast cfloat to float"),
            .cfloat => return .{
                .re = @floatCast(value.re),
                .im = @floatCast(value.im),
            },
            .integer => @compileError("Not implemented yet: casting from cfloat to integer"),
            .rational => @compileError("Not implemented yet: casting from cfloat to rational"),
            .real => @compileError("Not implemented yet: casting from cfloat to real"),
            .complex => @compileError("Not implemented yet: casting from cfloat to complex"),
            .expression => @compileError("Not implemented yet: casting from cfloat to expression"),
            .unsupported => unreachable,
        },
        .integer => switch (t1numeric) {
            .bool => @compileError("Not implemented yet: casting from integer to bool"),
            .int => @compileError("Not implemented yet: casting from integer to int"),
            .float => @compileError("Not implemented yet: casting from integer to float"),
            .cfloat => @compileError("Not implemented yet: casting from integer to cfloat"),
            .integer => return value,
            .rational => @compileError("Not implemented yet: casting from integer to rational"),
            .real => @compileError("Not implemented yet: casting from integer to real"),
            .complex => @compileError("Not implemented yet: casting from integer to complex"),
            .expression => @compileError("Not implemented yet: casting from integer to expression"),
            .unsupported => unreachable,
        },
        .rational => switch (t1numeric) {
            .bool => @compileError("Not implemented yet: casting from rational to bool"),
            .int => @compileError("Not implemented yet: casting from rational to int"),
            .float => @compileError("Not implemented yet: casting from rational to float"),
            .cfloat => @compileError("Not implemented yet: casting from rational to cfloat"),
            .integer => @compileError("Not implemented yet: casting from rational to integer"),
            .rational => return value,
            .real => @compileError("Not implemented yet: casting from rational to real"),
            .complex => @compileError("Not implemented yet: casting from rational to complex"),
            .expression => @compileError("Not implemented yet: casting from rational to expression"),
            .unsupported => unreachable,
        },
        .real => switch (t1numeric) {
            .bool => @compileError("Not implemented yet: casting from real to bool"),
            .int => @compileError("Not implemented yet: casting from real to int"),
            .float => @compileError("Not implemented yet: casting from real to float"),
            .cfloat => @compileError("Not implemented yet: casting from real to cfloat"),
            .integer => @compileError("Not implemented yet: casting from real to integer"),
            .rational => @compileError("Not implemented yet: casting from real to rational"),
            .real => return value,
            .complex => @compileError("Not implemented yet: casting from real to complex"),
            .expression => @compileError("Not implemented yet: casting from real to expression"),
            .unsupported => unreachable,
        },
        .complex => switch (t1numeric) {
            .bool => @compileError("Not implemented yet: casting from complex to bool"),
            .int => @compileError("Not implemented yet: casting from complex to int"),
            .float => @compileError("Not implemented yet: casting from complex to float"),
            .cfloat => @compileError("Not implemented yet: casting from complex to cfloat"),
            .integer => @compileError("Not implemented yet: casting from complex to integer"),
            .rational => @compileError("Not implemented yet: casting from complex to rational"),
            .real => @compileError("Not implemented yet: casting from complex to real"),
            .complex => return value,
            .expression => @compileError("Not implemented yet: casting from complex to expression"),
            .unsupported => unreachable,
        },
        .expression => switch (t1numeric) {
            .bool => @compileError("Not implemented yet: casting from expression to bool"),
            .int => @compileError("Not implemented yet: casting from expression to int"),
            .float => @compileError("Not implemented yet: casting from expression to float"),
            .cfloat => @compileError("Not implemented yet: casting from expression to cfloat"),
            .integer => @compileError("Not implemented yet: casting from expression to integer"),
            .rational => @compileError("Not implemented yet: casting from expression to rational"),
            .real => @compileError("Not implemented yet: casting from expression to real"),
            .complex => @compileError("Not implemented yet: casting from expression to complex"),
            .expression => return value,
            .unsupported => unreachable,
        },
        else => unreachable,
    }
}
