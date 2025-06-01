const std = @import("std");
const Array = @import("array.zig").Array;

pub const default_uint = usize;
pub const default_int = isize;
pub const default_float = f64;

const cfloat = @import("cfloat.zig");
const cf16 = @import("cfloat.zig").cf16;
const cf32 = @import("cfloat.zig").cf32;
const cf64 = @import("cfloat.zig").cf64;
const cf80 = @import("cfloat.zig").cf80;
const cf128 = @import("cfloat.zig").cf128;
const comptime_complex = @import("cfloat.zig").comptime_complex;
const integer = @import("integer.zig");
const Integer = integer.Integer;
const rational = @import("rational.zig");
const Rational = rational.Rational;
const real = @import("real.zig");
const Real = real.Real;
const complex = @import("complex.zig");
const Complex = complex.Complex;
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

pub const Order = enum {
    lt,
    eq,
    gt,
};

pub inline fn numericType(comptime T: type) NumericType {
    @setEvalBranchQuota(1000000);

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
            } else if (T == Integer) {
                return .integer;
            } else if (T == Rational) {
                return .rational;
            } else if (T == Real) {
                return .real;
            } else if (T == Complex(Integer) or T == Complex(Rational) or T == Complex(Real)) {
                return .complex;
                //} else if (T == Expression or T == Expression) {
                //    return .expression;
            } else {
                @compileError("Unsupported numeric type: " ++ @typeName(T));
            }
        },
    }
}

/// Checks if the input type is a one-item pointer.
pub fn isPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .one) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is a constant one-item pointer.
pub fn isConstPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .one) return false;

            if (info.is_const) {
                return true;
            } else {
                return false;
            }
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

/// Checks if the input type is an instance of an Array.
pub fn isArray(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .@"struct" => |tinfo| {
            if (tinfo.layout == .@"extern" or tinfo.layout == .@"packed") {
                return false;
            }

            if (!@hasDecl(T, "empty")) return false;

            const N = Array(f64);
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
                        .optional => continue, // Checking if child type is an Array would lead to infinite recursion
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

/// Checks if the input is a one-item pointer.
pub fn isOneItemPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .one) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is of fixed precision.
pub fn isFixedPrecision(comptime T: type) bool {
    switch (numericType(T)) {
        .bool => return true,
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

/// Checks if the input type need allocation for its data.
pub fn needsAllocator(comptime T: type) bool {
    return isArray(T) or isArbitraryPrecision(T);
}

/// Coerces the input types to the smallest type that can represent both types.
/// If at least one of the types is an Array or a slice, the result will be
/// an Array of the coerced type.
pub fn Coerce(comptime K: type, comptime V: type) type {
    if (K == V) {
        return K;
    }

    comptime if (isArray(K)) {
        if (isArray(V)) {
            return Array(Coerce(Numeric(K), Numeric(V)));
        } else if (isSlice(V)) {
            return Array(Coerce(Numeric(K), Numeric(V)));
        } else {
            return Array(Coerce(Numeric(K), V));
        }
    } else if (isSlice(K)) {
        if (isArray(V)) {
            return Array(Coerce(Numeric(K), Numeric(V)));
        } else if (isSlice(V)) { // Should two slices return another slice?
            return Array(Coerce(Numeric(K), Numeric(V)));
        } else {
            return Array(Coerce(Numeric(K), V));
        }
    } else {
        if (isArray(V)) {
            return Array(Coerce(K, Numeric(V)));
        } else if (isSlice(V)) {
            return Array(Coerce(K, Numeric(V)));
        }
        // Else two numeric types
    };

    const T1: type = K;
    const T2: type = V;

    const t1numeric = numericType(T1);
    const t2numeric = numericType(T2);

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
                    const t2info = @typeInfo(Scalar(T2));

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
                .float => {
                    const t1info = @typeInfo(Scalar(T1));
                    const t2info = @typeInfo(T2);

                    if (t1info.float.bits > t2info.float.bits) {
                        return T1;
                    } else {
                        return cfloat.Cfloat(T2);
                    }
                },
                .cfloat => {
                    const t1: T1 = .{ .re = 0, .im = 0 };
                    const t2: T2 = .{ .re = 0, .im = 0 };
                    return cfloat.Cfloat(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
                },
                .integer => return Complex(Rational),
                .rational => return Complex(Rational),
                .real => return Complex(Real),
                .complex => {
                    if (T2 == Complex(Integer)) {
                        return Complex(Rational);
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
            .cfloat => return Complex(Rational),
            .integer => return T1,
            .unsupported => unreachable,
            else => return T2,
        },
        .rational => switch (t2numeric) {
            .bool => return T1,
            .int => return T1,
            .float => return T1,
            .cfloat => return Complex(Rational),
            .integer => return T1,
            .rational => return T1,
            .complex => {
                if (T2 == Complex(Integer)) {
                    return Complex(Rational);
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
            .cfloat => return Complex(Real),
            .integer => return T1,
            .rational => return T1,
            .real => return T1,
            .complex => {
                if (T2 == Complex(Integer)) {
                    return Complex(Real);
                } else if (T2 == Complex(Rational)) {
                    return Complex(Real);
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
                if (T1 == Complex(Integer)) {
                    return Complex(Rational);
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
                return Complex(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
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
    comptime if (isArray(K)) {
        if (isArray(V)) {
            return canCoerce(Numeric(K), Numeric(V));
        } else if (isSlice(V)) {
            return false; // Should csting Array return a slice?
        } else {
            return false; // Only 1 element Arrays can be coerced to scalar types
        }
    } else if (isSlice(K)) {
        if (isArray(V)) {
            return canCoerce(Numeric(K), Numeric(V));
        } else if (isSlice(V)) {
            return canCoerce(Numeric(K), Numeric(V));
        } else {
            return false; // Only 1 element slices can be coerced to scalar types
        }
    } else {
        if (isArray(V)) {
            return canCoerce(K, Numeric(V));
        } else if (isSlice(V)) {
            return canCoerce(K, Numeric(V));
        }
        // Else two numeric types
    };

    const T1: type = K;
    const T2: type = V;

    const t1numeric = numericType(T1);
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
                .integer => T3 = Complex(Rational),
                .rational => T3 = Complex(Rational),
                .real => T3 = Complex(Real),
                .complex => {
                    if (T2 == Complex(Integer)) {
                        T3 = Complex(Rational);
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
            .cfloat => T3 = Complex(Rational),
            .integer => T3 = T1,
            .unsupported => unreachable,
            else => T3 = T2,
        },
        .rational => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = Complex(Rational),
            .integer => T3 = T1,
            .rational => T3 = T1,
            .complex => {
                if (T2 == Complex(Integer)) {
                    T3 = Complex(Rational);
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
            .cfloat => T3 = Complex(Real),
            .integer => T3 = T1,
            .rational => T3 = T1,
            .real => T3 = T1,
            .complex => {
                if (T2 == Complex(Integer)) {
                    T3 = Complex(Real);
                } else if (T2 == Complex(Rational)) {
                    T3 = Complex(Real);
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
                if (T1 == Complex(Integer)) {
                    T3 = Complex(Rational);
                } else if (T1 == Complex(Integer)) {
                    T3 = Complex(Rational);
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
                    T3 = Complex(Coerce(@TypeOf(@field(t1, "re")), @TypeOf(@field(t2, "re"))));
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

/// Coerces the second type to an array of itself if the first type is an
/// Array or a slice. If the first type is not an Array or a slice, it
/// returns the second type as is.
pub fn CoerceToArray(comptime K: type, comptime V: type) type {
    _ = numericType(V);

    if (isArray(K) or isSlice(K)) {
        return Array(V);
    } else {
        return V;
    }
}

/// Checks if t`L` can be cast to `R` safely, meaning that the cast
/// will not result in a runtime panic. For instance, casting signed integers to
/// unsigned integers is not safe, as it can result in a panic if the value is
/// negative.
pub fn canCastSafely(comptime L: type, comptime R: type) bool {
    if (L == R) {
        return true;
    }

    const lnumeric = numericType(L);
    const rnumeric = numericType(R);

    switch (lnumeric) {
        .bool => switch (rnumeric) {
            .bool => return true,
            .int => return true,
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .int => switch (rnumeric) {
            .bool => return true,
            .int => {
                const linfo = @typeInfo(L);
                const rinfo = @typeInfo(R);

                if (linfo.int.signedness == .unsigned and rinfo.int.signedness == .signed) {
                    return false; // Casting unsigned to signed is not safe
                }

                if (linfo.int.signedness == .unsigned and rinfo.int.signedness == .signed) {
                    if (linfo.int.bits >= rinfo.int.bits) {
                        return false; // Casting unsigned to signed is not safe if the unsigned type has more or equal bits
                    }
                }

                if (linfo.int.bits > rinfo.int.bits) {
                    return false; // Casting to a smaller integer type is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .float => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting float to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .cfloat => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting cfloat to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .integer => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting integer to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .rational => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting rational to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .real => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting real to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .complex => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting complex to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .expression => switch (rnumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(R);

                if (rinfo.int.signedness == .unsigned) {
                    return false; // Casting expression to unsigned int is not safe
                }

                return true;
            },
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
            .unsupported => unreachable,
        },
        .unsupported => unreachable,
    }
}

/// Coerces the input type to a floating point type if it is not already a
/// higher range type.
pub fn EnsureFloat(comptime T: type) type {
    switch (numericType(T)) {
        .bool => return default_float,
        .int => return default_float,
        .float => return T,
        .cfloat => return T,
        .integer => return Rational,
        .rational => return T,
        .real => return T,
        .complex => return T,
        else => unreachable,
    }
}

/// Returns the scalar type of a given numeric type, slice, or Array:
/// - Atomic types, such as `u8`, `i16`, `f32`, etc., are returned as-is.
/// - Complex types, such as `cf32`, `cf64`, etc., are returned as their atomic type.
/// - Real arbitrary precision types, such as `Integer`, `Rational`, etc., are returned as-is.
/// - Complex arbitrary precision types are returned as their atomic type.
/// - Slices and Array types are returned as their element type.
pub fn Scalar(comptime T: type) type {
    if (isArray(T)) {
        const t: T = .empty;
        return @typeInfo(@TypeOf(@field(t, "data"))).pointer.child;
    } else if (isSlice(T)) {
        return @typeInfo(T).pointer.child;
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
            Complex(Integer) => return Integer,
            Complex(Rational) => return Rational,
            Complex(Real) => return Real,
            else => unreachable,
        },
        //.expression => return Expression,
        else => unreachable,
    }
}

/// Returns the underlying numeric type of a given numeric type, slice, or
/// Array:
/// - Numeric types, such as floats, integers, complex, etc., are returned as-is.
/// - Slices and Array types are returned as their element type.
pub fn Numeric(comptime T: type) type {
    if (isArray(T)) {
        const t: T = .empty;
        return @typeInfo(@TypeOf(@field(t, "data"))).pointer.child;
    } else if (isSlice(T)) {
        const K = @typeInfo(T).pointer.child;
        _ = numericType(K); // Will raise compile error if K is not a supported numeric type
        return K;
    }

    _ = numericType(T); // Will raise compile error if T is not a supported numeric type

    return T;
}

/// Casts a value of any numeric type to any fixed precision numeric type. It is
/// a more concise version of `cast` that does not require any options and cannot
/// err.
pub inline fn scast(
    comptime T: type,
    value: anytype,
) T {
    const I: type = @TypeOf(value);
    const O: type = T;

    comptime if (!isFixedPrecision(O))
        @compileError("Expected a fixed precision type, but got " ++ @typeName(O));

    if (I == O) {
        return value;
    }

    const onumeric = numericType(O);
    const inumeric = numericType(I);

    switch (inumeric) {
        .bool => switch (onumeric) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1 else 0,
            .cfloat => return .{
                .re = if (value) 1 else 0,
                .im = 0,
            },
            else => unreachable,
        },
        .int => switch (onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0,
            },
            else => unreachable,
        },
        .float => switch (onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .cfloat => return if (I == Scalar(O)) .{
                .re = value,
                .im = 0,
            } else .{
                .re = @floatCast(value),
                .im = 0,
            },
            else => unreachable,
        },
        .cfloat => switch (onumeric) {
            .bool => return if (value.re != 0 or value.im != 0) true else false,
            .int => return @intFromFloat(value.re),
            .float => return if (Scalar(I) == O) value.re else @floatCast(value.re),
            .cfloat => return .{
                .re = @floatCast(value.re),
                .im = @floatCast(value.im),
            },
            else => unreachable,
        },
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
        allocator: if (needsAllocator(T)) std.mem.Allocator else void = {},
        copy: bool = false,
    },
) !T {
    const I: type = @TypeOf(value);
    const O: type = T;

    if (I == O) {
        switch (numericType(O)) {
            .bool, .int, .float, .cfloat => return value,
            .integer => if (options.copy) {
                // return try integer.copy(options.allocator, value);
                return value;
            } else {
                return value;
            },
            .rational => if (options.copy) {
                // return try rational.copy(options.allocator, value);
                return value;
            } else {
                return value;
            },
            .real => if (options.copy) {
                // return try real.copy(options.allocator, value);
                return value;
            } else {
                return value;
            },
            .complex => if (options.copy) {
                // return try complex.copy(options.allocator, value);
                return value;
            } else {
                return value;
            },
            .expression => if (options.copy) {
                // return try expression.copy(options.allocator, value);
                return value;
            } else {
                return value;
            },
            .unsupported => unreachable,
        }
        return value;
    }

    const onumeric = numericType(O);
    const inumeric = numericType(I);

    switch (inumeric) {
        .bool => switch (onumeric) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1 else 0,
            .cfloat => return .{
                .re = if (value) 1 else 0,
                .im = 0,
            },
            .integer => @compileError("Not implemented yet: casting from bool to integer"),
            .rational => @compileError("Not implemented yet: casting from bool to rational"),
            .real => @compileError("Not implemented yet: casting from bool to real"),
            .complex => @compileError("Not implemented yet: casting from bool to complex"),
            .expression => @compileError("Not implemented yet: casting from bool to expression"),
            .unsupported => unreachable,
        },
        .int => switch (onumeric) {
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
        .float => switch (onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intFromFloat(value),
            .float => return @floatCast(value),
            .cfloat => return if (I == Scalar(O)) .{
                .re = value,
                .im = 0,
            } else .{
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
        .cfloat => switch (onumeric) {
            .bool => return if (value.re != 0 or value.im != 0) true else false,
            .int => return @intFromFloat(value.re),
            .float => return if (Scalar(I) == O) value.re else @floatCast(value.re),
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
        .integer => switch (onumeric) {
            .bool => @compileError("Not implemented yet: casting from integer to bool"),
            .int => @compileError("Not implemented yet: casting from integer to int"),
            .float => @compileError("Not implemented yet: casting from integer to float"),
            .cfloat => @compileError("Not implemented yet: casting from integer to cfloat"),
            .integer => unreachable,
            .rational => @compileError("Not implemented yet: casting from integer to rational"),
            .real => @compileError("Not implemented yet: casting from integer to real"),
            .complex => @compileError("Not implemented yet: casting from integer to complex"),
            .expression => @compileError("Not implemented yet: casting from integer to expression"),
            .unsupported => unreachable,
        },
        .rational => switch (onumeric) {
            .bool => @compileError("Not implemented yet: casting from rational to bool"),
            .int => @compileError("Not implemented yet: casting from rational to int"),
            .float => @compileError("Not implemented yet: casting from rational to float"),
            .cfloat => @compileError("Not implemented yet: casting from rational to cfloat"),
            .integer => @compileError("Not implemented yet: casting from rational to integer"),
            .rational => unreachable,
            .real => @compileError("Not implemented yet: casting from rational to real"),
            .complex => @compileError("Not implemented yet: casting from rational to complex"),
            .expression => @compileError("Not implemented yet: casting from rational to expression"),
            .unsupported => unreachable,
        },
        .real => switch (onumeric) {
            .bool => @compileError("Not implemented yet: casting from real to bool"),
            .int => @compileError("Not implemented yet: casting from real to int"),
            .float => @compileError("Not implemented yet: casting from real to float"),
            .cfloat => @compileError("Not implemented yet: casting from real to cfloat"),
            .integer => @compileError("Not implemented yet: casting from real to integer"),
            .rational => @compileError("Not implemented yet: casting from real to rational"),
            .real => unreachable,
            .complex => @compileError("Not implemented yet: casting from real to complex"),
            .expression => @compileError("Not implemented yet: casting from real to expression"),
            .unsupported => unreachable,
        },
        .complex => switch (onumeric) {
            .bool => @compileError("Not implemented yet: casting from complex to bool"),
            .int => @compileError("Not implemented yet: casting from complex to int"),
            .float => @compileError("Not implemented yet: casting from complex to float"),
            .cfloat => @compileError("Not implemented yet: casting from complex to cfloat"),
            .integer => @compileError("Not implemented yet: casting from complex to integer"),
            .rational => @compileError("Not implemented yet: casting from complex to rational"),
            .real => @compileError("Not implemented yet: casting from complex to real"),
            .complex => return value, // Will have to check the type held by the Complex
            .expression => @compileError("Not implemented yet: casting from complex to expression"),
            .unsupported => unreachable,
        },
        .expression => switch (onumeric) {
            .bool => @compileError("Not implemented yet: casting from expression to bool"),
            .int => @compileError("Not implemented yet: casting from expression to int"),
            .float => @compileError("Not implemented yet: casting from expression to float"),
            .cfloat => @compileError("Not implemented yet: casting from expression to cfloat"),
            .integer => @compileError("Not implemented yet: casting from expression to integer"),
            .rational => @compileError("Not implemented yet: casting from expression to rational"),
            .real => @compileError("Not implemented yet: casting from expression to real"),
            .complex => @compileError("Not implemented yet: casting from expression to complex"),
            .expression => unreachable,
            .unsupported => unreachable,
        },
        else => unreachable,
    }
}

/// Returns the pointer child type of a given pointer type.
pub fn Child(comptime T: type) type {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            return info.child;
        },
        else => @compileError("Expected a pointer type, but got " ++ @typeName(T)),
    }
}
