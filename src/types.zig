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

/// `NumericType` is an enum that represents the different numeric types
/// supported by the library. It is used to categorize types based on their
/// properties and capabilities, such as whether they are integers, floats,
/// complex numbers, etc.
///
/// This enum is used in various places in the library to determine how to
/// handle different types of numeric data. It allows for type checking and
/// coercion between different numeric types, ensuring that operations are
/// performed correctly and efficiently.
///
/// Values
/// ------
/// - `bool`: Represents the boolean type (`bool`).
/// - `int`: Represents integer types:
///   - `usize`, `u8`, `u16`, `u32`, `u64`, `u128`
///   - `isize`, `i8`, `i16`, `i32`, `i64`, `i128`
///   - `comptime_int`
/// - `float`: Represents floating-point types:
///   - `f16`, `f32`, `f64`, `f80`, `f128`
///   - `comptime_float`
/// - `cfloat`: Represents complex floating-point types:
///   - `cf16`, `cf32`, `cf64`, `cf80`, `cf128`
///   - `comptime_complex`
/// - `integer`: Represents arbitrary precision integer type (`Integer`).
/// - `rational`: Represents arbitrary precision rational type (`Rational`).
/// - `real`: Represents arbitrary precision real type (`Real`).
/// - `complex`: Represents complex arbitrary precision types:
///   - `Complex(Integer)`
///   - `Complex(Rational)`
///   - `Complex(Real)`
/// - `expression`: Represents symbolic expressions (Expression).
pub const NumericType = enum {
    bool,
    int,
    float,
    cfloat,
    integer,
    rational,
    real,
    complex,
    expression,
};

pub const Order = enum {
    lt,
    eq,
    gt,
};

pub const useless_allocator: std.mem.Allocator = .{
    .ptr = undefined,
    .vtable = &vtable,
};

pub const vtable: std.mem.Allocator.VTable = .{
    .alloc = alloc,
    .resize = resize,
    .remap = remap,
    .free = free,
};

fn alloc(context: *anyopaque, len: usize, alignment: std.mem.Alignment, ra: usize) ?[*]u8 {
    _ = context;
    _ = len;
    _ = alignment;
    _ = ra;

    return null;
}

fn resize(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ra: usize) bool {
    _ = context;
    _ = memory;
    _ = alignment;
    _ = new_len;
    _ = ra;

    return false;
}

fn remap(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, new_len: usize, ra: usize) ?[*]u8 {
    _ = context;
    _ = memory;
    _ = alignment;
    _ = new_len;
    _ = ra;

    return null;
}

fn free(context: *anyopaque, memory: []u8, alignment: std.mem.Alignment, ra: usize) void {
    _ = context;
    _ = memory;
    _ = alignment;
    _ = ra;

    return;
}

/// Checks the the input type `T` and returns the corresponding `NumericType`.
///
/// Checks that the input type is a supported numeric type and returns the
/// corresponding `NumericType` enum value. If the type is not supported, it
/// will raise a compile error.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `NumericType`: The corresponding `NumericType` enum value.
///
/// Raises
/// ------
/// `@compileError`: If the type is not supported.
pub inline fn numericType(comptime T: type) NumericType {
    // Without inline functions calling this fail miserably. I have no idea why.

    @setEvalBranchQuota(1000000);

    switch (@typeInfo(T)) {
        .bool => return .bool,
        .int, .comptime_int => {
            if (T != u8 and T != u16 and T != u32 and T != u64 and T != u128 and T != usize and T != c_uint and T != i8 and T != i16 and T != i32 and T != i64 and T != i128 and T != isize and T != c_int and T != comptime_int)
                @compileError("Unsupported integer type: " ++ @typeName(T));

            return .int;
        },
        .float, .comptime_float => return .float,
        else => {
            if (T == cf16 or T == cf32 or T == cf64 or T == cf80 or T == cf128 or T == comptime_complex or T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_float)) {
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

pub fn isNumeric(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .bool => return true,
        .int, .comptime_int => {
            if (T != u8 and T != u16 and T != u32 and T != u64 and T != u128 and T != usize and T != c_uint and T != i8 and T != i16 and T != i32 and T != i64 and T != i128 and T != isize and T != c_int and T != comptime_int)
                return false;

            return true;
        },
        .float, .comptime_float => return true,
        else => {
            if (T == cf16 or T == cf32 or T == cf64 or T == cf80 or T == cf128 or T == comptime_complex or T == std.math.Complex(f16) or T == std.math.Complex(f32) or T == std.math.Complex(f64) or T == std.math.Complex(f80) or T == std.math.Complex(f128) or T == std.math.Complex(comptime_float)) {
                return true;
            } else if (T == Integer) {
                return true;
            } else if (T == Rational) {
                return true;
            } else if (T == Real) {
                return true;
            } else if (T == Complex(Integer) or T == Complex(Rational) or T == Complex(Real)) {
                return true;
                //} else if (T == Expression or T == Expression) {
                //    return .expression;
            } else {
                return false;
            }
        },
    }
}

/// Checks if the input type is a one-item pointer.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a one-item pointer, `false` otherwise.
pub fn isPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .one and
                info.size != .c) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is a many-item pointer.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a many-item pointer, `false` otherwise.
pub fn isManyPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .many and
                info.size != .c) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is a constant pointer. Works for one-item pointers,
/// many-item pointers, and slices.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a constant one-item pointer, `false`
/// otherwise.
pub fn isConstPointer(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.is_const) {
                return true;
            } else {
                return false;
            }
        },
        else => return false,
    }
}

/// Checks if the input type is a slice.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is a slice, `false` otherwise.
pub fn isSlice(comptime T: type) bool {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            if (info.size != .slice) return false;

            return true;
        },
        else => return false,
    }
}

/// Checks if the input type is an instance of an `Array`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check.
///
/// Returns
/// -------
/// `bool`: `true` if the type is an `Array`, `false` otherwise.
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

/// Checks if the input numeric type is of fixed precision.
///
/// This function checks if the input type is a fixed precision numeric type,
/// which includes types labeled as `bool`, `int`, `float`, and `cfloat`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is of fixed precision, `false` otherwise.
///
/// Raises
/// ------
/// `@compileError`: If the type is not a supported numeric type.
pub fn isFixedPrecision(comptime T: type) bool {
    switch (numericType(T)) {
        .bool => return true,
        .int => return true,
        .float => return true,
        .cfloat => return true,
        else => return false,
    }
}

/// Checks if the input numeric type is of arbitrary precision.
///
/// This function checks if the input type is an arbitrary precision numeric
/// type, which includes types labeled as `integer`, `rational`, `real`,
/// `complex`, and `expression`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is of arbitrary precision, `false` otherwise.
///
/// Raises
/// ------
/// `@compileError`: If the type is not a supported numeric type.
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

/// Checks if the input type is complex.
///
/// This function checks if the input type is a complex numeric type, which
/// includes types labeled as `cfloat`, `complex`, and `expression`.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type is complex, `false` otherwise.
///
/// Raises
/// ------
/// `@compileError`: If the type is not a supported numeric type.
pub fn isComplex(comptime T: type) bool {
    switch (numericType(T)) {
        .cfloat => return true,
        .complex => return true,
        else => return false,
    }
}

/// Checks if the input type needs allocation for its data.
///
/// This function checks if the input type requires an allocator for its data,
/// which inludes `Array`s and arbitrary precision types.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to check. Must be a supported numeric type or
/// an `Array`.
///
/// Returns
/// -------
/// `bool`: `true` if the type needs an allocator, `false` otherwise.
///
/// Raises
/// ------
/// `@compileError`: If the type is neither an `Array` nor a supported numeric
/// type.
pub fn needsAllocator(comptime T: type) bool {
    return isArray(T) or isArbitraryPrecision(T);
}

/// Coerces the input types to the smallest type that can represent both types.
///
/// This function takes two types `X` and `Y` and returns the smallest type that
/// can represent both types without loss of information. If at least one of the
/// types is an `Array` or a slice, the result will be an `Array` of the coerced
/// `Numeric` types.
///
/// Parameters
/// ----------
/// comptime X (`type`): The first type to coerce. Must be a supported numeric
/// type, an `Array`, or a slice.
///
/// comptime Y (`type`): The second type to coerce. Must be a supported numeric
/// type, an `Array`, or a slice.
///
/// Returns
/// -------
/// `type`: The coerced type that can represent both `X` and `Y`.
///
/// Raises
/// ------
/// `@compileError`: If the types cannot be coerced, or if either type is not an
/// `Array`, a slice, or a supported numeric type.
pub fn Coerce(comptime X: type, comptime Y: type) type {
    if (X == Y) {
        return X;
    }

    comptime if (isArray(X)) {
        if (isArray(Y)) {
            return Array(Coerce(Numeric(X), Numeric(Y)));
        } else if (isSlice(Y)) {
            return Array(Coerce(Numeric(X), Numeric(Y)));
        } else {
            return Array(Coerce(Numeric(X), Y));
        }
    } else if (isSlice(X)) {
        if (isArray(Y)) {
            return Array(Coerce(Numeric(X), Numeric(Y)));
        } else if (isSlice(Y)) { // Should two slices return another slice?
            return Array(Coerce(Numeric(X), Numeric(Y)));
        } else {
            return Array(Coerce(Numeric(X), Y));
        }
    } else {
        if (isArray(Y)) {
            return Array(Coerce(X, Numeric(Y)));
        } else if (isSlice(Y)) {
            return Array(Coerce(X, Numeric(Y)));
        }
    };

    // Else two numeric types
    const xnumeric = numericType(X);
    const ynumeric = numericType(Y);

    switch (xnumeric) {
        .bool => switch (ynumeric) {
            .bool => return bool,
            else => return Y,
        },
        .int => switch (ynumeric) {
            .bool => return X,
            .int => {
                if (X == comptime_int) {
                    return Y;
                } else if (Y == comptime_int) {
                    return X;
                }

                comptime var xinfo = @typeInfo(X);
                comptime var yinfo = @typeInfo(Y);

                if (xinfo.int.signedness == .unsigned) {
                    if (yinfo.int.signedness == .unsigned) {
                        if (xinfo.int.bits > yinfo.int.bits) {
                            return X;
                        } else {
                            return Y;
                        }
                    } else {
                        if (xinfo.int.bits > yinfo.int.bits) {
                            if (xinfo.int.bits == 128) {
                                @compileError("Cannot coerce " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " as " ++ @typeName(X) ++ " already has the maximum amount of bits.");
                            } else {
                                yinfo.int.bits = xinfo.int.bits * 2;
                                return @Type(yinfo);
                            }
                        } else if (xinfo.int.bits == yinfo.int.bits) {
                            yinfo.int.bits *= 2;
                            return @Type(yinfo);
                        } else {
                            return Y;
                        }
                    }
                } else {
                    if (yinfo.int.signedness == .unsigned) {
                        if (yinfo.int.bits > xinfo.int.bits) {
                            if (yinfo.int.bits == 128) {
                                @compileError("Cannot coerce " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " as " ++ @typeName(Y) ++ " already has the maximum amount of bits.");
                            } else {
                                xinfo.int.bits = yinfo.int.bits * 2;
                                return @Type(xinfo);
                            }
                        } else if (xinfo.int.bits == yinfo.int.bits) {
                            xinfo.int.bits *= 2;
                            return @Type(xinfo);
                        } else {
                            return X;
                        }
                    } else {
                        if (xinfo.int.bits > yinfo.int.bits) {
                            return X;
                        } else {
                            return Y;
                        }
                    }
                }
            },
            .float => {
                if (Y == comptime_float) {
                    return EnsureFloat(X);
                }

                return Y;
            },
            // .cfloat -> internal float management
            else => return Y,
        },
        .float => {
            switch (ynumeric) {
                .bool, .int => {
                    if (X == comptime_float) {
                        return EnsureFloat(Y);
                    }

                    return X;
                },
                .float => {
                    if (X == comptime_float) {
                        return Y;
                    } else if (Y == comptime_float) {
                        return X;
                    }

                    const xinfo = @typeInfo(X);
                    const yinfo = @typeInfo(Y);

                    if (xinfo.float.bits > yinfo.float.bits) {
                        return X;
                    } else {
                        return Y;
                    }
                },
                .cfloat => {
                    const xinfo = @typeInfo(X);
                    const yinfo = @typeInfo(Scalar(Y));

                    if (xinfo.float.bits > yinfo.float.bits) {
                        return cfloat.Cfloat(X);
                    } else {
                        return Y;
                    }
                },
                else => return Y,
            }
        },
        .cfloat => {
            switch (ynumeric) {
                .bool => return X,
                .int => return X,
                .float => {
                    const xinfo = @typeInfo(Scalar(X));
                    const yinfo = @typeInfo(Y);

                    if (xinfo.float.bits > yinfo.float.bits) {
                        return X;
                    } else {
                        return cfloat.Cfloat(Y);
                    }
                },
                .cfloat => {
                    const x: X = .{ .re = 0, .im = 0 };
                    const y: Y = .{ .re = 0, .im = 0 };
                    return cfloat.Cfloat(Coerce(@TypeOf(@field(x, "re")), @TypeOf(@field(y, "re"))));
                },
                .integer => return Complex(Rational),
                .rational => return Complex(Rational),
                .real => return Complex(Real),
                .complex => {
                    if (Y == Complex(Integer)) {
                        return Complex(Rational);
                    } else {
                        return Y;
                    }
                },
                else => return Y,
            }
        },
        .integer => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return Complex(Rational),
            .integer => return X,
            else => return Y,
        },
        .rational => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return Complex(Rational),
            .integer => return X,
            .rational => return X,
            .complex => {
                if (Y == Complex(Integer)) {
                    return Complex(Rational);
                } else {
                    return Y;
                }
            },
            else => return Y,
        },
        .real => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return Complex(Real),
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => {
                if (Y == Complex(Integer)) {
                    return Complex(Real);
                } else if (Y == Complex(Rational)) {
                    return Complex(Real);
                } else {
                    return Y;
                }
            },
            else => return Y,
        },
        .complex => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => {
                if (X == Complex(Integer)) {
                    return Complex(Rational);
                } else {
                    return X;
                }
            },
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => {
                const x: X = .empty;
                const y: Y = .empty;
                return Complex(Coerce(@TypeOf(@field(x, "re")), @TypeOf(@field(y, "re"))));
            },
            .expression => return Y,
        },
        .expression => switch (ynumeric) {
            .bool => return X,
            .int => return X,
            .float => return X,
            .cfloat => return X,
            .integer => return X,
            .rational => return X,
            .real => return X,
            .complex => return X,
            .expression => return X,
        },
    }
}

/// Checks if if `K` can be coerced to `V` without loss of information. This is
/// a more flexible version of `Coerce`, as it does not require `V` to be to the
/// smallest type that can represent both types. The only requirement is that
/// `V` can represent all values of the first two types.
///
/// Unused right now, kept for possible future use.
pub fn canCoerce(comptime K: type, comptime V: type) bool {
    comptime if (isArray(K)) {
        if (isArray(V)) {
            return canCoerce(Numeric(K), Numeric(V));
        } else if (isSlice(V)) {
            return false; // Should casting Array return a slice?
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
                else => T3 = T2,
            }
        },
        .integer => switch (t2numeric) {
            .bool => T3 = T1,
            .int => T3 = T1,
            .float => T3 = T1,
            .cfloat => T3 = Complex(Rational),
            .integer => T3 = T1,
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
    }

    return T2 == T3;
}

/// Coerces the second type to an array of itself if the first type is an
/// `Array` or a slice. If the first type is not an Array or a slice, it
/// returns the second type as is.
/// Coerces `Y` to `Array(Y)` if `X` is an `Array` or a slice, otherwise it
/// returns `Y` as is.
///
/// This function is useful for ensuring that the second type is always an
/// `Array` when the first type is an `Array` or a slice, allowing for
/// consistent handling of array-like types.
///
/// Parameters
/// ----------
/// comptime X (`type`): The type to check. Must be a supported numeric type, an
/// `Array`, or a slice.
///
/// comptime Y (`type`): The type to coerce. Must be a supported numeric type.
///
/// Returns
/// -------
/// `type`: The coerced type, which will be `Array(Y)` if `X` is an `Array` or a
/// slice, or `Y` otherwise.
///
/// Raises
/// ------
/// `@compileError`: If `Y` is not a supported numeric type.
pub fn CoerceToArray(comptime X: type, comptime Y: type) type {
    _ = numericType(Y);

    if (isArray(X) or isSlice(X)) {
        return Array(Y);
    } else {
        return Y;
    }
}

/// Checks if `X` can be safely cast to `Y`, meaning that the cast will not
/// result in a runtime panic.
///
/// This function checks if the type `X` can be safely cast to type `Y` without
/// losing information or causing a runtime panic. It considers various numeric
/// types and their properties, such as signedness and bit width, to determine
/// if the cast is safe.
///
/// For example, casting a signed integer to an unsigned integer is not
/// considered safe, as it may lead to unexpected results if the value is
/// negative.
///
/// Parameters
/// ----------
/// comptime X (`type`): The type to check if it can be cast safely. Must be a
/// supported numeric type.
///
/// comptime Y (`type`): The type to check if `X` can be cast to. Must be a
/// supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the cast is safe, `false` otherwise.
///
/// Raises
/// ------
/// `@compileError`: If either type is not a supported numeric type.
pub fn canCastSafely(comptime X: type, comptime Y: type) bool {
    if (X == Y) {
        return true;
    }

    const xnumeric = numericType(X);
    const ynumeric = numericType(Y);

    switch (xnumeric) {
        .bool => switch (ynumeric) {
            .bool => return true,
            .int => return true,
            .float => return true,
            .cfloat => return true,
            .integer => return true,
            .rational => return true,
            .real => return true,
            .complex => return true,
            .expression => return true,
        },
        .int => switch (ynumeric) {
            .bool => return true,
            .int => {
                const linfo = @typeInfo(X);
                const rinfo = @typeInfo(Y);

                if (linfo.int.signedness == .snsigned and rinfo.int.signedness == .unsigned) {
                    // Casting signed to unsigned is not safe
                    return false;
                }

                if (linfo.int.signedness == .unsigned and rinfo.int.signedness == .signed) {
                    if (linfo.int.bits >= rinfo.int.bits) {
                        // Casting unsigned to signed is not safe if the unsigned type has more or equal bits
                        return false;
                    }
                }

                if (linfo.int.bits > rinfo.int.bits) {
                    // Casting to a smaller integer type is not safe
                    return false;
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
        },
        .float => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
        .cfloat => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
        .integer => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
        .rational => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
        .real => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
        .complex => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
        .expression => switch (ynumeric) {
            .bool => return true,
            .int => {
                const rinfo = @typeInfo(Y);

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
        },
    }
}

/// Coerces the input type to a floating point type if it is not already a
/// higher range type.
pub fn EnsureFloat(comptime T: type) type {
    if (isArray(T) or isSlice(T))
        return Array(EnsureFloat(Numeric(T)));

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

/// Returns the scalar type of a given numeric type, slice, or `Array`.
///
/// This function returns the scalar type of a given numeric type, slice, or
/// `Array`. If the input type is an `Array` or a slice, it returns the element
/// type. If the input type is a numeric type, it returns the type itself,
/// unless it is a complex type, in which case it returns the scalar type of the
/// complex type.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to get the scalar type of. Must be a supported
/// numeric type, slice, or `Array`.
///
/// Returns
/// -------
/// `type`: The scalar type of the input type.
///
/// Raises
/// ------
/// `@compileError`: If the type is not a supported numeric type, slice, or
/// `Array`.
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
/// `Array`.
///
/// This function returns the underlying numeric type of a given numeric type,
/// slice, or `Array`. If the input type is an `Array` or a slice, it returns
/// the element type. If the input type is a numeric type, it returns the type
/// itself.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to get the numeric type of. Must be a
/// supported numeric type, slice, or `Array`.
///
/// Returns
/// -------
/// `type`: The underlying numeric type of the input type.
///
/// Raises
/// ------
/// `@compileError`: If the type is not a supported numeric type, slice, or
/// `Array`.
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

/// Casts a value of any numeric type to any fixed precision numeric type.
///
/// It is a more concise version of `cast` that does not require any options and
/// cannot return an error.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to cast to. Must be a fixed precision numeric
/// type.
///
/// value (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`): The value to cast.
///
/// Returns
/// -------
/// `T`: The value casted to the type `T`.
///
/// Raises
/// ------
/// `@compileError`: If the type `T` is not a fixed precision numeric type or if
/// the type of `value` is not a supported numeric type.
///
/// Notes
/// -----
/// This function does not check if the cast is safe, which could lead to
/// runtime panics if the value cannot be represented in the target type.
pub inline fn scast(
    comptime T: type,
    value: anytype,
) T {
    const I: type = @TypeOf(value);
    const O: type = T;

    comptime if (!isFixedPrecision(O))
        @compileError("Expected a fixed precision type, but got " ++ @typeName(O));

    if (comptime I == O) {
        return value;
    }

    const onumeric = numericType(O);
    const inumeric = numericType(I);

    switch (comptime inumeric) {
        .bool => switch (comptime onumeric) {
            .bool => unreachable,
            .int => return if (value) 1 else 0,
            .float => return if (value) 1 else 0,
            .cfloat => return .{
                .re = if (value) 1 else 0,
                .im = 0,
            },
            else => unreachable,
        },
        .int => switch (comptime onumeric) {
            .bool => return if (value != 0) true else false,
            .int => return @intCast(value),
            .float => return @floatFromInt(value),
            .cfloat => return .{
                .re = @floatFromInt(value),
                .im = 0,
            },
            else => unreachable,
        },
        .float => switch (comptime onumeric) {
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
        .cfloat => switch (comptime onumeric) {
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

/// Casts a value of any numeric type to any other numeric type.
///
/// Parameters
/// ----------
/// comptime `T` (`type`): The type to cast to. Must be a supported numeric type.
///
/// `value` (`bool`, `int`, `float`, `cfloat`, `integer`, `rational`, `real`,
/// `complex` or `expression`):The value to cast.
///
/// `ctx` (`struct`): The context for the cast operation.
///
/// Returns
/// -------
/// `T`: The value casted to the type `T`.
///
/// Errors
/// ------
/// `std.mem.Allocator.Error.OutOfMemory`: If the allocator fails to allocate
/// memory for the output value.
///
/// Raises
/// ------
/// `@compileError`: If the type `T` is not a supported numeric type or if the
/// type of `value` is not a supported numeric type.
///
/// Notes
/// -----
/// This function does not check if the cast is safe, which could lead to
/// runtime panics.
pub inline fn cast(
    comptime T: type,
    value: anytype,
    ctx: anytype,
) !T {
    const I: type = @TypeOf(value);
    const O: type = T;

    comptime if (isArbitraryPrecision(O)) {
        if (I == O) {
            validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                    .copy = .{ .type = bool, .required = false },
                },
            );
        } else {
            validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        }
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    if (I == O) {
        switch (comptime numericType(O)) {
            .bool, .int, .float, .cfloat => return value,
            .integer => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return integer.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .rational => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return rational.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .real => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return real.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .complex => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return complex.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
            .expression => if (getFieldOrDefault(ctx, "copy", bool, false)) {
                // return expression.copy(ctx.allocator, value);
                return value;
            } else {
                return value;
            },
        }
        return value;
    }

    const onumeric = numericType(O);
    const inumeric = comptime numericType(I);

    switch (inumeric) {
        .bool => switch (comptime onumeric) {
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
        },
        .int => switch (comptime onumeric) {
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
        },
        .float => switch (comptime onumeric) {
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
        },
        .cfloat => switch (comptime onumeric) {
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
        },
        .integer => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from integer to bool"),
            .int => @compileError("Not implemented yet: casting from integer to int"),
            .float => @compileError("Not implemented yet: casting from integer to float"),
            .cfloat => @compileError("Not implemented yet: casting from integer to cfloat"),
            .integer => unreachable,
            .rational => @compileError("Not implemented yet: casting from integer to rational"),
            .real => @compileError("Not implemented yet: casting from integer to real"),
            .complex => @compileError("Not implemented yet: casting from integer to complex"),
            .expression => @compileError("Not implemented yet: casting from integer to expression"),
        },
        .rational => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from rational to bool"),
            .int => @compileError("Not implemented yet: casting from rational to int"),
            .float => @compileError("Not implemented yet: casting from rational to float"),
            .cfloat => @compileError("Not implemented yet: casting from rational to cfloat"),
            .integer => @compileError("Not implemented yet: casting from rational to integer"),
            .rational => unreachable,
            .real => @compileError("Not implemented yet: casting from rational to real"),
            .complex => @compileError("Not implemented yet: casting from rational to complex"),
            .expression => @compileError("Not implemented yet: casting from rational to expression"),
        },
        .real => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from real to bool"),
            .int => @compileError("Not implemented yet: casting from real to int"),
            .float => @compileError("Not implemented yet: casting from real to float"),
            .cfloat => @compileError("Not implemented yet: casting from real to cfloat"),
            .integer => @compileError("Not implemented yet: casting from real to integer"),
            .rational => @compileError("Not implemented yet: casting from real to rational"),
            .real => unreachable,
            .complex => @compileError("Not implemented yet: casting from real to complex"),
            .expression => @compileError("Not implemented yet: casting from real to expression"),
        },
        .complex => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from complex to bool"),
            .int => @compileError("Not implemented yet: casting from complex to int"),
            .float => @compileError("Not implemented yet: casting from complex to float"),
            .cfloat => @compileError("Not implemented yet: casting from complex to cfloat"),
            .integer => @compileError("Not implemented yet: casting from complex to integer"),
            .rational => @compileError("Not implemented yet: casting from complex to rational"),
            .real => @compileError("Not implemented yet: casting from complex to real"),
            .complex => return value, // Will have to check the type held by the Complex
            .expression => @compileError("Not implemented yet: casting from complex to expression"),
        },
        .expression => switch (comptime onumeric) {
            .bool => @compileError("Not implemented yet: casting from expression to bool"),
            .int => @compileError("Not implemented yet: casting from expression to int"),
            .float => @compileError("Not implemented yet: casting from expression to float"),
            .cfloat => @compileError("Not implemented yet: casting from expression to cfloat"),
            .integer => @compileError("Not implemented yet: casting from expression to integer"),
            .rational => @compileError("Not implemented yet: casting from expression to rational"),
            .real => @compileError("Not implemented yet: casting from expression to real"),
            .complex => @compileError("Not implemented yet: casting from expression to complex"),
            .expression => unreachable,
        },
    }
}

/// Returns the pointer child type of a given pointer type.
///
/// Parameters
/// ----------
/// comptime T (`type`): The type to get the child type of. Must be a pointer
/// type.
///
/// Returns
/// -------
/// `type`: The child type of the input pointer type.
///
/// Raises
/// ------
/// `@compileError`: If the type is not a pointer type.
pub fn Child(comptime T: type) type {
    switch (@typeInfo(T)) {
        .pointer => |info| {
            return info.child;
        },
        else => @compileError("Expected a pointer type, but got " ++ @typeName(T)),
    }
}

/// Returns the return type of a function with one parameter when called with a
/// value of type `X`.
///
/// This function is useful when the return type of the function depends on
/// the type of the parameter passed to it. The function may have a second
/// parameter that is a context struct.
pub fn ReturnType1(comptime op: anytype, comptime X: type) type {
    const opinfo = @typeInfo(@TypeOf(op));

    comptime if (opinfo.@"fn".params.len != 1 and opinfo.@"fn".params.len != 2)
        @compileError("ReturnType1: op must be a function with one parameter, or two if the second is an context struct");

    const val: X = if (isArbitraryPrecision(X))
        .empty
    else if (isComplex(X))
        .{ .re = 1, .im = 1 }
    else
        1;

    const result_type: type = if (opinfo.@"fn".params.len == 1)
        @TypeOf(op(val))
    else // We must pass an allocator, although it is not used
        if (needsAllocator(X))
            @TypeOf(op(val, .{ .allocator = useless_allocator }))
        else
            @TypeOf(op(val, .{}));

    const resinfo = @typeInfo(result_type);
    switch (resinfo) {
        .error_union => return resinfo.error_union.payload,
        else => return result_type,
    }
}

/// Returns the return type of a function with two parameters when called with
/// values of types `X` and `Y`.
///
/// This function is useful when the return type of the function depends on
/// the types of the parameters passed to it. The function may have a third
/// parameter that is a context struct.
pub fn ReturnType2(comptime op: anytype, comptime X: type, comptime Y: type) type {
    const opinfo = @typeInfo(@TypeOf(op));

    comptime if (opinfo.@"fn".params.len != 2 and opinfo.@"fn".params.len != 3)
        @compileError("ReturnType2: op must be a function with two parameters, or three if the third is an context struct");

    const val1: X = if (isArbitraryPrecision(X))
        .empty
    else if (isComplex(X))
        .{ .re = 1, .im = 1 }
    else
        1;
    const val2: Y = if (isArbitraryPrecision(Y))
        .empty
    else if (isComplex(Y))
        .{ .re = 1, .im = 1 }
    else
        1;

    const result_type: type = if (opinfo.@"fn".params.len == 2)
        @TypeOf(op(val1, val2))
    else // We must pass an allocator, although it is not used
        if (needsAllocator(X) or needsAllocator(Y))
            @TypeOf(op(val1, val2, .{ .allocator = useless_allocator }))
        else
            @TypeOf(op(val1, val2, .{}));

    const resinfo = @typeInfo(result_type);
    switch (resinfo) {
        .error_union => return resinfo.error_union.payload,
        else => return result_type,
    }
}

/// Validates a context struct against a specification struct.
///
/// This function checks that the context struct matches the specification
/// struct, ensuring that all required fields are present and have the correct
/// types. It also checks for unexpected fields in the context struct.
///
/// Parameters
/// ----------
/// ctx (`anytype`): The context struct to validate. Must be a struct type.
/// spec (`anytype`): The specification struct that defines the expected fields
/// and their types. Must have the following structure:
/// ```zig
/// .{
///     .field_name = .{ .type = type, .required = bool },
///     ...
/// }
/// ```
/// where `field_name` is the name of the field, `type` is the expected type of
/// the field, and `required` is a boolean indicating whether the field is
/// required or not.
///
/// Returns
/// -------
/// `void`: If the context struct is valid according to the specification.
///
/// Raises
/// ------
/// `@compileError`: If the context struct does not match the specification.
///
/// Notes
/// -----
/// If the specification struct has no fields, the context struct is allowed to
/// be `void`.
pub fn validateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (specinfo.@"struct".fields.len == 0 and Ctx == void)
        return; // No context needed, nothing to validate

    if (ctxinfo != .@"struct")
        @compileError("Expected struct for context, got " ++ @typeName(Ctx));

    if (ctxinfo.@"struct".fields.len != 0 and specinfo.@"struct".fields.len == 0)
        @compileError("Context struct must be empty");

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = @hasField(@TypeOf(field_spec), "required") and field_spec.required;
        const expected_type = field_spec.type;

        if (@hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);
            const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
                const type_info = @typeInfo(expected_type);
                break :blk type_info == .@"enum";
            } else actual_type == expected_type;

            if (!types_match)
                @compileError("Field '" ++ field.name ++ "' has type " ++ @typeName(actual_type) ++ ", expected " ++ @typeName(expected_type));
        } else if (required) {
            @compileError("Required field '" ++ field.name ++ "' missing from context");
        }
    }

    if (ctxinfo.@"struct".fields.len != 0 and specinfo.@"struct".fields.len == 0)
        @compileError("Context struct must be empty");

    // Check for unexpected fields in context
    inline for (ctxinfo.@"struct".fields) |ctx_field| {
        if (!@hasField(SpecType, ctx_field.name)) {
            @compileError("Unexpected field '" ++ ctx_field.name ++ "' in context");
        }
    }
}

pub fn getFieldOrDefault(ctx: anytype, comptime field_name: []const u8, comptime FieldType: type, default_value: FieldType) FieldType {
    const T = @TypeOf(ctx);
    if (@hasField(T, field_name)) {
        if (@FieldType(T, field_name) != FieldType)
            @compileError("Field '" ++ field_name ++ "' has type " ++ @typeName(@FieldType(T, field_name)) ++ ", expected " ++ @typeName(FieldType));

        return @field(ctx, field_name);
    }

    return default_value;
}

pub fn MixStructs(comptime S1: type, comptime S2: type) type {
    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("MixStructs: both types must be structs");

    for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("Field '" ++ field.name ++ "' already exists in the second struct");
    }

    return @Type(.{
        .@"struct" = .{
            .layout = .auto,
            .fields = info1.@"struct".fields ++ info2.@"struct".fields,
            .decls = &.{},
            .is_tuple = false,
        },
    });
}

pub fn mixStructs(s1: anytype, s2: anytype) MixStructs(@TypeOf(s1), @TypeOf(s2)) {
    const S1 = @TypeOf(s1);
    const S2 = @TypeOf(s2);

    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("mixStructs: both types must be structs");

    var result: MixStructs(S1, S2) = .{};
    inline for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("Field '" ++ field.name ++ "' already exists in the second struct");

        @field(result, field.name) = @field(s1, field.name);
    }
    inline for (info2.@"struct".fields) |field| {
        @field(result, field.name) = @field(s2, field.name);
    }

    return result;
}

pub fn StripStruct(comptime S: type, comptime fields_to_remove: []const []const u8) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    // Calculate how many fields will remain
    comptime var remaining_count: usize = 0;
    inline for (info.@"struct".fields) |field| {
        comptime var should_remove: bool = false;
        inline for (fields_to_remove) |field_to_remove| {
            if (std.mem.eql(u8, field.name, field_to_remove)) {
                should_remove = true;
                break;
            }
        }

        if (!should_remove) {
            remaining_count += 1;
        }
    }

    // Create new fields array
    comptime var new_fields: [remaining_count]std.builtin.Type.StructField = undefined;
    comptime var new_field_index: usize = 0;
    inline for (info.@"struct".fields) |field| {
        var should_remove = false;
        inline for (fields_to_remove) |field_to_remove| {
            if (std.mem.eql(u8, field.name, field_to_remove)) {
                should_remove = true;
                break;
            }
        }

        if (!should_remove) {
            new_fields[new_field_index] = field;
            new_field_index += 1;
        }
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

pub fn stripStruct(s: anytype, comptime fields_to_remove: []const []const u8) StripStruct(@TypeOf(s), fields_to_remove) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    var result: StripStruct(S, fields_to_remove) = .{};
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}
