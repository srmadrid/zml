const std = @import("std");
const types = @import("types.zig");
const float = @import("float.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const scast = types.scast;

pub const cf16 = Cfloat(f16);
pub const cf32 = Cfloat(f32);
pub const cf64 = Cfloat(f64);
pub const cf80 = Cfloat(f80);
pub const cf128 = Cfloat(f128);
pub const comptime_complex = Cfloat(comptime_float);

pub fn Cfloat(comptime T: type) type {
    if (types.numericType(T) != .float) @compileError("Unsupported type for Cfloat: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,

        pub fn init(re: T, im: T) Cfloat(T) {
            return .{
                .re = re,
                .im = im,
            };
        }

        pub fn initReal(re: T) Cfloat(T) {
            return .{
                .re = re,
                .im = 0,
            };
        }

        pub fn initImag(im: T) Cfloat(T) {
            return .{
                .re = 0,
                .im = im,
            };
        }

        pub fn initPolar(r: T, theta: T) Cfloat(T) {
            return .{
                .re = r * float.cos(theta),
                .im = r * float.sin(theta),
            };
        }

        pub fn add(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            return .{
                .re = left.re + right.re,
                .im = left.im + right.im,
            };
        }

        pub fn addReal(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.re + right,
                .im = left.im,
            };
        }

        pub fn addImag(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.re,
                .im = left.im + right,
            };
        }

        pub fn sub(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            return .{
                .re = left.re - right.re,
                .im = left.im - right.im,
            };
        }

        pub fn subReal(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.re - right,
                .im = left.im,
            };
        }

        pub fn subImag(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.re,
                .im = left.im - right,
            };
        }

        pub fn mul(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            return .{
                .re = @mulAdd(T, left.re, right.re, -left.im * right.im),
                .im = @mulAdd(T, left.re, right.im, left.im * right.re),
            };
        }

        pub fn mulReal(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.re * right,
                .im = left.im * right,
            };
        }

        pub fn mulImag(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = -left.im * right,
                .im = left.re * right,
            };
        }

        pub fn div(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            if (float.abs(right.im) < float.abs(right.re)) {
                const tmp1 = right.im / right.re;
                const tmp2 = 1 / (right.re + tmp1 * right.im);
                return .{
                    .re = (left.re + left.im * tmp1) * tmp2,
                    .im = (left.im - left.re * tmp1) * tmp2,
                };
            } else {
                const tmp1 = right.re / right.im;
                const tmp2 = 1 / (right.im + tmp1 * right.re);
                return .{
                    .re = (left.re * tmp1 + left.im) * tmp2,
                    .im = (left.im * tmp1 - left.re) * tmp2,
                };
            }
        }

        pub fn divReal(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.re / right,
                .im = left.im / right,
            };
        }

        pub fn divImag(left: Cfloat(T), right: T) Cfloat(T) {
            return .{
                .re = left.im / right,
                .im = -left.re / right,
            };
        }

        pub fn conjugate(self: Cfloat(T)) Cfloat(T) {
            return .{
                .re = self.re,
                .im = -self.im,
            };
        }

        pub fn negative(self: Cfloat(T)) Cfloat(T) {
            return .{
                .re = -self.re,
                .im = -self.im,
            };
        }

        pub fn inverse(self: Cfloat(T)) Cfloat(T) {
            const s = 1 / float.hypot(self.re, self.im);
            return .{
                .re = self.re * s * s,
                .im = -self.im * s * s,
            };
        }
    };
}

pub inline fn add(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.add requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    return scast(C, right).addReal(scast(Scalar(C), left));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    return scast(C, left).addReal(scast(Scalar(C), right));
                },
                .cfloat => {
                    return scast(C, left).add(scast(C, right));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn add_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("cfloat.add_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.add_ requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    out.* = scast(O, scast(C, right).addReal(scast(Scalar(C), left)));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    out.* = scast(O, scast(C, left).addReal(scast(Scalar(C), right)));
                },
                .cfloat => {
                    out.* = scast(O, scast(C, left).add(scast(C, right)));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn sub(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.sub requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    return .{
                        .re = scast(Scalar(C), left) - scast(Scalar(C), right.re),
                        .im = -scast(Scalar(C), right.im),
                    };
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    return scast(C, left).subReal(scast(Scalar(C), right));
                },
                .cfloat => {
                    return scast(C, left).sub(scast(C, right));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn sub_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("cfloat.sub_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.sub_ requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    out.* = scast(O, .{
                        .re = scast(Scalar(C), left) - scast(Scalar(C), right.re),
                        .im = -scast(Scalar(C), right.im),
                    });
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    out.* = scast(O, scast(C, left).subReal(scast(Scalar(C), right)));
                },
                .cfloat => {
                    out.* = scast(O, scast(C, left).sub(scast(C, right)));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn mul(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.mul requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    return scast(C, right).mulReal(scast(Scalar(C), left));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    return scast(C, left).mulReal(scast(Scalar(C), right));
                },
                .cfloat => {
                    return scast(C, left).mul(scast(C, right));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn mul_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("cfloat.mul_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.mul_ requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    out.* = scast(O, scast(C, right).mulReal(scast(Scalar(C), left)));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    out.* = scast(O, scast(C, left).mulReal(scast(Scalar(C), right)));
                },
                .cfloat => {
                    out.* = scast(O, scast(C, left).mul(scast(C, right)));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn div(
    left: anytype,
    right: anytype,
) Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.div requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    return scast(C, right).inverse().mulReal(scast(Scalar(C), left));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    return scast(C, left).divReal(scast(Scalar(C), right));
                },
                .cfloat => {
                    return scast(C, left).div(scast(C, right));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn div_(
    out: anytype,
    left: anytype,
    right: anytype,
) void {
    comptime var O: type = @TypeOf(out);
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);
    const C: type = Coerce(L, R);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("cfloat.div_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if ((types.numericType(L) != .bool and types.numericType(L) != .int and types.numericType(L) != .float and types.numericType(L) != .cfloat) or
        (types.numericType(R) != .bool and types.numericType(R) != .int and types.numericType(R) != .float and types.numericType(R) != .cfloat) or (types.numericType(L) != .cfloat and types.numericType(R) != .cfloat))
        @compileError("cfloat.div_ requires at least one of left or right to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(L) ++ " and " ++ @typeName(R));

    comptime if (!types.canCastSafely(C, O))
        @compileError("Cannot cast " ++ @typeName(C) ++ " to " ++
            @typeName(O) ++ " safely");

    switch (types.numericType(L)) {
        .int, .float => {
            switch (types.numericType(R)) {
                .cfloat => {
                    out.* = scast(O, scast(C, right).inverse().mulReal(scast(Scalar(C), left)));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .int, .float => {
                    out.* = scast(O, scast(C, left).divReal(scast(Scalar(C), right)));
                },
                .cfloat => {
                    out.* = scast(O, scast(C, left).div(scast(C, right)));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub const arg = @import("cfloat/arg.zig").arg; // 11/215 tests fail: 11 for cf128
pub const abs = @import("cfloat/abs.zig").abs; // 0/106 tests fail
pub const abs2 = @import("cfloat/abs2.zig").abs2;
pub const sqrt = @import("cfloat/sqrt.zig").sqrt; // 446/1695 tests fail: 17 for cf32, 68 for cf64, 141 for cf80, 220 for cf128
pub const exp = @import("cfloat/exp.zig").exp; // 28/309 tests fail: 4 for cf32, 4 for cf64, 8 for cf80, 12 for cf128
// pub const pow
pub const log = @import("cfloat/log.zig").log; // 1064/4869 tests fail: 74 for cf32, 133 for cf64, 217 for cf80, 640 for cf128
pub const log10 = @import("cfloat/log10.zig").log10; // 2427/4861 tests fail: 278 for cf32, 596 for cf64, 635 for cf80, 918 for cf128
pub const sin = @import("cfloat/sin.zig").sin; // 36/260 tests fail: 6 for cf32, 4 for cf64, 3 for cf80, 23 for cf128
pub const cos = @import("cfloat/cos.zig").cos; // 37/176 tests fail: 6 for cf32, 5 for cf64, 2 for cf80, 24 for cf128
pub const tan = @import("cfloat/tan.zig").tan; // 122/280 tests fail: 12 for cf32, 16 for cf64, 33 for cf80, 61 for cf128
pub const asin = @import("cfloat/asin.zig").asin; // 2540/7087 tests fail: 332 for cf32, 396 for cf64, 713 for cf80, 1126 for cf128
pub const acos = @import("cfloat/acos.zig").acos; // 2244/7087 tests fail: 312 for cf32, 375 for cf64, 537 for cf80, 1020 for cf128
pub const atan = @import("cfloat/atan.zig").atan; // 449/5836 tests fail: 49 for cf32, 105 for cf64, 34 for cf80, 261 for cf128
pub const sinh = @import("cfloat/sinh.zig").sinh; // 37/260 tests fail: 7 for cf32, 5 for cf64, 3 for cf80, 22 for cf128
pub const cosh = @import("cfloat/cosh.zig").cosh; // 38/176 tests fail: 7 for cf32, 5 for cf64, 2 for cf80, 24 for cf128
pub const tanh = @import("cfloat/tanh.zig").tanh; // 129/280 tests fail: 14 for cf32, 18 for cf64, 35 for cf80, 62 for cf128
pub const asinh = @import("cfloat/asinh.zig").asinh; // 2543/7087 tests fail: 333 for cf32, 396 for cf64, 714 for cf80, 1100 for cf128
pub const acosh = @import("cfloat/acosh.zig").acosh; // 2244/7087 tests fail: 312 for cf32, 375 for cf64, 537 for cf80, 1020 for cf128
pub const atanh = @import("cfloat/atanh.zig").atanh; // 449/5836 tests fail: 48 for cf32, 105 for cf64, 34 for cf80, 262 for cf128
pub const pow = @import("cfloat/pow.zig").pow; // 175/184 tests fail (MANY CATASTROFIC): 8 for cf32, 19 for cf64, 42 for cf80, 106 for cf128
