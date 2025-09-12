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

        pub fn add(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            return .{
                .re = x.re + y.re,
                .im = x.im + y.im,
            };
        }

        pub fn addReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re + y,
                .im = x.im,
            };
        }

        pub fn addImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re,
                .im = x.im + y,
            };
        }

        pub fn sub(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            return .{
                .re = x.re - y.re,
                .im = x.im - y.im,
            };
        }

        pub fn subReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re - y,
                .im = x.im,
            };
        }

        pub fn subImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re,
                .im = x.im - y,
            };
        }

        pub fn mul(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            return .{
                .re = @mulAdd(T, x.re, y.re, -x.im * y.im),
                .im = @mulAdd(T, x.re, y.im, x.im * y.re),
            };
        }

        pub fn mulReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re * y,
                .im = x.im * y,
            };
        }

        pub fn mulImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = -x.im * y,
                .im = x.re * y,
            };
        }

        pub fn div(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            if (float.abs(y.im) < float.abs(y.re)) {
                const tmp1 = y.im / y.re;
                const tmp2 = 1 / (y.re + tmp1 * y.im);
                return .{
                    .re = (x.re + x.im * tmp1) * tmp2,
                    .im = (x.im - x.re * tmp1) * tmp2,
                };
            } else {
                const tmp1 = y.re / y.im;
                const tmp2 = 1 / (y.im + tmp1 * y.re);
                return .{
                    .re = (x.re * tmp1 + x.im) * tmp2,
                    .im = (x.im * tmp1 - x.re) * tmp2,
                };
            }
        }

        pub fn divReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re / y,
                .im = x.im / y,
            };
        }

        pub fn divImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.im / y,
                .im = -x.re / y,
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
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.add requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (types.numericType(X)) {
        .bool, .int, .float => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return .{
                        .re = scast(Scalar(C), x) + scast(Scalar(C), y.re),
                        .im = scast(Scalar(C), y.im),
                    };
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => {
                    return .{
                        .re = scast(Scalar(C), x.re) + scast(Scalar(C), y),
                        .im = scast(Scalar(C), x.im),
                    };
                },
                .cfloat => {
                    return .{
                        .re = scast(Scalar(C), x.re) + scast(Scalar(C), y.re),
                        .im = scast(Scalar(C), x.im) + scast(Scalar(C), y.im),
                    };
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn sub(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.sub requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (types.numericType(X)) {
        .bool, .int, .float => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return .{
                        .re = scast(Scalar(C), x) - scast(Scalar(C), y.re),
                        .im = -scast(Scalar(C), y.im),
                    };
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => {
                    return .{
                        .re = scast(Scalar(C), x.re) - scast(Scalar(C), y),
                        .im = scast(Scalar(C), x.im),
                    };
                },
                .cfloat => {
                    return .{
                        .re = scast(Scalar(C), x.re) - scast(Scalar(C), y.re),
                        .im = scast(Scalar(C), x.im) - scast(Scalar(C), y.im),
                    };
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn mul(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.mul requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (types.numericType(X)) {
        .bool, .int, .float => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return .{
                        .re = scast(Scalar(C), x) * scast(Scalar(C), y.re),
                        .im = scast(Scalar(C), x) * scast(Scalar(C), y.im),
                    };
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => {
                    return .{
                        .re = scast(Scalar(C), x.re) * scast(Scalar(C), y),
                        .im = scast(Scalar(C), x.im) * scast(Scalar(C), y),
                    };
                },
                .cfloat => {
                    return scast(C, x).mul(scast(C, y));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn div(
    x: anytype,
    y: anytype,
) Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.div requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (types.numericType(X)) {
        .bool, .int, .float => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return scast(C, y).inverse().mulReal(scast(Scalar(C), x));
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => {
                    return scast(C, x).divReal(scast(Scalar(C), y));
                },
                .cfloat => {
                    return scast(C, x).div(scast(C, y));
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn eq(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(X, Y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.eq requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    switch (types.numericType(X)) {
        .bool, .int, .float => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return y.im == 0 and scast(Scalar(C), x) == scast(Scalar(C), y.re);
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => {
                    return x.im == 0 and scast(Scalar(C), x.re) == scast(Scalar(C), y);
                },
                .cfloat => {
                    return scast(C, x) == scast(C, y);
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub inline fn ne(
    x: anytype,
    y: anytype,
) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if ((types.numericType(X) != .bool and types.numericType(X) != .int and types.numericType(X) != .float and types.numericType(X) != .cfloat) or
        (types.numericType(Y) != .bool and types.numericType(Y) != .int and types.numericType(Y) != .float and types.numericType(Y) != .cfloat) or (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("cfloat.ne requires at least one of x or y to be an int, the other must be a bool, int, float or cfloat, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    return !eq(x, y);
}

pub inline fn abs1(
    z: anytype,
) Scalar(@TypeOf(z)) {
    const Z: type = @TypeOf(z);

    comptime if (types.numericType(Z) != .cfloat)
        @compileError("cfloat.abs1 requires a cfloat argument, got " ++ @typeName(Z));

    return float.abs(z.re) + float.abs(z.im);
}

pub const arg = @import("cfloat/arg.zig").arg; // 11/215 tests fail: 11 for cf128
pub const abs = @import("cfloat/abs.zig").abs; // 0/106 tests fail
pub const abs2 = @import("cfloat/abs2.zig").abs2;
pub const sqrt = @import("cfloat/sqrt.zig").sqrt; // 446/1695 tests fail: 17 for cf32, 68 for cf64, 141 for cf80, 220 for cf128
pub const exp = @import("cfloat/exp.zig").exp; // 28/309 tests fail: 4 for cf32, 4 for cf64, 8 for cf80, 12 for cf128
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
pub const pow = @import("cfloat/pow.zig").pow; // 175/184 tests fail (MANY CATASTYOFIC): 8 for cf32, 19 for cf64, 42 for cf80, 106 for cf128
