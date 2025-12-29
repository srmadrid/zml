//! Namespace for cfloat operations.

const std = @import("std");
const types = @import("types.zig");
const float = @import("float.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const scast = types.scast;

/// Fixed precision complex type, represented as two 16-bit floats.
pub const cf16 = Cfloat(f16);
/// Fixed precision complex type, represented as two 32-bit floats.
pub const cf32 = Cfloat(f32);
/// Fixed precision complex type, represented as two 64-bit floats.
pub const cf64 = Cfloat(f64);
/// Fixed precision complex type, represented as two 80-bit floats.
pub const cf80 = Cfloat(f80);
/// Fixed precision complex type, represented as two 128-bit floats.
pub const cf128 = Cfloat(f128);
/// Compile-time complex type, represented as two compile-time floats.
pub const comptime_cfloat = Cfloat(comptime_float);

pub fn Cfloat(comptime T: type) type {
    if (types.numericType(T) != .float and types.numericType(T) != .dyadic)
        @compileError("Unsupported type for Cfloat: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,

        /// Type signature
        pub const is_cfloat = {};

        /// Scalar type
        pub const Scalar = T;

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

        pub fn conj(self: Cfloat(T)) Cfloat(T) {
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

pub const fma = @import("cfloat/fma.zig").fma;
pub const arg = @import("cfloat/arg.zig").arg;
pub const abs = @import("cfloat/abs.zig").abs;
pub const abs2 = @import("cfloat/abs2.zig").abs2;
pub const neg = @import("cfloat/neg.zig").neg;
pub const exp = @import("cfloat/exp.zig").exp;
pub const log = @import("cfloat/log.zig").log;
pub const pow = @import("cfloat/pow.zig").pow;
pub const sqrt = @import("cfloat/sqrt.zig").sqrt;
pub const sin = @import("cfloat/sin.zig").sin;
pub const cos = @import("cfloat/cos.zig").cos;
pub const tan = @import("cfloat/tan.zig").tan;
pub const asin = @import("cfloat/asin.zig").asin;
pub const acos = @import("cfloat/acos.zig").acos;
pub const atan = @import("cfloat/atan.zig").atan;
pub const sinh = @import("cfloat/sinh.zig").sinh;
pub const cosh = @import("cfloat/cosh.zig").cosh;
pub const tanh = @import("cfloat/tanh.zig").tanh;
pub const asinh = @import("cfloat/asinh.zig").asinh;
pub const acosh = @import("cfloat/acosh.zig").acosh;
pub const atanh = @import("cfloat/atanh.zig").atanh;
