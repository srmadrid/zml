const std = @import("std");
const types = @import("types.zig");
const float = @import("float.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;

pub const cf16 = Cfloat(f16);
pub const cf32 = Cfloat(f32);
pub const cf64 = Cfloat(f64);
pub const cf80 = Cfloat(f80);
pub const cf128 = Cfloat(f128);
pub const comptime_complex = Cfloat(comptime_float);

pub const add = @import("cfloat/add.zig").add;
pub const sub = @import("cfloat/sub.zig").sub;
pub const mul = @import("cfloat/mul.zig").mul;
pub const div = @import("cfloat/div.zig").div;

pub const arg = @import("cfloat/arg.zig").arg; // 11/215 tests fail: 11 for cf128
pub const abs = @import("cfloat/abs.zig").abs; // 0/106 tests fail
pub const abs2 = @import("cfloat/abs2.zig").abs2;
pub const logabs = @import("cfloat/logabs.zig").logabs;
pub const sqrt = @import("cfloat/sqrt.zig").sqrt; // 446/1695 tests fail: 17 for cf32, 68 for cf64, 141 for cf80, 220 for cf128
pub const exp = @import("cfloat/exp.zig").exp; // 28/309 tests fail: 4 for cf32, 4 for cf64, 8 for cf80, 12 for cf128
// pub const pow
pub const log = @import("cfloat/log.zig").log; // 1064/4869 tests fail: 74 for f32, 133 for f64, 217 for f80, 640 for f128
pub const log10 = @import("cfloat/log10.zig").log10; // 2427/4861 tests fail: 278 for f32, 596 for f64, 635 for f80, 918 for f128
pub const sin = @import("cfloat/sin.zig").sin; // 36/260 tests fail: 6 for f32, 4 for f64, 3 for f80, 23 for f128
pub const cos = @import("cfloat/cos.zig").cos; // 37/176 tests fail: 6 for f32, 5 for f64, 2 for f80, 24 for f128
pub const tan = @import("cfloat/tan.zig").tan; // 122/280 tests fail: 12 for f32, 16 for f64, 33 for f80, 61 for f128
pub const asin = @import("cfloat/asin.zig").asin; // 2540/7087 tests fail: 332 for f32, 396 for f64, 713 for f80, 1126 for f128
pub const acos = @import("cfloat/acos.zig").acos; // 2244/7087 tests fail: 312 for f32, 375 for f64, 537 for f80, 1020 for f128
pub const atan = @import("cfloat/atan.zig").atan; // 449/5836 tests fail: 49 for f32, 105 for f64, 34 for f80, 261 for f128
pub const sinh = @import("cfloat/sinh.zig").sinh; // 37/260 tests fail: 7 for f32, 5 for f64, 3 for f80, 22 for f128
pub const cosh = @import("cfloat/cosh.zig").cosh; // 38/176 tests fail: 7 for f32, 5 for f64, 2 for f80, 24 for f128
pub const tanh = @import("cfloat/tanh.zig").tanh; // 129/280 tests fail: 14 for f32, 18 for f64, 35 for f80, 62 for f128
pub const asinh = @import("cfloat/asinh.zig").asinh; // 2543/7087 tests fail: 333 for f32, 396 for f64, 714 for f80, 1100 for f128
pub const acosh = @import("cfloat/acosh.zig").acosh; // 2244/7087 tests fail: 312 for f32, 375 for f64, 537 for f80, 1020 for f128
pub const atanh = @import("cfloat/atanh.zig").atanh; // 449/5836 tests fail: 48 for f32, 105 for f64, 34 for f80, 262 for f128

pub fn Cfloat(comptime T: type) type {
    if (types.numericType(T) != .float) @compileError("Unsupported type for cfloat: " ++ @typeName(T));

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

pub fn pow(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)) or types.numericType(@TypeOf(left)) == .int or types.numericType(@TypeOf(right)) == .int or (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right))))
        @compileError("At least one of left or right must be cfloat, the other must be a float or cfloat");

    switch (types.numericType(@TypeOf(left))) {
        .float => {
            switch (types.numericType(@TypeOf(right))) {
                .float => unreachable,
                .cfloat => {
                    return pow(types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{}), right);
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(@TypeOf(right))) {
                .float => {
                    const l = types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{});
                    const r = types.cast(Scalar(Coerce(@TypeOf(left), @TypeOf(right))), right, .{});

                    if (l.re == 0 and l.im == 0) {
                        if (r == 0) {
                            return .{
                                .re = 1,
                                .im = 0,
                            };
                        } else {
                            return .{
                                .re = 0,
                                .im = 0,
                            };
                        }
                    } else {
                        const logr = logabs(l);
                        const theta = arg(l);
                        const rho = @exp(logr * r);
                        const beta = theta * r;

                        return .{
                            .re = rho * float.cos(beta),
                            .im = rho * float.sin(beta),
                        };
                    }
                },
                .cfloat => {
                    const l = types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{});
                    const r = types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{});

                    if (l.re == 0 and l.im == 0) {
                        if (r.re == 0 and r.im == 0) {
                            return .{
                                .re = 1,
                                .im = 0,
                            };
                        } else {
                            return .{
                                .re = 0,
                                .im = 0,
                            };
                        }
                    } else if (r.re == 1 and r.im == 0) {
                        return l;
                    } else if (r.re == -1 and r.im == 0) {
                        return r.inverse();
                    } else {
                        const logr = logabs(l);
                        const theta = arg(l);

                        const rho = @exp(logr * r.re - theta * r.im);
                        const beta = logr * r.im + theta * r.re;

                        return .{
                            .re = rho * float.cos(beta),
                            .im = rho * float.sin(beta),
                        };
                    }
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub fn logBase(z: anytype, base: anytype) Coerce(@TypeOf(z), @TypeOf(base)) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or !types.isFixedPrecision(@TypeOf(base)) or (!types.isComplex(@TypeOf(z)) and !types.isComplex(@TypeOf(base))))
        @compileError("At least one of z or base must be cfloat, the other must be an int, float or cfloat");

    const zz: Coerce(@TypeOf(z), @TypeOf(base)) = types.cast(Coerce(@TypeOf(z), @TypeOf(base)), z, .{});
    const bb: Coerce(@TypeOf(z), @TypeOf(base)) = types.cast(Coerce(@TypeOf(z), @TypeOf(base)), base, .{});

    return log(zz).div(log(bb));
}

pub fn sec(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return cos(z).inverse();
}

pub fn csc(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return sin(z).inverse();
}

pub fn cot(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return tan(z).inverse();
}

pub fn asec(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (z <= -1 or z >= 1) {
                return .{
                    .re = float.acos(1 / z),
                    .im = 0,
                };
            } else {
                if (z >= 0) {
                    return .{
                        .re = 0,
                        .im = float.acosh(1 / z),
                    };
                } else {
                    return .{
                        .re = std.math.pi,
                        .im = -float.acosh(-1 / z),
                    };
                }
            }
        },
        .cfloat => {
            return acos(z.inverse());
        },
        else => unreachable,
    }
}

pub fn acsc(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (z <= -1 or z >= 1) {
                return .{
                    .re = float.asin(1 / z),
                    .im = 0,
                };
            } else {
                if (z >= 0) {
                    return .{
                        .re = std.math.pi / 2,
                        .im = -float.acosh(1 / z),
                    };
                } else {
                    return .{
                        .re = -std.math.pi / 2,
                        .im = float.acosh(-1 / z),
                    };
                }
            }
        },
        .cfloat => {
            return asin(z.inverse());
        },
        else => unreachable,
    }
}

pub fn acot(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (z.re == 0 and z.im == 0) {
        return .{
            .re = std.math.pi / 2,
            .im = 0,
        };
    } else {
        return atan(z.inverse());
    }
}

pub fn sech(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return cosh(z).inverse();
}

pub fn csch(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return sinh(z).inverse();
}

pub fn coth(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return tanh(z).inverse();
}

pub fn asech(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return acosh(z.inverse());
}

pub fn acsch(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return asinh(z.inverse());
}

pub fn acoth(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return atanh(z.inverse());
}
