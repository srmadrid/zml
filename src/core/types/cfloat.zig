const std = @import("std");
const types = @import("../types.zig");
const Coerce = types.Coerce;

pub fn Cfloat(comptime T: type) type {
    if (types.numericType(T) != .float) @compileError("Unsupported type for cfloat: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,

        pub fn init(re: T, im: T) Cfloat(T) {
            return Cfloat(T){
                .re = re,
                .im = im,
            };
        }

        pub fn initReal(re: T) Cfloat(T) {
            return Cfloat(T){
                .re = re,
                .im = 0,
            };
        }

        pub fn initImag(im: T) Cfloat(T) {
            return Cfloat(T){
                .re = 0,
                .im = im,
            };
        }

        pub fn initPolar(r: T, theta: T) Cfloat(T) {
            return Cfloat(T){
                .re = r * @cos(theta),
                .im = r * @sin(theta),
            };
        }

        pub fn add(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            return Cfloat(T){
                .re = left.re + right.re,
                .im = left.im + right.im,
            };
        }

        pub fn addReal(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.re + right,
                .im = left.im,
            };
        }

        pub fn addImag(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.re,
                .im = left.im + right,
            };
        }

        pub fn sub(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            return Cfloat(T){
                .re = left.re - right.re,
                .im = left.im - right.im,
            };
        }

        pub fn subReal(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.re - right,
                .im = left.im,
            };
        }

        pub fn subImag(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.re,
                .im = left.im - right,
            };
        }

        pub fn mul(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            return Cfloat(T){
                .re = left.re * right.re - left.im * right.im,
                .im = left.re * right.im + left.im * right.re,
            };
        }

        pub fn mulReal(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.re * right,
                .im = left.im * right,
            };
        }

        pub fn mulImag(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = -left.im * right,
                .im = left.re * right,
            };
        }

        pub fn div(left: Cfloat(T), right: Cfloat(T)) Cfloat(T) {
            if (@abs(right.im) < @abs(right.re)) {
                const tmp1 = right.im / right.re;
                const tmp2 = 1 / (right.re + tmp1 * right.im);
                return Cfloat(T){
                    .re = (left.re + left.im * tmp1) * tmp2,
                    .im = (left.im - left.re * tmp1) * tmp2,
                };
            } else {
                const tmp1 = right.re / right.im;
                const tmp2 = 1 / (right.im + tmp1 * right.re);
                return Cfloat(T){
                    .re = (left.re * tmp1 + left.im) * tmp2,
                    .im = (left.im * tmp1 - left.re) * tmp2,
                };
            }
        }

        pub fn divReal(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.re / right,
                .im = left.im / right,
            };
        }

        pub fn divImag(left: Cfloat(T), right: T) Cfloat(T) {
            return Cfloat(T){
                .re = left.im / right,
                .im = -left.re / right,
            };
        }

        pub fn conjugate(self: Cfloat(T)) Cfloat(T) {
            return Cfloat(T){
                .re = self.re,
                .im = -self.im,
            };
        }

        pub fn negative(self: Cfloat(T)) Cfloat(T) {
            return Cfloat(T){
                .re = -self.re,
                .im = -self.im,
            };
        }

        pub fn inverse(self: Cfloat(T)) Cfloat(T) {
            const s = 1 / std.math.hypot(self.re, self.im);
            return Cfloat(T){
                .re = self.re * s * s,
                .im = -self.im * s * s,
            };
        }
    };
}

pub fn add(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)))
        @compileError("Unsupported types for cfloat.add: " ++ @typeName(@TypeOf(left)) ++ " and " ++ @typeName(@TypeOf(right)));

    comptime if (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right)))
        @compileError("Either left or right must be cfloat");

    std.debug.print("Not implemented yet\n", .{});
}

pub fn arg(z: anytype) @TypeOf(z.re, z.im) {
    if (z.re == 0 and z.im == 0) {
        return 0;
    }

    return std.math.atan2(z.im, z.re);
}

pub fn abs(z: anytype) @TypeOf(z.re, z.im) {
    return std.math.hypot(z.re, z.im);
}

pub fn abs2(z: anytype) @TypeOf(z.re, z.im) {
    return z.re * z.re + z.im * z.im;
}

pub fn logabs(z: anytype) @TypeOf(z.re, z.im) {
    const rabs = @abs(z.re);
    const iabs = @abs(z.im);
    var max: @TypeOf(z.re, z.im) = undefined;
    var u: @TypeOf(z.re, z.im) = undefined;

    if (rabs >= iabs) {
        max = rabs;
        u = iabs / rabs;
    } else {
        max = iabs;
        u = rabs / iabs;
    }

    return std.math.log(@TypeOf(z.re, z.im), std.math.e, max) + std.math.log1p(u * u) / 2;
}

pub fn sqrt(z: anytype) @TypeOf(z) {
    if (z.re == 0 and z.im == 0) {
        return .{
            .re = 0,
            .im = 0,
        };
    }

    const rabs = @abs(z.re);
    const iabs = @abs(z.im);
    var w: @TypeOf(z.re, z.im) = undefined;

    if (rabs >= iabs) {
        const t = iabs / rabs;
        w = @sqrt(rabs) * 0.5 * (1 + @sqrt(1 + t * t));
    } else {
        const t = rabs / iabs;
        w = @sqrt(iabs) * 0.5 * (t + @sqrt(1 + t * t));
    }

    if (z.re >= 0) {
        return .{
            .re = w,
            .im = z.im / (2 * w),
        };
    } else {
        const v = if (z.im >= 0) w else -w;
        return .{
            .re = z.im / (2 * v),
            .im = v,
        };
    }
}

pub fn exp(z: anytype) @TypeOf(z) {
    const rho = @exp(z.re);
    const theta = z.im;

    return .{
        .re = rho * @cos(theta),
        .im = rho * @sin(theta),
    };
}

pub fn pow(left: anytype, right: anytype) Coerce(@TypeOf(left), @TypeOf(right)) {
    comptime if (!types.isFixedPrecision(@TypeOf(left)) or !types.isFixedPrecision(@TypeOf(right)))
        @compileError("Unsupported types for cfloat.add: " ++ @typeName(@TypeOf(left)) ++ " and " ++ @typeName(@TypeOf(right)));

    comptime if (!types.isComplex(@TypeOf(left)) and !types.isComplex(@TypeOf(right)))
        @compileError("Either left or right must be cfloat");

    const lnumeric = types.numericType(@TypeOf(left));
    const rnumeric = types.numericType(@TypeOf(right));

    switch (lnumeric) {
        .int => {
            switch (rnumeric) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return pow(Coerce(@TypeOf(left), @TypeOf(right)){
                        .re = types.cast(types.Numeric(Coerce(@TypeOf(left), @TypeOf(right))), left),
                        .im = 0,
                    }, right);
                },
                else => unreachable,
            }
        },
        .float => {
            switch (rnumeric) {
                .int => unreachable,
                .float => unreachable,
                .cfloat => {
                    return pow(Coerce(@TypeOf(left), @TypeOf(right)){
                        .re = types.cast(types.Numeric(Coerce(@TypeOf(left), @TypeOf(right))), left),
                        .im = 0,
                    }, right);
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (rnumeric) {
                .int => {
                    return pow(left, Coerce(@TypeOf(left), @TypeOf(right)){
                        .re = types.cast(types.Numeric(Coerce(@TypeOf(left), @TypeOf(right))), right),
                        .im = 0,
                    });
                },
                .float => {
                    return pow(left, Coerce(@TypeOf(left), @TypeOf(right)){
                        .re = types.cast(types.Numeric(Coerce(@TypeOf(left), @TypeOf(right))), right),
                        .im = 0,
                    });
                },
                .cfloat => {
                    if (left.re == 0 and left.im == 0) {
                        if (right.re == 0 and right.im == 0) {
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
                    } else if (right.re == 1 and right.im == 0) {
                        return left;
                    } else if (right.re == -1 and right.im == 0) {
                        return right.inverse();
                    } else {
                        const logr = logabs(left);
                        const theta = arg(left);

                        const rho = @exp(logr * right.re - theta * right.im);
                        const beta = logr * right.im + theta * right.re;

                        return .{
                            .re = rho * @cos(beta),
                            .im = rho * @sin(beta),
                        };
                    }
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub const cf16 = Cfloat(f16);
pub const cf32 = Cfloat(f32);
pub const cf64 = Cfloat(f64);
pub const cf80 = Cfloat(f80);
pub const cf128 = Cfloat(f128);
pub const comptime_complex = Cfloat(comptime_float);
