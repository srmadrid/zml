const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
const Coerce = types.Coerce;

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
                .re = r * @cos(theta),
                .im = r * @sin(theta),
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
            if (T == f80 or T == f128) {
                return .{
                    .re = left.re * right.re - left.im * right.im,
                    .im = left.re * right.im + left.im * right.re,
                };
            } else {
                return .{
                    .re = @mulAdd(T, left.re, right.re, -left.im * right.im),
                    .im = @mulAdd(T, left.re, right.im, left.im * right.re),
                };
            }
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
            if (@abs(right.im) < @abs(right.re)) {
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
            const s = 1 / std.math.hypot(self.re, self.im);
            return .{
                .re = self.re * s * s,
                .im = -self.im * s * s,
            };
        }
    };
}

// Check all functions for any necessary casts

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

    if (@TypeOf(z.re, z.im) == f128) {
        return math.atan2_128(z.im, z.re);
    } else {
        return std.math.atan2(z.im, z.re);
    }
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

    if (@TypeOf(z.re, z.im) == f128) {
        return @log(max) + math.log1p128(u * u) / 2;
    } else {
        return @log(max) + std.math.log1p(u * u) / 2;
    }
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

    const l: Coerce(@TypeOf(left), @TypeOf(right)) = types.cast(Coerce(@TypeOf(left), @TypeOf(right)), left, .{});
    const r: Coerce(@TypeOf(left), @TypeOf(right)) = types.cast(Coerce(@TypeOf(left), @TypeOf(right)), right, .{});

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
            .re = rho * @cos(beta),
            .im = rho * @sin(beta),
        };
    }
}

pub fn log(z: anytype) @TypeOf(z) {
    return .{
        .re = logabs(z),
        .im = arg(z),
    };
}

pub fn log10(z: anytype) @TypeOf(z) {
    return log(z).mulReal(1 / @log(10));
}

pub fn logBase(z: anytype, base: anytype) Coerce(@TypeOf(z), @TypeOf(base)) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or !types.isFixedPrecision(@TypeOf(base)))
        @compileError("Unsupported types for cfloat.add: " ++ @typeName(@TypeOf(z)) ++ " and " ++ @typeName(@TypeOf(base)));

    comptime if (!types.isComplex(@TypeOf(z)) and !types.isComplex(@TypeOf(z)))
        @compileError("Either z or base must be cfloat");

    const zz: Coerce(@TypeOf(z), @TypeOf(base)) = types.cast(Coerce(@TypeOf(z), @TypeOf(base)), z, .{});
    const bb: Coerce(@TypeOf(z), @TypeOf(base)) = types.cast(Coerce(@TypeOf(z), @TypeOf(base)), base, .{});

    return log(zz).div(log(bb));
}

pub fn sin(z: anytype) @TypeOf(z) {
    if (z.im == 0) {
        return .{
            .re = @sin(z.re),
            .im = 0,
        };
    } else {
        if (@TypeOf(z.re, z.im) == f128) {
            return .{
                .re = @sin(z.re) * math.cosh128(z.im),
                .im = @cos(z.re) * math.sinh128(z.im),
            };
        } else {
            return .{
                .re = @sin(z.re) * std.math.cosh(z.im),
                .im = @cos(z.re) * std.math.sinh(z.im),
            };
        }
    }
}

pub fn cos(z: anytype) @TypeOf(z) {
    if (z.im == 0) {
        return .{
            .re = @cos(z.re),
            .im = 0,
        };
    } else {
        if (@TypeOf(z.re, z.im) == f128) {
            return .{
                .re = @cos(z.re) * math.cosh128(z.im),
                .im = @sin(z.re) * math.sinh128(-z.im),
            };
        } else {
            return .{
                .re = @cos(z.re) * std.math.cosh(z.im),
                .im = @sin(z.re) * std.math.sinh(-z.im),
            };
        }
    }
}

pub fn tan(z: anytype) @TypeOf(z) {
    if (@abs(z.im) < 1) {
        if (@TypeOf(z.re, z.im) == f128) {
            const D = 1 / (math.pow128(z.re, 2) + math.pow128(math.cosh128(z.im), 2));

            return .{
                .re = 0.5 * D * @sin(2 * z.re),
                .im = 0.5 * D * math.sinh128(2 * z.im),
            };
        } else {
            const D = 1 / (std.math.pow(z.re, 2) + std.math.pow(std.math.cosh(z.im), 2));

            return .{
                .re = 0.5 * D * @sin(2 * z.re),
                .im = 0.5 * D * std.math.sinh(2 * z.im),
            };
        }
    } else {
        if (@TypeOf(z.re, z.im) == f128) {
            const D = 1 / (math.pow128(z.re, 2) + math.pow128(math.sinh128(z.im), 2));
            const F = 1 + math.pow128(@cos(z.re) / math.sinh128(z.im), 2);

            return .{
                .re = 0.5 * D * @sin(2 * z.re),
                .im = 1 / (math.tanh128(z.im) * F),
            };
        } else {
            const D = 1 / (std.math.pow(z.re, 2) + std.math.pow(std.math.sinh(z.im), 2));
            const F = 1 + std.math.pow(@cos(z.re) / std.math.sinh(z.im), 2);

            return .{
                .re = 0.5 * D * @sin(2 * z.re),
                .im = 1 / (std.math.tanh(z.im) * F),
            };
        }
    }
}

pub fn sec(z: anytype) @TypeOf(z) {
    return cos(z).inverse();
}

pub fn csc(z: anytype) @TypeOf(z) {
    return sin(z).inverse();
}

pub fn cot(z: anytype) @TypeOf(z) {
    return tan(z).inverse();
}

pub const cf16 = Cfloat(f16);
pub const cf32 = Cfloat(f32);
pub const cf64 = Cfloat(f64);
pub const cf80 = Cfloat(f80);
pub const cf128 = Cfloat(f128);
pub const comptime_complex = Cfloat(comptime_float);
