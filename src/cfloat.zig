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

pub fn cos(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (z.im == 0) {
        return .{
            .re = float.cos(z.re),
            .im = 0,
        };
    } else {
        return .{
            .re = float.cos(z.re) * float.cosh(z.im),
            .im = float.sin(z.re) * float.sinh(-z.im),
        };
    }
}

pub fn tan(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (@abs(z.im) < 1) {
        const D = 1 / (float.pow(z.re, 2) + float.pow(float.cosh(z.im), 2));

        return .{
            .re = 0.5 * D * float.sin(2 * z.re),
            .im = 0.5 * D * float.sinh(2 * z.im),
        };
    } else {
        const D = 1 / (float.pow(z.re, 2) + float.pow(float.sinh(z.im), 2));
        const F = 1 + float.pow(float.cos(z.re) / float.sinh(z.im), 2);

        return .{
            .re = 0.5 * D * float.sin(2 * z.re),
            .im = 1 / (float.tanh(z.im) * F),
        };
    }
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

pub fn asin(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (@abs(z) <= 1) {
                return .{
                    .re = float.asin(z),
                    .im = 0,
                };
            } else {
                const pi2: @TypeOf(z) = std.math.pi / 2;

                if (z < 0) {
                    return .{
                        .re = -pi2,
                        .im = float.acosh(-z),
                    };
                } else {
                    return .{
                        .re = pi2,
                        .im = -float.acosh(z),
                    };
                }
            }
        },
        .cfloat => {
            if (z.re == 0) {
                return asin(z.re);
            } else {
                const x = @abs(z.re);
                const y = @abs(z.im);
                const r = float.hypot(x + 1, y);
                const s = float.hypot(x - 1, y);
                const A = 0.5 * (r + s);
                const B = x / A;
                const y2 = y * y;

                const A_crossover: Scalar(@TypeOf(z)) = 1.5;
                const B_crossover: Scalar(@TypeOf(z)) = 0.6417;

                var real: Scalar(@TypeOf(z)) = undefined;
                var imag: Scalar(@TypeOf(z)) = undefined;

                if (B <= B_crossover) {
                    real = float.asin(B);
                } else {
                    if (x <= 1) {
                        const D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
                        real = float.atan(x / @sqrt(D));
                    } else {
                        const Apx = A + x;
                        const D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
                        real = float.atan(x / (y * @sqrt(D)));
                    }
                }

                if (A <= A_crossover) {
                    var Am1: Scalar(@TypeOf(z)) = undefined;

                    if (x < 1) {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
                    } else {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
                    }

                    imag = float.log1p(Am1 + @sqrt(Am1 * (A + 1)));
                } else {
                    imag = @log(A + @sqrt(A * A - 1));
                }

                if (z.re >= 0) {
                    if (z.im >= 0) {
                        return .{
                            .re = real,
                            .im = imag,
                        };
                    } else {
                        return .{
                            .re = real,
                            .im = -imag,
                        };
                    }
                } else {
                    if (z.im >= 0) {
                        return .{
                            .re = -real,
                            .im = imag,
                        };
                    } else {
                        return .{
                            .re = -real,
                            .im = -imag,
                        };
                    }
                }
            }
        },
        else => unreachable,
    }
}

pub fn acos(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (@abs(z) <= 1) {
                return .{
                    .re = float.acos(z),
                    .im = 0,
                };
            } else {
                if (z < 0) {
                    return .{
                        .re = std.math.pi,
                        .im = -float.acosh(-z),
                    };
                } else {
                    return .{
                        .re = 0,
                        .im = float.acosh(z),
                    };
                }
            }
        },
        .cfloat => {
            if (z.re == 0) {
                return acos(z.re);
            } else {
                const x = @abs(z.re);
                const y = @abs(z.im);
                const r = float.hypot(x + 1, y);
                const s = float.hypot(x - 1, y);
                const A = 0.5 * (r + s);
                const B = x / A;
                const y2 = y * y;

                const A_crossover: Scalar(@TypeOf(z)) = 1.5;
                const B_crossover: Scalar(@TypeOf(z)) = 0.6417;

                var real: Scalar(@TypeOf(z)) = undefined;
                var imag: Scalar(@TypeOf(z)) = undefined;

                if (B <= B_crossover) {
                    real = float.acos(B);
                } else {
                    if (x <= 1) {
                        const D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
                        real = float.atan(@sqrt(D) / x);
                    } else {
                        const Apx = A + x;
                        const D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
                        real = float.atan((y * @sqrt(D)) / x);
                    }
                }

                if (A <= A_crossover) {
                    var Am1: Scalar(@TypeOf(z)) = undefined;

                    if (x < 1) {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
                    } else {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
                    }

                    imag = float.log1p(Am1 + @sqrt(Am1 * (A + 1)));
                } else {
                    imag = @log(A + @sqrt(A * A - 1));
                }

                if (z.re >= 0) {
                    if (z.im >= 0) {
                        return .{
                            .re = real,
                            .im = -imag,
                        };
                    } else {
                        return .{
                            .re = real,
                            .im = imag,
                        };
                    }
                } else {
                    if (z.im >= 0) {
                        return .{
                            .re = std.math.pi - real,
                            .im = -imag,
                        };
                    } else {
                        return .{
                            .re = std.math.pi - real,
                            .im = imag,
                        };
                    }
                }
            }
        },
        else => unreachable,
    }
}

pub fn atan(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (z.im == 0) {
        return .{
            .re = float.atan(z.re),
            .im = 0,
        };
    } else {
        const r = float.hypot(z.re, z.im);
        const u = 2 * z.im / (1 + r * r);
        var imag: Scalar(@TypeOf(z)) = undefined;

        if (@abs(u) < 0.1) {
            imag = 0.25 * (float.log1p(u) - float.log1p(-u));
        } else {
            const A = float.hypot(z.re, z.im + 1);
            const B = float.hypot(z.re, z.im - 1);
            imag = 0.5 * @log(A / B);
        }

        if (z.re == 0) {
            if (z.im > 1) {
                return .{
                    .re = std.math.pi / 2,
                    .im = imag,
                };
            } else if (z.im < -1) {
                return .{
                    .re = -std.math.pi / 2,
                    .im = imag,
                };
            } else {
                return .{
                    .re = 0,
                    .im = imag,
                };
            }
        } else {
            return .{
                .re = 0.5 * float.atan2(2 * z.re, ((1 + r) * (1 - r))),
                .im = imag,
            };
        }
    }
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

pub fn sinh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return .{
        .re = float.sinh(z.re) * float.cosh(z.im),
        .im = float.cosh(z.re) * float.sinh(z.im),
    };
}

pub fn cosh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return .{
        .re = float.cosh(z.re) * float.cosh(z.im),
        .im = float.sinh(z.re) * float.sinh(z.im),
    };
}

pub fn tanh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (@abs(z.re) < 1) {
        const D = 1 / (float.pow(float.cos(z.re), 2) + float.pow(float.sinh(z.im), 2));

        return .{
            .re = D * float.sinh(z.re) * float.cosh(z.re),
            .im = D * 0.5 * float.sin(2 * z.im),
        };
    } else {
        const D = 1 / (float.pow(float.cos(z.re), 2) + float.pow(float.sinh(z.im), 2));
        const F = 1 + float.pow(float.cos(z.re) / float.sinh(z.im), 2);

        return .{
            .re = 1 / (float.tanh(z.re) * F),
            .im = D * 0.5 * float.sin(2 * z.im),
        };
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

pub fn asinh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return asin(z.mulImag(1)).mulImag(-1);
}

pub fn acosh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (z >= 1) {
                return .{
                    .re = float.acosh(z),
                    .im = 0,
                };
            } else {
                if (z >= -1) {
                    return .{
                        .re = 0,
                        .im = float.acos(z),
                    };
                } else {
                    return .{
                        .re = float.acosh(-z),
                        .im = std.math.pi,
                    };
                }
            }
        },
        .cfloat => {
            const tmp = acos(z);

            return tmp.mulImag(if (z.re > 0) -1 else 1);
        },
        else => unreachable,
    }
}

pub fn atanh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (z > -1 and z < 1) {
                return .{
                    .re = float.atanh(z),
                    .im = 0,
                };
            } else {
                return .{
                    .re = float.atanh(1 / z),
                    .im = if (z < 0) std.math.pi / 2 else -std.math.pi / 2,
                };
            }
        },
        .cfloat => {
            if (z.im == 0) {
                return atanh(z.re);
            } else {
                return atan(z.mulImag(1)).mulImag(-1);
            }
        },
        else => unreachable,
    }
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
