const std = @import("std");
const types = @import("../types.zig");
const math = @import("../math.zig");
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
                            .re = rho * @cos(beta),
                            .im = rho * @sin(beta),
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

pub fn log(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return .{
        .re = logabs(z),
        .im = arg(z),
    };
}

pub fn log10(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    return log(z).mulReal(1 / @log(10));
}

pub fn logBase(z: anytype, base: anytype) Coerce(@TypeOf(z), @TypeOf(base)) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or !types.isFixedPrecision(@TypeOf(base)) or (!types.isComplex(@TypeOf(z)) and !types.isComplex(@TypeOf(base))))
        @compileError("At least one of z or base must be cfloat, the other must be an int, float or cfloat");

    const zz: Coerce(@TypeOf(z), @TypeOf(base)) = types.cast(Coerce(@TypeOf(z), @TypeOf(base)), z, .{});
    const bb: Coerce(@TypeOf(z), @TypeOf(base)) = types.cast(Coerce(@TypeOf(z), @TypeOf(base)), base, .{});

    return log(zz).div(log(bb));
}

pub fn sin(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (z.im == 0) {
        return .{
            .re = @sin(z.re),
            .im = 0,
        };
    } else {
        if (Scalar(@TypeOf(z)) == f128) {
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

pub fn cos(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (z.im == 0) {
        return .{
            .re = @cos(z.re),
            .im = 0,
        };
    } else {
        if (Scalar(@TypeOf(z)) == f128) {
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

pub fn tan(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (@abs(z.im) < 1) {
        if (Scalar(@TypeOf(z)) == f128) {
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
        if (Scalar(@TypeOf(z)) == f128) {
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
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.asin128(z),
                        .im = 0,
                    };
                } else {
                    return .{
                        .re = std.math.asin(z),
                        .im = 0,
                    };
                }
            } else {
                const pi2: @TypeOf(z) = std.math.pi / 2;

                if (z < 0) {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = -pi2,
                            .im = math.acosh128(-z),
                        };
                    } else {
                        return .{
                            .re = -pi2,
                            .im = std.math.acosh(-z),
                        };
                    }
                } else {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = pi2,
                            .im = -math.acosh128(z),
                        };
                    } else {
                        return .{
                            .re = pi2,
                            .im = -std.math.acosh(z),
                        };
                    }
                }
            }
        },
        .cfloat => {
            if (z.re == 0) {
                return asin(z.re);
            } else {
                const x = @abs(z.re);
                const y = @abs(z.im);
                const r = std.math.hypot(x + 1, y);
                const s = std.math.hypot(x - 1, y);
                const A = 0.5 * (r + s);
                const B = x / A;
                const y2 = y * y;

                const A_crossover: Scalar(@TypeOf(z)) = 1.5;
                const B_crossover: Scalar(@TypeOf(z)) = 0.6417;

                var real: Scalar(@TypeOf(z)) = undefined;
                var imag: Scalar(@TypeOf(z)) = undefined;

                if (B <= B_crossover) {
                    if (Scalar(@TypeOf(z)) == f128) {
                        real = math.asin128(B);
                    } else {
                        real = std.math.asin(B);
                    }
                } else {
                    if (x <= 1) {
                        const D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
                        if (Scalar(@TypeOf(z)) == f128) {
                            real = math.atan128(x / @sqrt(D));
                        } else {
                            real = std.math.atan(x / @sqrt(D));
                        }
                    } else {
                        const Apx = A + x;
                        const D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
                        if (Scalar(@TypeOf(z)) == f128) {
                            real = math.atan128(x / (y * @sqrt(D)));
                        } else {
                            real = std.math.atan(x / (y * @sqrt(D)));
                        }
                    }
                }

                if (A <= A_crossover) {
                    var Am1: Scalar(@TypeOf(z)) = undefined;

                    if (x < 1) {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
                    } else {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
                    }

                    if (Scalar(@TypeOf(z)) == f128) {
                        imag = math.log1p128(Am1 + @sqrt(Am1 * (A + 1)));
                    } else {
                        imag = std.math.log1p(Am1 + @sqrt(Am1 * (A + 1)));
                    }
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
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.acos128(z),
                        .im = 0,
                    };
                } else {
                    return .{
                        .re = std.math.acos(z),
                        .im = 0,
                    };
                }
            } else {
                if (z < 0) {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = std.math.pi,
                            .im = -math.acosh128(-z),
                        };
                    } else {
                        return .{
                            .re = std.math.pi,
                            .im = -std.math.acosh(-z),
                        };
                    }
                } else {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = 0,
                            .im = math.acosh128(z),
                        };
                    } else {
                        return .{
                            .re = 0,
                            .im = std.math.acosh(z),
                        };
                    }
                }
            }
        },
        .cfloat => {
            if (z.re == 0) {
                return acos(z.re);
            } else {
                const x = @abs(z.re);
                const y = @abs(z.im);
                const r = std.math.hypot(x + 1, y);
                const s = std.math.hypot(x - 1, y);
                const A = 0.5 * (r + s);
                const B = x / A;
                const y2 = y * y;

                const A_crossover: Scalar(@TypeOf(z)) = 1.5;
                const B_crossover: Scalar(@TypeOf(z)) = 0.6417;

                var real: Scalar(@TypeOf(z)) = undefined;
                var imag: Scalar(@TypeOf(z)) = undefined;

                if (B <= B_crossover) {
                    if (Scalar(@TypeOf(z)) == f128) {
                        real = math.acos128(B);
                    } else {
                        real = std.math.acos(B);
                    }
                } else {
                    if (x <= 1) {
                        const D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
                        if (Scalar(@TypeOf(z)) == f128) {
                            real = math.atan128(@sqrt(D) / x);
                        } else {
                            real = std.math.atan(@sqrt(D) / x);
                        }
                    } else {
                        const Apx = A + x;
                        const D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
                        if (Scalar(@TypeOf(z)) == f128) {
                            real = math.atan128((y * @sqrt(D)) / x);
                        } else {
                            real = std.math.atan((y * @sqrt(D)) / x);
                        }
                    }
                }

                if (A <= A_crossover) {
                    var Am1: Scalar(@TypeOf(z)) = undefined;

                    if (x < 1) {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
                    } else {
                        Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
                    }

                    if (Scalar(@TypeOf(z)) == f128) {
                        imag = math.log1p128(Am1 + @sqrt(Am1 * (A + 1)));
                    } else {
                        imag = std.math.log1p(Am1 + @sqrt(Am1 * (A + 1)));
                    }
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
        if (Scalar(@TypeOf(z)) == f128) {
            return .{
                .re = math.atan128(z.re),
                .im = 0,
            };
        } else {
            return .{
                .re = std.math.atan(z.re),
                .im = 0,
            };
        }
    } else {
        const r = std.math.hypot(z.re, z.im);
        const u = 2 * z.im / (1 + r * r);
        var imag: Scalar(@TypeOf(z)) = undefined;

        if (@abs(u) < 0.1) {
            if (Scalar(@TypeOf(z)) == f128) {
                imag = 0.25 * (math.log1p128(u) - math.log1p128(-u));
            } else {
                imag = 0.25 * (std.math.log1p(u) - std.math.log1p(-u));
            }
        } else {
            const A = std.math.hypot(z.re, z.im + 1);
            const B = std.math.hypot(z.re, z.im - 1);
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
            if (Scalar(@TypeOf(z)) == f128) {
                return .{
                    .re = 0.5 * math.atan2_128(2 * z.re, ((1 + r) * (1 - r))),
                    .im = imag,
                };
            } else {
                return .{
                    .re = 0.5 * std.math.atan2(2 * z.re, ((1 + r) * (1 - r))),
                    .im = imag,
                };
            }
        }
    }
}

pub fn asec(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int)
        @compileError("z must be a float or cfloat");

    switch (types.numericType(@TypeOf(z))) {
        .float => {
            if (z <= -1 or z >= 1) {
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.acos128(1 / z),
                        .im = 0,
                    };
                } else {
                    return .{
                        .re = std.math.acos(1 / z),
                        .im = 0,
                    };
                }
            } else {
                if (z >= 0) {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = 0,
                            .im = math.acosh128(1 / z),
                        };
                    } else {
                        return .{
                            .re = 0,
                            .im = std.math.acosh(1 / z),
                        };
                    }
                } else {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = std.math.pi,
                            .im = -math.acosh128(-1 / z),
                        };
                    } else {
                        return .{
                            .re = std.math.pi,
                            .im = -std.math.acosh(-1 / z),
                        };
                    }
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
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.asin128(1 / z),
                        .im = 0,
                    };
                } else {
                    return .{
                        .re = std.math.asin(1 / z),
                        .im = 0,
                    };
                }
            } else {
                if (z >= 0) {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = std.math.pi / 2,
                            .im = -math.acosh128(1 / z),
                        };
                    } else {
                        return .{
                            .re = std.math.pi / 2,
                            .im = -std.math.acosh(1 / z),
                        };
                    }
                } else {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = -std.math.pi / 2,
                            .im = math.acosh128(-1 / z),
                        };
                    } else {
                        return .{
                            .re = -std.math.pi / 2,
                            .im = std.math.acosh(-1 / z),
                        };
                    }
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

    if (Scalar(@TypeOf(z)) == f128) {
        return .{
            .re = math.sinh128(z.re) * @cos(z.im),
            .im = math.cosh128(z.re) * @sin(z.im),
        };
    } else {
        return .{
            .re = std.math.sinh(z.re) * std.math.cosh(z.im),
            .im = std.math.cosh(z.re) * std.math.sinh(z.im),
        };
    }
}

pub fn cosh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (Scalar(@TypeOf(z)) == f128) {
        return .{
            .re = math.cosh128(z.re) * @cos(z.im),
            .im = math.sinh128(z.re) * @sin(z.im),
        };
    } else {
        return .{
            .re = std.math.cosh(z.re) * std.math.cosh(z.im),
            .im = std.math.sinh(z.re) * std.math.sinh(z.im),
        };
    }
}

pub fn tanh(z: anytype) Cfloat(Scalar(@TypeOf(z))) {
    comptime if (!types.isFixedPrecision(@TypeOf(z)) or types.numericType(@TypeOf(z)) == .int or types.numericType(@TypeOf(z)) == .float)
        @compileError("z must be a cfloat");

    if (@abs(z.re) < 1) {
        if (Scalar(@TypeOf(z)) == f128) {
            const D = 1 / (math.pow128(@cos(z.re), 2) + math.pow128(math.sinh128(z.im), 2));

            return .{
                .re = D * math.sinh128(z.re) * math.cosh128(z.im),
                .im = D * 0.5 * @sin(2 * z.im),
            };
        } else {
            const D = 1 / (std.math.pow(@cos(z.re), 2) + std.math.pow(std.math.sinh(z.im), 2));

            return .{
                .re = D * std.math.sinh(z.re) * std.math.cosh(z.re),
                .im = D * 0.5 * @sin(2 * z.im),
            };
        }
    } else {
        if (Scalar(@TypeOf(z)) == f128) {
            const D = 1 / (math.pow128(@cos(z.re), 2) + math.pow128(math.sinh128(z.im), 2));
            const F = 1 + math.pow128(@cos(z.re) / math.sinh128(z.im), 2);

            return .{
                .re = 1 / (math.tanh128(z.re) * F),
                .im = D * 0.5 * @sin(2 * z.im),
            };
        } else {
            const D = 1 / (std.math.pow(@cos(z.re), 2) + std.math.pow(std.math.sinh(z.im), 2));
            const F = 1 + std.math.pow(@cos(z.re) / std.math.sinh(z.im), 2);

            return .{
                .re = 1 / (std.math.tanh(z.re) * F),
                .im = D * 0.5 * @sin(2 * z.im),
            };
        }
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
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.acosh128(z),
                        .im = 0,
                    };
                } else {
                    return .{
                        .re = std.math.acosh(z),
                        .im = 0,
                    };
                }
            } else {
                if (z >= -1) {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = 0,
                            .im = math.acos128(z),
                        };
                    } else {
                        return .{
                            .re = 0,
                            .im = std.math.acos(z),
                        };
                    }
                } else {
                    if (@TypeOf(z) == f128) {
                        return .{
                            .re = math.acosh128(-z),
                            .im = std.math.pi,
                        };
                    } else {
                        return .{
                            .re = std.math.acosh(-z),
                            .im = std.math.pi,
                        };
                    }
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
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.atanh128(z),
                        .im = 0,
                    };
                } else {
                    return .{
                        .re = std.math.atanh(z),
                        .im = 0,
                    };
                }
            } else {
                if (@TypeOf(z) == f128) {
                    return .{
                        .re = math.atanh128(1 / z),
                        .im = if (z < 0) std.math.pi / 2 else -std.math.pi / 2,
                    };
                } else {
                    return .{
                        .re = std.math.atanh(1 / z),
                        .im = if (z < 0) std.math.pi / 2 else -std.math.pi / 2,
                    };
                }
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

test {
    std.testing.refAllDeclsRecursive(@This());
}
