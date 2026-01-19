//! Namespace for cfloat operations.

const std = @import("std");
const types = @import("types.zig");
const ops = @import("ops.zig");
const constants = @import("constants.zig");

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
    if (!types.isNumeric(T) or (types.numericType(T) != .float and types.numericType(T) != .dyadic))
        @compileError("zml.Cfloat: T must be a float or dyadic type, got \n\tT: " ++ @typeName(T) ++ "\n");

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
                .im = constants.zero(T, .{}) catch unreachable,
            };
        }

        pub fn initImag(im: T) Cfloat(T) {
            return .{
                .re = constants.zero(T, .{}) catch unreachable,
                .im = im,
            };
        }

        pub fn initPolar(r: T, theta: T) Cfloat(T) {
            return .{
                .re = ops.mul(
                    r,
                    ops.cos(theta, .{}) catch unreachable,
                    .{},
                ) catch unreachable,
                .im = ops.mul(
                    r,
                    ops.sin(theta, .{}) catch unreachable,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn add(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            return .{
                .re = ops.add(
                    x.re,
                    y.re,
                    .{},
                ) catch unreachable,
                .im = ops.add(
                    x.im,
                    y.im,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn addReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = ops.add(
                    x.re,
                    y,
                    .{},
                ) catch unreachable,
                .im = x.im,
            };
        }

        pub fn addImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re,
                .im = ops.add(
                    x.im,
                    y,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn sub(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            return .{
                .re = ops.sub(
                    x.re,
                    y.re,
                    .{},
                ) catch unreachable,
                .im = ops.sub(
                    x.im,
                    y.im,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn subReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = ops.sub(
                    x.re,
                    y,
                    .{},
                ) catch unreachable,
                .im = x.im,
            };
        }

        pub fn subImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = x.re,
                .im = ops.sub(
                    x.im,
                    y,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn mul(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            return .{
                .re = ops.fma(
                    T,
                    x.re,
                    y.re,
                    ops.mul(
                        ops.neg(x.im, .{}) catch unreachable,
                        y.im,
                        .{},
                    ) catch unreachable,
                    .{},
                ) catch unreachable,
                .im = ops.fma(
                    T,
                    x.re,
                    y.im,
                    ops.mul(
                        x.im,
                        y.re,
                        .{},
                    ) catch unreachable,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn mulReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = ops.mul(
                    x.re,
                    y,
                    .{},
                ) catch unreachable,
                .im = ops.mul(
                    x.im,
                    y,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn mulImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = ops.mul(
                    ops.neg(x.im, .{}) catch unreachable,
                    y,
                    .{},
                ) catch unreachable,
                .im = ops.mul(
                    x.re,
                    y,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn div(x: Cfloat(T), y: Cfloat(T)) Cfloat(T) {
            if (ops.lt(
                ops.abs(y.im, .{}) catch unreachable,
                ops.abs(y.re, .{}) catch unreachable,
                .{},
            ) catch unreachable) {
                const tmp1 = ops.div(
                    y.im,
                    y.re,
                    .{},
                ) catch unreachable;
                const tmp2 = ops.div(
                    1,
                    ops.fma(
                        tmp1,
                        y.im,
                        y.re,
                        .{},
                    ) catch unreachable,
                    .{},
                ) catch unreachable;

                return .{
                    .re = ops.mul(
                        ops.fma(
                            x.im,
                            tmp1,
                            x.re,
                            .{},
                        ) catch unreachable,
                        tmp2,
                        .{},
                    ) catch unreachable,
                    .im = ops.mul(
                        ops.fma(
                            ops.neg(x.re, .{}) catch unreachable,
                            tmp1,
                            x.im,
                            .{},
                        ) catch unreachable,
                        tmp2,
                        .{},
                    ) catch unreachable,
                };
            } else {
                const tmp1 = ops.div(
                    y.re,
                    y.im,
                    .{},
                ) catch unreachable;
                const tmp2 = ops.div(
                    1,
                    ops.fma(
                        tmp1,
                        y.re,
                        y.im,
                        .{},
                    ) catch unreachable,
                    .{},
                ) catch unreachable;

                return .{
                    .re = ops.mul(
                        ops.fma(
                            x.re,
                            tmp1,
                            x.im,
                            .{},
                        ) catch unreachable,
                        tmp2,
                        .{},
                    ) catch unreachable,
                    .im = ops.mul(
                        ops.fma(
                            x.im,
                            tmp1,
                            ops.neg(x.re, .{}) catch unreachable,
                            .{},
                        ) catch unreachable,
                        tmp2,
                        .{},
                    ) catch unreachable,
                };
            }
        }

        pub fn divReal(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = ops.div(
                    x.re,
                    y,
                    .{},
                ) catch unreachable,
                .im = ops.div(
                    x.im,
                    y,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn divImag(x: Cfloat(T), y: T) Cfloat(T) {
            return .{
                .re = ops.div(
                    x.im,
                    y,
                    .{},
                ) catch unreachable,
                .im = ops.div(
                    ops.neg(x.re, .{}) catch unreachable,
                    y,
                    .{},
                ) catch unreachable,
            };
        }

        pub fn conj(self: Cfloat(T)) Cfloat(T) {
            return .{
                .re = self.re,
                .im = ops.neg(self.im, .{}) catch unreachable,
            };
        }

        pub fn neg(self: Cfloat(T)) Cfloat(T) {
            return .{
                .re = ops.neg(self.re, .{}) catch unreachable,
                .im = ops.neg(self.im, .{}) catch unreachable,
            };
        }

        pub fn inverse(self: Cfloat(T)) Cfloat(T) {
            const s = ops.div(
                1,
                ops.hypot(self.re, self.im, .{}) catch unreachable,
                .{},
            ) catch unreachable;
            const s2 = ops.mul(
                s,
                s,
                .{},
            ) catch unreachable;

            return .{
                .re = ops.mul(
                    self.re,
                    s2,
                    .{},
                ) catch unreachable,
                .im = ops.mul(
                    ops.neg(self.im, .{}) catch unreachable,
                    s2,
                    .{},
                ) catch unreachable,
            };
        }
    };
}

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.add: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn add(
    x: anytype,
    y: anytype,
) Add(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Add(X, Y);

    switch (types.numericType(X)) {
        .bool, .int, .float, .dyadic => {
            switch (types.numericType(Y)) {
                .cfloat => return types.scast(R, y).addReal(types.scast(types.Scalar(R), x)),
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .dyadic => return types.scast(R, x).addReal(types.scast(types.Scalar(R), y)),
                .cfloat => return types.scast(R, x).add(types.scast(R, y)),
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.sub: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn sub(
    x: anytype,
    y: anytype,
) Sub(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Sub(X, Y);

    switch (types.numericType(X)) {
        .bool, .int, .float, .dyadic => {
            switch (types.numericType(Y)) {
                .cfloat => return types.scast(R, y).neg().addReal(types.scast(types.Scalar(R), x)),
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .dyadic => return types.scast(R, x).subReal(types.scast(types.Scalar(R), y)),
                .cfloat => return types.scast(R, x).sub(types.scast(R, y)),
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.mul: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn mul(
    x: anytype,
    y: anytype,
) Mul(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Mul(X, Y);

    switch (types.numericType(X)) {
        .bool, .int, .float, .dyadic => {
            switch (types.numericType(Y)) {
                .cfloat => return types.scast(R, y).mulReal(types.scast(types.Scalar(R), x)),
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .dyadic => return types.scast(R, x).mulReal(types.scast(types.Scalar(R), y)),
                .cfloat => return types.scast(R, x).mul(types.scast(R, y)),
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.div: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return types.Coerce(X, Y);
}

pub inline fn div(
    x: anytype,
    y: anytype,
) Div(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Div(X, Y);

    switch (types.numericType(X)) {
        .bool, .int, .float, .dyadic => {
            switch (types.numericType(Y)) {
                .cfloat => return types.scast(R, y).divReal(types.scast(types.Scalar(R), x)),
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .dyadic => return types.scast(R, x).divReal(types.scast(types.Scalar(R), y)),
                .cfloat => return types.scast(R, x).div(types.scast(R, y)),
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

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.eq: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    switch (types.numericType(X)) {
        .bool, .int, .float, .dyadic => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return ops.eq(x, y.re, .{}) catch unreachable and
                        ops.eq(0, y.im, .{}) catch unreachable;
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .dyadic => {
                    return ops.eq(x.re, y, .{}) catch unreachable and
                        ops.eq(x.im, 0, .{}) catch unreachable;
                },
                .cfloat => {
                    return ops.eq(x.re, y.re, .{}) catch unreachable and
                        ops.eq(x.im, y.im, .{}) catch unreachable;
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

    comptime if (!types.isNumeric(X) or !types.isNumeric(Y) or
        !types.numericType(X).le(.cfloat) or !types.numericType(Y).le(.cfloat) or
        (types.numericType(X) != .cfloat and types.numericType(Y) != .cfloat))
        @compileError("zml.cfloat.ne: at least one of x or y to be a cfloat, the other must be a bool, an int, a float or a cfloat, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    switch (types.numericType(X)) {
        .bool, .int, .float, .dyadic => {
            switch (types.numericType(Y)) {
                .cfloat => {
                    return ops.ne(x, y.re, .{}) catch unreachable or
                        ops.ne(0, y.im, .{}) catch unreachable;
                },
                else => unreachable,
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .dyadic => {
                    return ops.ne(x.re, y, .{}) catch unreachable or
                        ops.ne(x.im, 0, .{}) catch unreachable;
                },
                .cfloat => {
                    return ops.ne(x.re, y.re, .{}) catch unreachable or
                        ops.ne(x.im, y.im, .{}) catch unreachable;
                },
                else => unreachable,
            }
        },
        else => unreachable,
    }
}

pub const Fma = @import("cfloat/fma.zig").Fma;
pub const fma = @import("cfloat/fma.zig").fma;
pub const Arg = @import("cfloat/arg.zig").Arg;
pub const arg = @import("cfloat/arg.zig").arg;
pub const Abs = @import("cfloat/abs.zig").Abs;
pub const abs = @import("cfloat/abs.zig").abs;
pub const Abs1 = @import("cfloat/abs1.zig").Abs1;
pub const abs1 = @import("cfloat/abs1.zig").abs1;
pub const Abs2 = @import("cfloat/abs2.zig").Abs2;
pub const abs2 = @import("cfloat/abs2.zig").abs2;
pub const exp = @import("cfloat/exp.zig").exp;
pub const Log = @import("cfloat/log.zig").Log;
pub const log = @import("cfloat/log.zig").log;
pub const Pow = @import("cfloat/pow.zig").Pow;
pub const pow = @import("cfloat/pow.zig").pow;
pub const Sqrt = @import("cfloat/sqrt.zig").Sqrt;
pub const sqrt = @import("cfloat/sqrt.zig").sqrt; // Adapt to dyadics
pub const sin = @import("cfloat/sin.zig").sin;
pub const cos = @import("cfloat/cos.zig").cos;
pub const tan = @import("cfloat/tan.zig").tan; // Adapt to dyadics
pub const asin = @import("cfloat/asin.zig").asin; // Adapt to dyadics
pub const acos = @import("cfloat/acos.zig").acos; // Adapt to dyadics
pub const atan = @import("cfloat/atan.zig").atan; // Adapt to dyadics
pub const sinh = @import("cfloat/sinh.zig").sinh;
pub const cosh = @import("cfloat/cosh.zig").cosh;
pub const tanh = @import("cfloat/tanh.zig").tanh;
pub const asinh = @import("cfloat/asinh.zig").asinh;
pub const acosh = @import("cfloat/acosh.zig").acosh;
pub const atanh = @import("cfloat/atanh.zig").atanh;
