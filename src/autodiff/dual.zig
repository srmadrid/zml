const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const autodiff = @import("../autodiff.zig");

pub fn isDual(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"struct" => return @hasDecl(T, "is_dual") and T.is_dual,
        else => return false,
    }
}

/// Represents a dual number `x + yϵ`, where `ϵ² = 0`.
pub fn Dual(comptime N: type) type {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.autodiff.Dual: N must be a numeric type, got \n\tT: " ++ @typeName(N) ++ "\n");

    return struct {
        val: N,
        eps: N,

        // Type signature
        pub const is_custom = true;
        pub const is_numeric = true;
        pub const is_dual = true;
        pub const is_complex = meta.isComplex(N);

        pub const Accumulator = Dual(meta.Accumulator(N));
        pub const Real = Dual(meta.Real(N));
        pub const Scalar = N;

        pub const empty: autodiff.Dual(N) = .{
            .val = undefined,
            .eps = undefined,
        };

        // Constants
        pub const zero: autodiff.Dual(N) = .{
            .val = numeric.zero(N),
            .eps = numeric.zero(N),
        };

        // Basic operations
        pub const Abs = autodiff.dual.Abs;
        pub const abs = autodiff.dual.abs;
        pub const Abs1 = autodiff.dual.Abs1;
        pub const abs1 = autodiff.dual.abs1;
        pub const Abs2 = autodiff.dual.Abs2;
        pub const abs2 = autodiff.dual.abs2;
        pub const Neg = autodiff.dual.Neg;
        pub const neg = autodiff.dual.neg;
        pub const Re = autodiff.dual.Re;
        pub const re = autodiff.dual.re;
        pub const Im = autodiff.dual.Im;
        pub const im = autodiff.dual.im;
        pub const Conj = autodiff.dual.Conj;
        pub const conj = autodiff.dual.conj;
        pub const Sign = autodiff.dual.Sign;
        pub const sign = autodiff.dual.sign;

        // Arithmetic operations
        pub const Add = autodiff.dual.Add;
        pub const add = autodiff.dual.add;
        pub const Sub = autodiff.dual.Sub;
        pub const sub = autodiff.dual.sub;
        pub const Mul = autodiff.dual.Mul;
        pub const mul = autodiff.dual.mul;
        pub const Div = autodiff.dual.Div;
        pub const div = autodiff.dual.div;

        // Comparison operations
        // pub const cmp = ops.cmp;
        pub const eq = autodiff.dual.eq;
        pub const ne = autodiff.dual.ne;
        pub const lt = autodiff.dual.lt;
        pub const le = autodiff.dual.le;
        pub const gt = autodiff.dual.gt;
        pub const ge = autodiff.dual.ge;
        pub const Max = autodiff.dual.Max;
        pub const max = autodiff.dual.max;
        pub const Min = autodiff.dual.Min;
        pub const min = autodiff.dual.min;

        // Exponential functions
        pub const Exp = autodiff.dual.Exp;
        pub const exp = autodiff.dual.exp;
        pub const Ln = autodiff.dual.Ln;
        pub const ln = autodiff.dual.ln;

        // Power functions
        pub const Pow = autodiff.dual.Pow;
        pub const pow = autodiff.dual.pow;
        pub const Sqrt = autodiff.dual.Sqrt;
        pub const sqrt = autodiff.dual.sqrt;
        pub const Cbrt = autodiff.dual.Cbrt;
        pub const cbrt = autodiff.dual.cbrt;
        pub const Hypot = autodiff.dual.Hypot;
        pub const hypot = autodiff.dual.hypot;

        // Trigonometric functions
        pub const Sin = autodiff.dual.Sin;
        pub const sin = autodiff.dual.sin;
        pub const Cos = autodiff.dual.Cos;
        pub const cos = autodiff.dual.cos;
        pub const Tan = autodiff.dual.Tan;
        pub const tan = autodiff.dual.tan;
        pub const Asin = autodiff.dual.Asin;
        pub const asin = autodiff.dual.asin;
        pub const Acos = autodiff.dual.Acos;
        pub const acos = autodiff.dual.acos;
        pub const Atan = autodiff.dual.Atan;
        pub const atan = autodiff.dual.atan;
        pub const Atan2 = autodiff.dual.Atan2;
        pub const atan2 = autodiff.dual.atan2;

        // Hyperbolic functions
        pub const Sinh = autodiff.dual.Sinh;
        pub const sinh = autodiff.dual.sinh;
        pub const Cosh = autodiff.dual.Cosh;
        pub const cosh = autodiff.dual.cosh;
        pub const Tanh = autodiff.dual.Tanh;
        pub const tanh = autodiff.dual.tanh;
        pub const Asinh = autodiff.dual.Asinh;
        pub const asinh = autodiff.dual.asinh;
        pub const Acosh = autodiff.dual.Acosh;
        pub const acosh = autodiff.dual.acosh;
        pub const Atanh = autodiff.dual.Atanh;
        pub const atanh = autodiff.dual.atanh;

        pub fn fromFloat(x: anytype) autodiff.Dual(N) {
            return .{
                .val = numeric.cast(N, x),
                .eps = numeric.zero(N),
            };
        }

        pub fn toFloat(self: autodiff.Dual(N), comptime Float: type) Float {
            return numeric.cast(Float, self.val);
        }

        pub fn toCfloat(self: Dual(N), comptime Cfloat: type) Cfloat {
            return meta.scast(Cfloat, self.val);
        }
    };
}

pub fn Abs(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Abs: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Abs(meta.Scalar(X)));
}

pub fn abs(x: anytype) autodiff.dual.Abs(@TypeOf(x)) {
    const absx = numeric.abs(x.val);

    return if (comptime !meta.isComplex(meta.Scalar(@TypeOf(x))))
        .{
            // |x|
            .val = absx,

            // sign(x) * y
            .eps = numeric.mul(numeric.sign(x.val), x.eps),
        }
    else
        .{
            // |x| = √(re(x)² + im(x)²)
            .val = absx,

            // if (x == 0)  0  else  (re(x) * re(y) + im(x) * im(y)) / |x|
            .eps = if (numeric.eq(x.val, numeric.zero(@TypeOf(x.val))))
                numeric.zero(@TypeOf(x.val))
            else
                numeric.div(
                    numeric.add(
                        numeric.mul(numeric.re(x.val), numeric.re(x.eps)),
                        numeric.mul(numeric.im(x.val), numeric.im(x.eps)),
                    ),
                    absx,
                ),
        };
}

pub fn Abs1(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Abs1: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Abs1(meta.Scalar(X)));
}

pub fn abs1(x: anytype) autodiff.dual.Abs1(@TypeOf(x)) {
    return if (comptime !meta.isComplex(meta.Scalar(@TypeOf(x))))
        .{
            // |x|
            .val = numeric.abs1(x.val),

            // sign(x) * y
            .eps = numeric.mul(numeric.sign(x.val), x.eps),
        }
    else
        .{
            // |re(x)| + |im(x)|
            .val = numeric.abs1(x.val),

            // sign(re(x)) * re(y) + sign(im(x)) * im(y)
            .eps = numeric.add(
                numeric.mul(numeric.sign(numeric.re(x.val)), numeric.re(x.eps)),
                numeric.mul(numeric.sign(numeric.im(x.val)), numeric.im(x.eps)),
            ),
        };
}

pub fn Abs2(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Abs2: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Abs2(meta.Scalar(X)));
}

pub fn abs2(x: anytype) autodiff.dual.Abs2(@TypeOf(x)) {
    return if (comptime !meta.isComplex(meta.Scalar(@TypeOf(x))))
        .{
            // |x|²
            .val = numeric.abs2(x.val),

            // 2 * x * y
            .eps = numeric.mul(numeric.two(meta.Scalar(@TypeOf(x))), numeric.mul(x.val, x.eps)),
        }
    else
        .{
            // |x|² = re(x)² + im(x)²
            .val = numeric.abs2(x.val),

            // 2 * (re(x) * re(y) + im(x) * im(y))
            .eps = numeric.mul(
                numeric.two(meta.Scalar(@TypeOf(x))),
                numeric.add(
                    numeric.mul(numeric.re(x.val), numeric.re(x.eps)),
                    numeric.mul(numeric.im(x.val), numeric.im(x.eps)),
                ),
            ),
        };
}

pub fn Neg(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Neg: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Neg(meta.Scalar(X)));
}

pub fn neg(x: anytype) autodiff.dual.Neg(@TypeOf(x)) {
    return .{
        // -x
        .val = numeric.neg(x.val),

        // -y
        .eps = numeric.neg(x.eps),
    };
}

pub fn Re(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Re: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Re(meta.Scalar(X)));
}

pub fn re(x: anytype) autodiff.dual.Re(@TypeOf(x)) {
    return .{
        // re(x)
        .val = numeric.re(x.val),

        // re(y)
        .eps = numeric.re(x.eps),
    };
}

pub fn Im(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Im: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Im(meta.Scalar(X)));
}

pub fn im(x: anytype) autodiff.dual.Im(@TypeOf(x)) {
    return .{
        // im(x)
        .val = numeric.im(x.val),

        // im(y)
        .eps = numeric.im(x.eps),
    };
}

pub fn Conj(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Conj: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Conj(meta.Scalar(X)));
}

pub fn conj(x: anytype) autodiff.dual.Conj(@TypeOf(x)) {
    return .{
        // conj(x)
        .val = numeric.conj(x.val),

        // conj(y)
        .eps = numeric.conj(x.eps),
    };
}

pub fn Sign(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Sign: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Abs(meta.Scalar(X)));
}

pub fn sign(x: anytype) autodiff.dual.Sign(@TypeOf(x)) {
    if (comptime !meta.isComplex(meta.Scalar(@TypeOf(x)))) {
        return .{
            // sign(x)
            .val = numeric.sign(x.val),

            // 0
            .eps = numeric.zero(meta.Scalar(@TypeOf(x))),
        };
    } else {
        if (numeric.eq(x.val, numeric.zero(@TypeOf(x.val)))) {
            return .{
                // 0
                .val = numeric.zero(@TypeOf(x.val)),

                // 0
                .eps = numeric.zero(@TypeOf(x.eps)),
            };
        }

        const absx = numeric.abs(x.val);
        const signx = numeric.div(x.val, absx);

        return .{
            // sign(x) = x / |x|
            .val = signx,

            // (y - sign(x) * (re(x) * re(y) + im(x) * im(y)) / |x|) / |x|
            .eps = numeric.div(
                numeric.sub(
                    x.eps,
                    numeric.mul(
                        signx,
                        numeric.div(
                            numeric.add(
                                numeric.mul(numeric.re(x.val), numeric.re(x.eps)),
                                numeric.mul(numeric.im(x.val), numeric.im(x.eps)),
                            ),
                            absx,
                        ),
                    ),
                ),
                absx,
            ),
        };
    }
}

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Add: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Add(SX, SY));
}

pub fn add(x: anytype, y: anytype) autodiff.dual.Add(@TypeOf(x), @TypeOf(y)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            return .{
                .val = numeric.add(x.val, y.val),
                .eps = numeric.add(x.eps, y.eps),
            };
        } else {
            return .{
                .val = numeric.add(x.val, y),
                .eps = x.eps,
            };
        }
    } else {
        return .{
            .val = numeric.add(x, y.val),
            .eps = y.eps,
        };
    }
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Sub: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Sub(SX, SY));
}

pub fn sub(x: anytype, y: anytype) autodiff.dual.Sub(@TypeOf(x), @TypeOf(y)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            return .{
                .val = numeric.sub(x.val, y.val),
                .eps = numeric.sub(x.eps, y.eps),
            };
        } else {
            return .{
                .val = numeric.sub(x.val, y),
                .eps = x.eps,
            };
        }
    } else {
        return .{
            .val = numeric.sub(x, y.val),
            .eps = numeric.neg(y.eps),
        };
    }
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Mul: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Mul(SX, SY));
}

pub fn mul(x: anytype, y: anytype) autodiff.dual.Mul(@TypeOf(x), @TypeOf(y)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            return .{
                .val = numeric.mul(x.val, y.val),
                .eps = numeric.fma(x.val, y.eps, numeric.mul(x.eps, y.val)),
            };
        } else {
            return .{
                .val = numeric.mul(x.val, y),
                .eps = numeric.mul(x.eps, y),
            };
        }
    } else {
        return .{
            .val = numeric.mul(x, y.val),
            .eps = numeric.mul(x, y.eps),
        };
    }
}

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z) or (!isDual(X) and !isDual(Y) and !isDual(Z)))
        @compileError("zsl.autodiff.dual.Fma: at least one of X, Y or Z must be a dual type, the others must be numeric or dual types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;
    const SZ: type = if (isDual(Z)) meta.Scalar(Z) else Z;

    return Dual(numeric.Fma(SX, SY, SZ));
}

pub fn fma(x: anytype, y: anytype, z: anytype) autodiff.dual.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            if (comptime isDual(@TypeOf(z))) {
                return .{
                    .val = numeric.fma(x.val, y.val, z.val),
                    .eps = numeric.fma(x.eps, y.val, numeric.fma(x.val, y.eps, z.eps)),
                };
            } else {
                return .{
                    .val = numeric.fma(x.val, y.val, z),
                    .eps = numeric.fma(x.eps, y.val, numeric.mul(x.val, y.eps)),
                };
            }
        } else {
            if (comptime isDual(@TypeOf(z))) {
                return .{
                    .val = numeric.fma(x.val, y, z.val),
                    .eps = numeric.fma(x.eps, y, z.eps),
                };
            } else {
                return .{
                    .val = numeric.fma(x.val, y, z),
                    .eps = numeric.mul(x.eps, y),
                };
            }
        }
    } else {
        if (comptime isDual(@TypeOf(y))) {
            if (comptime isDual(@TypeOf(z))) {
                return .{
                    .val = numeric.fma(x, y.val, z.val),
                    .eps = numeric.fma(x, y.eps, z.eps),
                };
            } else {
                return .{
                    .val = numeric.fma(x, y.val, z),
                    .eps = numeric.mul(x, y.eps),
                };
            }
        } else {
            if (comptime isDual(@TypeOf(z))) {
                return .{
                    .val = numeric.fma(x, y, z.val),
                    .eps = numeric.cast(meta.Scalar(autodiff.dual.Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z))), z.eps),
                };
            }
        }
    }
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Div: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Div(SX, SY));
}

pub fn div(x: anytype, y: anytype) autodiff.dual.Div(@TypeOf(x), @TypeOf(y)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            const invy = numeric.div(numeric.one(meta.Scalar(@TypeOf(y))), y.val);

            return .{
                .val = numeric.mul(x.val, invy),
                .eps = numeric.mul(numeric.fma(x.eps, y.val, numeric.neg(numeric.mul(x.val, y.eps))), numeric.mul(invy, invy)),
            };
        } else {
            const invy = numeric.div(numeric.one(@TypeOf(y)), y.val);

            return .{
                .val = numeric.mul(x.val, invy),
                .eps = numeric.mul(x.eps, invy),
            };
        }
    } else {
        const invy = numeric.div(numeric.one(meta.Scalar(@TypeOf(y))), y.val);

        return .{
            .val = numeric.mul(x, invy),
            .eps = numeric.neg(numeric.mul(numeric.mul(x, y.eps), numeric.mul(invy, invy))),
        };
    }
}

pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.eq: at least one of x or y must be a dual, the other must be a numeric or a dual, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.eq(
        if (comptime isDual(X)) x.val else x,
        if (comptime isDual(Y)) y.val else y,
    );
}

pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.ne: at least one of x or y must be a dual, the other must be a numeric or a dual, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.ne(
        if (comptime isDual(X)) x.val else x,
        if (comptime isDual(Y)) y.val else y,
    );
}

pub fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.lt: at least one of x or y must be a dual, the other must be a numeric or a dual, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.lt(
        if (comptime isDual(X)) x.val else x,
        if (comptime isDual(Y)) y.val else y,
    );
}

pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.le: at least one of x or y must be a dual, the other must be a numeric or a dual, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.le(
        if (comptime isDual(X)) x.val else x,
        if (comptime isDual(Y)) y.val else y,
    );
}

pub fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.gt: at least one of x or y must be a dual, the other must be a numeric or a dual, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.gt(
        if (comptime isDual(X)) x.val else x,
        if (comptime isDual(Y)) y.val else y,
    );
}

pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.ge: at least one of x or y must be a dual, the other must be a numeric or a dual, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.ge(
        if (comptime isDual(X)) x.val else x,
        if (comptime isDual(Y)) y.val else y,
    );
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Max: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Max(SX, SY));
}

pub fn max(x: anytype, y: anytype) autodiff.dual.Max(@TypeOf(x), @TypeOf(y)) {
    const R: type = autodiff.dual.Max(@TypeOf(x), @TypeOf(y));

    return if (numeric.gt(x, y)) numeric.cast(R, x) else numeric.cast(R, y);
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Min: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Min(SX, SY));
}

pub fn min(x: anytype, y: anytype) autodiff.dual.Min(@TypeOf(x), @TypeOf(y)) {
    const R: type = autodiff.dual.Min(@TypeOf(x), @TypeOf(y));

    return if (numeric.lt(x, y)) numeric.cast(R, x) else numeric.cast(R, y);
}

pub fn Exp(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Exp: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Exp(meta.Scalar(X)));
}

pub fn exp(x: anytype) autodiff.dual.Exp(@TypeOf(x)) {
    const expx = numeric.exp(x.val);

    return .{
        // eˣ
        .val = expx,

        // y * eˣ
        .eps = numeric.mul(x.eps, expx),
    };
}

pub fn Ln(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Ln: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Ln(meta.Scalar(X)));
}

pub fn ln(x: anytype) autodiff.dual.Ln(@TypeOf(x)) {
    return .{
        // ln(x)
        .val = numeric.ln(x.val),

        // y / x
        .eps = numeric.div(x.eps, x.val),
    };
}

pub fn Pow(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Pow: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Pow(SX, SY));
}

pub fn pow(x: anytype, y: anytype) autodiff.dual.Pow(@TypeOf(x), @TypeOf(y)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            const xpowy = numeric.pow(x.val, y.val);

            return .{
                .val = xpowy,
                .eps = numeric.fma(
                    numeric.mul(
                        y.val,
                        numeric.pow(
                            x.val,
                            numeric.sub(
                                y.val,
                                numeric.one(meta.Scalar(@TypeOf(y))),
                            ),
                        ),
                    ),
                    x.eps,
                    numeric.mul(
                        xpowy,
                        numeric.mul(
                            numeric.ln(x.val),
                            y.eps,
                        ),
                    ),
                ),
            };
        } else {
            const xpowy = numeric.pow(x.val, y);

            return .{
                .val = xpowy,
                .eps = numeric.mul(
                    numeric.mul(
                        y,
                        numeric.pow(
                            x.val,
                            numeric.sub(
                                y,
                                numeric.one(meta.Scalar(@TypeOf(y))),
                            ),
                        ),
                    ),
                    x.eps,
                ),
            };
        }
    } else {
        const xpowy = numeric.pow(x, y.val);

        return .{
            .val = xpowy,
            .eps = numeric.mul(
                xpowy,
                numeric.mul(
                    numeric.ln(x.val),
                    y.eps,
                ),
            ),
        };
    }
}

pub fn Sqrt(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Sqrt: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Sqrt(meta.Scalar(X)));
}

pub fn sqrt(x: anytype) autodiff.dual.Sqrt(@TypeOf(x)) {
    const sqrtx = numeric.sqrt(x.val);

    return .{
        // √x
        .val = sqrtx,

        // y / (2 * √x)
        .eps = numeric.div(x.eps, numeric.mul(numeric.two(@TypeOf(x.val)), sqrtx)),
    };
}

pub fn Cbrt(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Cbrt: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Cbrt(meta.Scalar(X)));
}

pub fn cbrt(x: anytype) autodiff.dual.Cbrt(@TypeOf(x)) {
    const cbrtx = numeric.cbrt(x.val);

    return .{
        // ∛x
        .val = cbrtx,

        // y / (3 * ∛x²)
        .eps = numeric.div(x.eps, numeric.mul(numeric.add(numeric.two(@TypeOf(x.val)), numeric.one(@TypeOf(x.val))), numeric.mul(cbrtx, cbrtx))),
    };
}

pub fn Hypot(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isDual(X) and !isDual(Y)))
        @compileError("zsl.autodiff.dual.Hypot: at least one of X or Y must be a dual type, the other must be a numeric or a dual type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isDual(X)) meta.Scalar(X) else X;
    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;

    return Dual(numeric.Hypot(SX, SY));
}

pub fn hypot(x: anytype, y: anytype) autodiff.dual.Hypot(@TypeOf(x), @TypeOf(y)) {
    if (comptime isDual(@TypeOf(x))) {
        if (comptime isDual(@TypeOf(y))) {
            const hypotxy = numeric.hypot(x.val, y.val);

            return .{
                .val = hypotxy,
                .eps = numeric.div(numeric.fma(x.val, x.eps, numeric.mul(y.val, y.eps)), hypotxy),
            };
        } else {
            const hypotxy = numeric.hypot(x.val, y);

            return .{
                .val = hypotxy,
                .eps = numeric.div(numeric.mul(x.val, x.eps), hypotxy),
            };
        }
    } else {
        const hypotxy = numeric.hypot(x, y.val);

        return .{
            .val = hypotxy,
            .eps = numeric.div(numeric.mul(y.val, y.eps), hypotxy),
        };
    }
}

pub fn Sin(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Sin: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Sin(meta.Scalar(X)));
}

pub fn sin(x: anytype) autodiff.dual.Sin(@TypeOf(x)) {
    return .{
        // sin(x)
        .val = numeric.sin(x.val),

        // y * cos(x)
        .eps = numeric.mul(x.eps, numeric.cos(x.val)),
    };
}

pub fn Cos(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Cos: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Cos(meta.Scalar(X)));
}

pub fn cos(x: anytype) autodiff.dual.Cos(@TypeOf(x)) {
    return .{
        // cos(x)
        .val = numeric.cos(x.val),

        // -y * sin(x)
        .eps = numeric.neg(numeric.mul(x.eps, numeric.sin(x.val))),
    };
}

pub fn Tan(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Tan: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Tan(meta.Scalar(X)));
}

pub fn tan(x: anytype) autodiff.dual.Tan(@TypeOf(x)) {
    const tanx = numeric.tan(x.val);

    return .{
        // tan(x)
        .val = tanx,

        // y * (1 + tan(x)²)
        .eps = numeric.mul(x.eps, numeric.add(numeric.one(@TypeOf(x.val)), numeric.mul(tanx, tanx))),
    };
}

pub fn Asin(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Asin: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Asin(meta.Scalar(X)));
}

pub fn asin(x: anytype) autodiff.dual.Asin(@TypeOf(x)) {
    return .{
        // asin(x)
        .val = numeric.asin(x.val),

        // y / √(1 - x²)
        .eps = numeric.div(x.eps, numeric.sqrt(numeric.sub(numeric.one(@TypeOf(x.val)), numeric.mul(x.val, x.val)))),
    };
}

pub fn Acos(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Acos: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Acos(meta.Scalar(X)));
}

pub fn acos(x: anytype) autodiff.dual.Acos(@TypeOf(x)) {
    return .{
        // acos(x)
        .val = numeric.acos(x.val),

        // -y / √(1 - x²)
        .eps = numeric.neg(numeric.div(x.eps, numeric.sqrt(numeric.sub(numeric.one(@TypeOf(x.val)), numeric.mul(x.val, x.val))))),
    };
}

pub fn Atan(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Atan: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Atan(meta.Scalar(X)));
}

pub fn atan(x: anytype) autodiff.dual.Atan(@TypeOf(x)) {
    return .{
        // atan(x)
        .val = numeric.atan(x.val),

        // y / (1 + x²)
        .eps = numeric.div(x.eps, numeric.add(numeric.one(@TypeOf(x.val)), numeric.mul(x.val, x.val))),
    };
}

pub fn Atan2(comptime Y: type, comptime X: type) type {
    comptime if (!meta.isNumeric(Y) or !meta.isNumeric(X) or (!isDual(Y) and !isDual(X)))
        @compileError("zsl.autodiff.dual.Atan2: at least one of Y or X must be a dual, the other must be a numeric or a dual type, got\n\tY = " ++
            @typeName(Y) ++ "\n\tX = " ++ @typeName(X) ++ "\n");

    const SY: type = if (isDual(Y)) meta.Scalar(Y) else Y;
    const SX: type = if (isDual(X)) meta.Scalar(X) else X;

    return Dual(numeric.Atan2(SY, SX));
}

pub fn atan2(y: anytype, x: anytype) autodiff.dual.Atan2(@TypeOf(y), @TypeOf(x)) {
    if (comptime isDual(@TypeOf(y))) {
        if (comptime isDual(@TypeOf(x))) {
            return .{
                .val = numeric.atan2(y.val, x.val),
                .eps = numeric.div(numeric.fma(x.val, y.eps, numeric.neg(numeric.mul(y.val, x.eps))), numeric.add(numeric.mul(x.val, x.val), numeric.mul(y.val, y.val))),
            };
        } else {
            return .{
                .val = numeric.atan2(y.val, x),
                .eps = numeric.div(numeric.mul(x, y.eps), numeric.add(numeric.mul(x, x), numeric.mul(y.val, y.val))),
            };
        }
    } else {
        return .{
            .val = numeric.atan2(y, x.val),
            .eps = numeric.div(numeric.neg(numeric.mul(y, x.eps)), numeric.add(numeric.mul(x.val, x.val), numeric.mul(y, y))),
        };
    }
}

pub fn Sinh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Sinh: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Sinh(meta.Scalar(X)));
}

pub fn sinh(x: anytype) autodiff.dual.Sinh(@TypeOf(x)) {
    return .{
        // sinh(x)
        .val = numeric.sinh(x.val),

        // y * cosh(x)
        .eps = numeric.mul(x.eps, numeric.cosh(x.val)),
    };
}

pub fn Cosh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Cosh: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Cosh(meta.Scalar(X)));
}

pub fn cosh(x: anytype) autodiff.dual.Cosh(@TypeOf(x)) {
    return .{
        // cosh(x)
        .val = numeric.cosh(x.val),

        // y * sinh(x)
        .eps = numeric.mul(x.eps, numeric.sinh(x.val)),
    };
}

pub fn Tanh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Tanh: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Tanh(meta.Scalar(X)));
}

pub fn tanh(x: anytype) autodiff.dual.Tanh(@TypeOf(x)) {
    const tanhx = numeric.tanh(x.val);

    return .{
        // tanh(x)
        .val = tanhx,

        // y * (1 - tanh(x)²)
        .eps = numeric.mul(x.eps, numeric.sub(numeric.one(@TypeOf(x.val)), numeric.mul(tanhx, tanhx))),
    };
}

pub fn Asinh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Asinh: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Asinh(meta.Scalar(X)));
}

pub fn asinh(x: anytype) autodiff.dual.Asinh(@TypeOf(x)) {
    return .{
        // asinh(x)
        .val = numeric.asinh(x.val),

        // y / √(x² + 1)
        .eps = numeric.div(x.eps, numeric.sqrt(numeric.add(numeric.mul(x.val, x.val), numeric.one(@TypeOf(x.val))))),
    };
}

pub fn Acosh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Acosh: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Acosh(meta.Scalar(X)));
}

pub fn acosh(x: anytype) autodiff.dual.Acosh(@TypeOf(x)) {
    const acoshx = numeric.acosh(x.val);

    return .{
        // acosh(x)
        .val = acoshx,

        // y / sinh(acosh(x))
        .eps = numeric.div(x.eps, numeric.sinh(acoshx)),
    };
}

pub fn Atanh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isDual(X))
        @compileError("zsl.autodiff.dual.Atanh: X must be a dual type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Dual(numeric.Atanh(meta.Scalar(X)));
}

pub fn atanh(x: anytype) autodiff.dual.Atanh(@TypeOf(x)) {
    return .{
        // atanh(x)
        .val = numeric.atanh(x.val),

        // y / (1 - x²)
        .eps = numeric.div(x.eps, numeric.sub(numeric.one(@TypeOf(x.val)), numeric.mul(x.val, x.val))),
    };
}
