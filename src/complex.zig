//! Namespace for complex operations.

const complex = @This();

const meta = @import("meta.zig");
const numeric = @import("numeric.zig");

/// 32-bit complex type.
pub const cf16 = Complex(f16);
/// 64-bit complex type.
pub const cf32 = Complex(f32);
/// 128-bit complex type.
pub const cf64 = Complex(f64);
/// 160-bit complex type.
pub const cf80 = Complex(f80);
/// 256-bit complex type.
pub const cf128 = Complex(f128);
/// Compile-time complex type.
pub const comptime_complex = Complex(comptime_float);

pub fn Complex(comptime N: type) type {
    if (!meta.isNumeric(N) or meta.isIntegral(N))
        @compileError("zsl.Complex: N must be a non-integral numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        re: N,
        im: N,

        // Type signature
        pub const is_numeric = true;
        pub const is_complex = true;

        pub const Accumulator = Complex(meta.Accumulator(N));
        pub const Real = N;
        pub const Scalar = N;

        /// Initializes a complex from any numeric value.
        ///
        /// ## Arguments
        /// * `value` (`anytype`): The value to set the complex to. Must be a
        ///   numeric.
        ///
        /// ## Returns
        /// `Complex(N)`: The new complex.
        pub fn initValue(value: anytype) Complex(N) {
            const V: type = @TypeOf(value);

            comptime if (!meta.isNumeric(V))
                @compileError("zsl.Complex(N).initValue: value must be a numeric, got \n\tvalue: " ++ @typeName(V) ++ "\n");

            switch (comptime meta.numericType(V)) {
                .bool, .int, .float, .dyadic => return .{
                    .re = numeric.cast(N, value),
                    .im = numeric.cast(N, 0),
                },
                .complex => return .{
                    .re = numeric.cast(N, value.re),
                    .im = numeric.cast(N, value.im),
                },
                .custom => return numeric.cast(Complex(N), value),
            }
        }

        pub fn toInt(self: Complex(N), comptime Int: type) Int {
            return numeric.cast(Int, self.re);
        }

        pub fn toFloat(self: Complex(N), comptime Float: type) Float {
            return numeric.cast(Float, self.re);
        }

        // fn parse

        pub fn add(x: Complex(N), y: Complex(N)) Complex(N) {
            return .{
                .re = numeric.add(x.re, y.re),
                .im = numeric.add(x.im, y.im),
            };
        }

        pub fn addReal(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = numeric.add(x.re, y),
                .im = x.im,
            };
        }

        pub fn addImag(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = x.re,
                .im = numeric.add(x.im, y),
            };
        }

        pub fn sub(x: Complex(N), y: Complex(N)) Complex(N) {
            return .{
                .re = numeric.sub(x.re, y.re),
                .im = numeric.sub(x.im, y.im),
            };
        }

        pub fn subReal(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = numeric.sub(x.re, y),
                .im = x.im,
            };
        }

        pub fn subImag(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = x.re,
                .im = numeric.sub(x.im, y),
            };
        }

        pub fn mul(x: Complex(N), y: Complex(N)) Complex(N) {
            return .{
                .re = numeric.fma(x.re, y.re, numeric.mul(numeric.neg(x.im), y.im)),
                .im = numeric.fma(x.re, y.im, numeric.mul(x.im, y.re)),
            };
        }

        pub fn mulReal(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = numeric.mul(x.re, y),
                .im = numeric.mul(x.im, y),
            };
        }

        pub fn mulImag(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = numeric.mul(numeric.neg(x.im), y),
                .im = numeric.mul(x.re, y),
            };
        }

        pub fn div(x: Complex(N), y: Complex(N)) Complex(N) {
            if (numeric.lt(numeric.abs(y.im), numeric.abs(y.re))) {
                const tmp1 = numeric.div(y.im, y.re);
                const tmp2 = numeric.div(1, numeric.fma(tmp1, y.im, y.re));

                return .{
                    .re = numeric.mul(numeric.fma(x.im, tmp1, x.re), tmp2),
                    .im = numeric.mul(numeric.fma(numeric.neg(x.re), tmp1, x.im), tmp2),
                };
            } else {
                const tmp1 = numeric.div(y.re, y.im);
                const tmp2 = numeric.div(1, numeric.fma(tmp1, y.re, y.im));

                return .{
                    .re = numeric.mul(numeric.fma(x.re, tmp1, x.im), tmp2),
                    .im = numeric.mul(numeric.fma(x.im, tmp1, numeric.neg(x.re)), tmp2),
                };
            }
        }

        pub fn divReal(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = numeric.div(x.re, y),
                .im = numeric.div(x.im, y),
            };
        }

        pub fn divImag(x: Complex(N), y: N) Complex(N) {
            return .{
                .re = numeric.div(x.im, y),
                .im = numeric.div(numeric.neg(x.re), y),
            };
        }

        pub fn conj(self: Complex(N)) Complex(N) {
            return .{
                .re = self.re,
                .im = numeric.neg(self.im),
            };
        }

        pub fn neg(self: Complex(N)) Complex(N) {
            return .{
                .re = numeric.neg(self.re),
                .im = numeric.neg(self.im),
            };
        }
    };
}

pub const Coerce = @import("complex/coerce.zig").Coerce;

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex))
        @compileError("zsl.complex.Add: at least one of X or Y must be a complex type, the other must be a bool, an int, a float, a dyadic or a complex type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return complex.Coerce(X, Y);
}

/// Performs addition between two operands of complex, dyadic, float, int or
/// bool types, where at least one operand must be of complex type. The result
/// type is determined by coercing the operand types, and the operation is
/// performed by casting both operands to the result type, then adding them.
///
/// ## Signature
/// ```zig
/// complex.add(x: X, y: Y) complex.Add(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `complex.Add(@TypeOf(x), @TypeOf(y))`: The result of the addition.
pub fn add(x: anytype, y: anytype) complex.Add(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Add(X, Y);

    switch (comptime meta.numericType(X)) {
        .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Y)) {
            .complex => return .addReal(numeric.cast(R, y), numeric.cast(meta.Scalar(R), x)),
            else => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return .addReal(numeric.cast(R, x), numeric.cast(meta.Scalar(R), y)),
            .complex => return .add(numeric.cast(R, x), numeric.cast(R, y)),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex))
        @compileError("zsl.complex.Sub: at least one of X or Y to be a complex type, the other must be a bool, an int, a float or a complex type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return complex.Coerce(X, Y);
}

/// Performs subtraction between two operands of complex, dyadic, float,
/// int or bool types, where at least one operand must be of complex type. The
/// result type is determined by coercing the operand types, and the operation
/// is performed by casting both operands to the result type, then subtracting
/// them.
///
/// ## Signature
/// ```zig
/// complex.sub(x: X, y: Y) complex.Sub(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `complex.Sub(@TypeOf(x), @TypeOf(y))`: The result of the subtraction.
pub fn sub(x: anytype, y: anytype) Sub(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Sub(X, Y);

    switch (comptime meta.numericType(X)) {
        .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Y)) {
            .complex => return .addReal(numeric.cast(R, y).neg(), numeric.cast(meta.Scalar(R), x)),
            else => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return .subReal(numeric.cast(R, x), numeric.cast(meta.Scalar(R), y)),
            .complex => return .sub(numeric.cast(R, x), numeric.cast(R, y)),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex))
        @compileError("zsl.complex.Mul: at least one of X or Y to be a complex type, the other must be a bool, an int, a float or a complex type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return complex.Coerce(X, Y);
}

/// Performs multiplication between two operands of complex, dyadic, float,
/// int or bool types, where at least one operand must be of complex type. The
/// result type is determined by coercing the operand types, and the operation
/// is performed by casting both operands to the result type, then multiplying
/// them.
///
/// ## Signature
/// ```zig
/// complex.mul(x: X, y: Y) complex.Mul(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `complex.Mul(@TypeOf(x), @TypeOf(y))`: The result of the multiplication.
pub fn mul(x: anytype, y: anytype) Mul(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Mul(X, Y);

    switch (comptime meta.numericType(X)) {
        .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Y)) {
            .complex => return .mulReal(numeric.cast(R, y), numeric.cast(meta.Scalar(R), x)),
            else => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return .mulReal(numeric.cast(R, x), numeric.cast(meta.Scalar(R), y)),
            .complex => return .mul(numeric.cast(R, x), numeric.cast(R, y)),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex))
        @compileError("zsl.complex.Div: at least one of X or Y to be a complex type, the other must be a bool, an int, a float or a complex type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return complex.Coerce(X, Y);
}

/// Performs division between two operands of complex, dyadic, float, int or
/// bool types, where at least one operand must be of complex type. The result
/// type is determined by coercing the operand types, and the operation is
/// performed by casting both operands to the result type, then dividing them.
///
/// ## Signature
/// ```zig
/// complex.div(x: X, y: Y) complex.Div(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `complex.Div(@TypeOf(x), @TypeOf(y))`: The result of the division.
pub fn div(x: anytype, y: anytype) Div(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Div(X, Y);

    switch (comptime meta.numericType(X)) {
        .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Y)) {
            .complex => return .div(numeric.cast(R, x), numeric.cast(R, y)),
            else => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return .divReal(numeric.cast(R, x), numeric.cast(meta.Scalar(R), y)),
            .complex => return .div(numeric.cast(R, x), numeric.cast(R, y)),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Compares two operands of complex, dyadic, float, int or bool types, where at
/// least one operand must be of complex type, for equality. The operation is
/// performed by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// complex.eq(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are equal, `false` otherwise.
pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex))
        @compileError("zsl.complex.eq: at least one of x or y to be a complex, the other must be a bool, an int, a float or a complex, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Y)) {
            .complex => return numeric.eq(x, y.re) and numeric.eq(y.im, 0),
            else => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return numeric.eq(x.re, y) and numeric.eq(x.im, 0),
            .complex => return numeric.eq(x.re, y.re) and numeric.eq(x.im, y.im),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

/// Compares two operands of complex, dyadic, float, int or bool types, where at
/// least one operand must be of complex type, for inequality. The operation is
/// performed by casting both operands to the coerced type, then comparing them.
///
/// ## Signature
/// ```zig
/// complex.ne(x: X, y: Y) bool
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left operand.
/// * `y` (`anytype`): The right operand.
///
/// ## Returns
/// `bool`: `true` if the operands are not equal, `false` otherwise.
pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or
        !meta.numericType(X).le(.complex) or !meta.numericType(Y).le(.complex) or
        (meta.numericType(X) != .complex and meta.numericType(Y) != .complex))
        @compileError("zsl.complex.ne: at least one of x or y to be a complex, the other must be a bool, an int, a float or a complex, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    switch (comptime meta.numericType(X)) {
        .bool, .int, .float, .dyadic => switch (comptime meta.numericType(Y)) {
            .complex => return numeric.ne(x, y.re) or numeric.ne(y.im, 0),
            else => unreachable,
        },
        .complex => switch (comptime meta.numericType(Y)) {
            .bool, .int, .float, .dyadic => return numeric.ne(x.re, y) or numeric.ne(x.im, 0),
            .complex => return numeric.ne(x.re, y.re) or numeric.ne(x.im, y.im),
            .custom => unreachable,
        },
        .custom => unreachable,
    }
}

pub const Fma = @import("complex/fma.zig").Fma;
pub const fma = @import("complex/fma.zig").fma;
pub const Arg = @import("complex/arg.zig").Arg;
pub const arg = @import("complex/arg.zig").arg;
pub const Abs = @import("complex/abs.zig").Abs;
pub const abs = @import("complex/abs.zig").abs;
pub const Abs1 = @import("complex/abs1.zig").Abs1;
pub const abs1 = @import("complex/abs1.zig").abs1;
pub const Abs2 = @import("complex/abs2.zig").Abs2;
pub const abs2 = @import("complex/abs2.zig").abs2;
pub const exp = @import("complex/exp.zig").exp;
pub const ln = @import("complex/ln.zig").ln;
pub const Pow = @import("complex/pow.zig").Pow;
pub const pow = @import("complex/pow.zig").pow;
pub const sqrt = @import("complex/sqrt.zig").sqrt;
pub const sin = @import("complex/sin.zig").sin;
pub const cos = @import("complex/cos.zig").cos;
pub const tan = @import("complex/tan.zig").tan;
pub const asin = @import("complex/asin.zig").asin;
pub const acos = @import("complex/acos.zig").acos;
pub const atan = @import("complex/atan.zig").atan;
pub const sinh = @import("complex/sinh.zig").sinh;
pub const cosh = @import("complex/cosh.zig").cosh;
pub const tanh = @import("complex/tanh.zig").tanh;
pub const asinh = @import("complex/asinh.zig").asinh;
pub const acosh = @import("complex/acosh.zig").acosh;
pub const atanh = @import("complex/atanh.zig").atanh;
