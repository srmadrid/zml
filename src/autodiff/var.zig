const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const int = @import("../int.zig");

const autodiff = @import("../autodiff.zig");
const stats = @import("../stats.zig");

pub fn isVar(comptime T: type) bool {
    switch (comptime @typeInfo(T)) {
        .@"union" => return @hasDecl(T, "is_var") and T.is_var,
        else => return false,
    }
}

pub fn Var(N: type) type {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.autodiff.Var: N must be a numeric type, got \n\tN: " ++ @typeName(N) ++ "\n");

    return union(enum) {
        constant: N,
        tracked: struct {
            tape: *autodiff.Tape(N),
            id: usize,
        },

        // Type signature
        pub const is_numeric = true;
        pub const is_var = true;

        pub const Accumulator = Var(meta.Accumulator(N));
        pub const Real = Var(meta.Real(N));
        pub const Scalar = N;

        // Constants
        pub const zero: autodiff.Var(N) = .{ .constant = numeric.zero(N) };
        pub const one: autodiff.Var(N) = .{ .constant = numeric.one(N) };
        pub const two: autodiff.Var(N) = .{ .constant = numeric.two(N) };

        pub fn init(tape: *autodiff.Tape(N), value: N) Var(N) {
            const id = tape.pushAssumeCapacity(.{
                .op = .@"var",
                .left = int.highest(usize),
                .right = int.highest(usize),
                .val = value,
            });

            return .{ .tracked = .{ .tape = tape, .id = id } };
        }

        // Basic operations
        pub const Abs = autodiff.@"var".Abs;
        pub const abs = autodiff.@"var".abs;
        pub const Abs1 = autodiff.@"var".Abs1;
        pub const abs1 = autodiff.@"var".abs1;
        pub const Abs2 = autodiff.@"var".Abs2;
        pub const abs2 = autodiff.@"var".abs2;
        pub const Neg = autodiff.@"var".Neg;
        pub const neg = autodiff.@"var".neg;
        pub const Re = autodiff.@"var".Re;
        pub const re = autodiff.@"var".re;
        pub const Im = autodiff.@"var".Im;
        pub const im = autodiff.@"var".im;
        pub const Conj = autodiff.@"var".Conj;
        pub const conj = autodiff.@"var".conj;
        pub const Sign = autodiff.@"var".Sign;
        pub const sign = autodiff.@"var".sign;

        // Arithmetic operations
        pub const Add = autodiff.@"var".Add;
        pub const add = autodiff.@"var".add;
        pub const Sub = autodiff.@"var".Sub;
        pub const sub = autodiff.@"var".sub;
        pub const Mul = autodiff.@"var".Mul;
        pub const mul = autodiff.@"var".mul;
        pub const Fma = autodiff.@"var".Fma;
        pub const fma = autodiff.@"var".fma;
        pub const Div = autodiff.@"var".Div;
        pub const div = autodiff.@"var".div;

        // Comparison operations
        // pub const cmp = ops.cmp;
        pub const eq = autodiff.@"var".eq;
        pub const ne = autodiff.@"var".ne;
        pub const lt = autodiff.@"var".lt;
        pub const le = autodiff.@"var".le;
        pub const gt = autodiff.@"var".gt;
        pub const ge = autodiff.@"var".ge;
        pub const Max = autodiff.@"var".Max;
        pub const max = autodiff.@"var".max;
        pub const Min = autodiff.@"var".Min;
        pub const min = autodiff.@"var".min;

        // Exponential functions
        pub const Exp = autodiff.@"var".Exp;
        pub const exp = autodiff.@"var".exp;
        pub const Ln = autodiff.@"var".Ln;
        pub const ln = autodiff.@"var".ln;

        // Power functions
        pub const Pow = autodiff.@"var".Pow;
        pub const pow = autodiff.@"var".pow;
        pub const Sqrt = autodiff.@"var".Sqrt;
        pub const sqrt = autodiff.@"var".sqrt;
        pub const Cbrt = autodiff.@"var".Cbrt;
        pub const cbrt = autodiff.@"var".cbrt;
        pub const Hypot = autodiff.@"var".Hypot;
        pub const hypot = autodiff.@"var".hypot;

        // Trigonometric functions
        pub const Sin = autodiff.@"var".Sin;
        pub const sin = autodiff.@"var".sin;
        pub const Cos = autodiff.@"var".Cos;
        pub const cos = autodiff.@"var".cos;
        pub const Tan = autodiff.@"var".Tan;
        pub const tan = autodiff.@"var".tan;
        pub const Asin = autodiff.@"var".Asin;
        pub const asin = autodiff.@"var".asin;
        pub const Acos = autodiff.@"var".Acos;
        pub const acos = autodiff.@"var".acos;
        pub const Atan = autodiff.@"var".Atan;
        pub const atan = autodiff.@"var".atan;
        pub const Atan2 = autodiff.@"var".Atan2;
        pub const atan2 = autodiff.@"var".atan2;

        // Hyperbolic functions
        pub const Sinh = autodiff.@"var".Sinh;
        pub const sinh = autodiff.@"var".sinh;
        pub const Cosh = autodiff.@"var".Cosh;
        pub const cosh = autodiff.@"var".cosh;
        pub const Tanh = autodiff.@"var".Tanh;
        pub const tanh = autodiff.@"var".tanh;
        pub const Asinh = autodiff.@"var".Asinh;
        pub const asinh = autodiff.@"var".asinh;
        pub const Acosh = autodiff.@"var".Acosh;
        pub const acosh = autodiff.@"var".acosh;
        pub const Atanh = autodiff.@"var".Atanh;
        pub const atanh = autodiff.@"var".atanh;

        pub fn backward(self: Var(N)) void {
            switch (self) {
                .constant => return,
                .tracked => |t| {
                    t.tape.nodes[t.id].grad = numeric.one(N);

                    var i: usize = t.id + 1;
                    while (i > 0) {
                        i -= 1;
                        const node = &t.tape.nodes[i];
                        const g = node.grad;

                        if (numeric.eq(g, 0))
                            continue;

                        switch (node.op) {
                            .@"var" => {},
                            .abs => {
                                // d/dx |x| = x / |x| = x / node.val (x != 0)
                                if (!numeric.eq(node.val, 0)) {
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        numeric.div(numeric.mul(g, t.tape.nodes[node.left].val), node.val),
                                    );
                                }
                            },
                            .abs1 => {
                                // d/dx abs1(x) = sign(x)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.sign(t.tape.nodes[node.left].val)),
                                );
                            },
                            .abs2 => {
                                // d/dx |x|² = 2x
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.mul(numeric.two(N), t.tape.nodes[node.left].val)),
                                );
                            },
                            .neg => {
                                // d/dx (-x) = -1
                                numeric.subInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    g,
                                );
                            },
                            .re => {
                                // d/dx re(x) = re(g)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.re(g),
                                );
                            },
                            .im => {
                                // d/dx im(x) = im(g)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.im(g),
                                );
                            },
                            .conj => {
                                // d/dx conj(x) = conj(g)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.conj(g),
                                );
                            },
                            .sign => {
                                // d/dx sign(x) = 0 (almost everywhere)
                            },
                            .add => {
                                // d/dx (x + y) = 1
                                if (node.left != int.highest(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    );

                                // d/dy (x + y) = 1
                                if (node.right != int.highest(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                            },
                            .sub => {
                                // d/dx (x - y) = 1
                                if (node.left != int.highest(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    );

                                // d/dy (x - y) = -1
                                if (node.right != int.highest(usize))
                                    numeric.subInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                            },
                            .mul => {
                                // d/dx (x * y) = y
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, t.tape.nodes[node.right].val),
                                );

                                // d/dy (x * y) = x
                                numeric.addInto(
                                    &t.tape.nodes[node.right].grad,
                                    t.tape.nodes[node.right].grad,
                                    numeric.mul(g, t.tape.nodes[node.left].val),
                                );
                            },
                            .div => {
                                // d/dx (x / y) = 1 / y
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        t.tape.nodes[node.right].val,
                                    ),
                                );

                                // d/dy (x / y) = -x / y²
                                numeric.subInto(
                                    &t.tape.nodes[node.right].grad,
                                    t.tape.nodes[node.right].grad,
                                    numeric.div(
                                        numeric.mul(g, t.tape.nodes[node.left].val),
                                        numeric.mul(t.tape.nodes[node.right].val, t.tape.nodes[node.right].val),
                                    ),
                                );
                            },
                            .max => {
                                // d/dx max(x, y) = if x > y  1  else  0
                                if (numeric.gt(t.tape.nodes[node.left].val, t.tape.nodes[node.right].val))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    )
                                else // d/dy max(x, y) = if y >= x  1  else  0
                                    numeric.addInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                            },
                            .min => {
                                // d/dx min(x, y) = if x < y  1  else  0
                                if (numeric.lt(t.tape.nodes[node.left].val, t.tape.nodes[node.right].val))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    )
                                else // d/dy min(x, y) = if y <= x  1  else  0
                                    numeric.addInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                            },
                            .exp => {
                                // d/dx exp(x) = exp(x) = node.val
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, node.val),
                                );
                            },
                            .ln => {
                                // d/dx ln(x) = 1 / x
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(g, t.tape.nodes[node.left].val),
                                );
                            },
                            .pow => {
                                // d/dx (x^y) = y * x^(y - 1) = (y * node.val) / x
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(
                                        g,
                                        numeric.div(
                                            numeric.mul(t.tape.nodes[node.right].val, node.val),
                                            t.tape.nodes[node.left].val,
                                        ),
                                    ),
                                );

                                // d/dy (x^y) = x^y * ln(x) = node.val * ln(x)
                                numeric.addInto(
                                    &t.tape.nodes[node.right].grad,
                                    t.tape.nodes[node.right].grad,
                                    numeric.mul(
                                        g,
                                        numeric.mul(
                                            node.val,
                                            numeric.ln(t.tape.nodes[node.left].val),
                                        ),
                                    ),
                                );
                            },
                            .sqrt => {
                                // d/dx sqrt(x) = 1 / (2 * sqrt(x)) = 1 / (2 * node.val)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(g, numeric.mul(numeric.two(N), node.val)),
                                );
                            },
                            .cbrt => {
                                // d/dx cbrt(x) = 1 / (3 * cbrt(x)²) = 1 / (3 * node.val²)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.mul(numeric.add(numeric.one(N), numeric.two(N)), numeric.mul(node.val, node.val)),
                                    ),
                                );
                            },
                            .hypot => {
                                // d/dx hypot(x, y) = x / hypot(x, y) = x / node.val
                                if (node.left != int.highest(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        numeric.div(
                                            numeric.mul(g, t.tape.nodes[node.left].val),
                                            node.val,
                                        ),
                                    );

                                // d/dy hypot(x, y) = y / hypot(x, y) = y / node.val
                                if (node.right != int.highest(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        numeric.div(
                                            numeric.mul(g, t.tape.nodes[node.right].val),
                                            node.val,
                                        ),
                                    );
                            },
                            .sin => {
                                // d/dx sin(x) = cos(x)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.cos(t.tape.nodes[node.left].val)),
                                );
                            },
                            .cos => {
                                // d/dx cos(x) = -sin(x)
                                numeric.subInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.sin(t.tape.nodes[node.left].val)),
                                );
                            },
                            .tan => {
                                // d/dx tan(x) = sec²(x) = 1 + tan²(x) = 1 + node.val²
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.add(numeric.one(N), numeric.mul(node.val, node.val))),
                                );
                            },
                            .asin => {
                                // d/dx asin(x) = 1 / sqrt(1 - x²)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.sqrt(
                                            numeric.sub(
                                                numeric.one(N),
                                                numeric.mul(
                                                    t.tape.nodes[node.left].val,
                                                    t.tape.nodes[node.left].val,
                                                ),
                                            ),
                                        ),
                                    ),
                                );
                            },
                            .acos => {
                                // d/dx acos(x) = -1 / sqrt(1 - x²)
                                numeric.subInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.sqrt(
                                            numeric.sub(
                                                numeric.one(N),
                                                numeric.mul(
                                                    t.tape.nodes[node.left].val,
                                                    t.tape.nodes[node.left].val,
                                                ),
                                            ),
                                        ),
                                    ),
                                );
                            },
                            .atan => {
                                // d/dx atan(x) = 1 / (1 + x²)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.add(
                                            numeric.one(N),
                                            numeric.mul(
                                                t.tape.nodes[node.left].val,
                                                t.tape.nodes[node.left].val,
                                            ),
                                        ),
                                    ),
                                );
                            },
                            .atan2 => {
                                const denom = numeric.add(
                                    numeric.mul(
                                        t.tape.nodes[node.left].val,
                                        t.tape.nodes[node.left].val,
                                    ),
                                    numeric.mul(
                                        t.tape.nodes[node.right].val,
                                        t.tape.nodes[node.right].val,
                                    ),
                                );

                                // d/dx atan2(x, y) = y / (x² + y²)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        numeric.mul(
                                            g,
                                            t.tape.nodes[node.right].val,
                                        ),
                                        denom,
                                    ),
                                );

                                // d/dy atan2(x, y) = -x / (x² + y²)
                                numeric.subInto(
                                    &t.tape.nodes[node.right].grad,
                                    t.tape.nodes[node.right].grad,
                                    numeric.div(
                                        numeric.mul(
                                            g,
                                            t.tape.nodes[node.left].val,
                                        ),
                                        denom,
                                    ),
                                );
                            },
                            .sinh => {
                                // d/dx sinh(x) = cosh(x)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(
                                        g,
                                        numeric.cosh(t.tape.nodes[node.left].val),
                                    ),
                                );
                            },
                            .cosh => {
                                // d/dx cosh(x) = sinh(x)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(
                                        g,
                                        numeric.sinh(t.tape.nodes[node.left].val),
                                    ),
                                );
                            },
                            .tanh => {
                                // d/dx tanh(x) = sech²(x) = 1 - tanh²(x) = 1 - node.val²
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(
                                        g,
                                        numeric.sub(
                                            numeric.one(N),
                                            numeric.mul(node.val, node.val),
                                        ),
                                    ),
                                );
                            },
                            .asinh => {
                                // d/dx asinh(x) = 1 / sqrt(x² + 1)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.sqrt(
                                            numeric.add(
                                                numeric.mul(
                                                    t.tape.nodes[node.left].val,
                                                    t.tape.nodes[node.left].val,
                                                ),
                                                numeric.one(N),
                                            ),
                                        ),
                                    ),
                                );
                            },
                            .acosh => {
                                // d/dx acosh(x) = 1 / sqrt(x² - 1)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.sqrt(
                                            numeric.sub(
                                                numeric.mul(
                                                    t.tape.nodes[node.left].val,
                                                    t.tape.nodes[node.left].val,
                                                ),
                                                numeric.one(N),
                                            ),
                                        ),
                                    ),
                                );
                            },
                            .atanh => {
                                // d/dx atanh(x) = 1 / (1 - x²)
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        numeric.sub(
                                            numeric.one(N),
                                            numeric.mul(
                                                t.tape.nodes[node.left].val,
                                                t.tape.nodes[node.left].val,
                                            ),
                                        ),
                                    ),
                                );
                            },
                        }
                    }
                },
            }
        }

        pub fn val(self: Var(N)) N {
            switch (self) {
                .constant => |c| return c,
                .tracked => |t| return t.tape.nodes[t.id].val,
            }
        }

        pub fn grad(self: Var(N)) N {
            switch (self) {
                .constant => return numeric.zero(N),
                .tracked => |t| return t.tape.nodes[t.id].grad,
            }
        }

        pub fn fromFloat(x: anytype) autodiff.Var(N) {
            return .{ .constant = numeric.cast(N, x) };
        }

        pub fn toFloat(self: autodiff.Var(N), comptime Float: type) Float {
            return numeric.cast(Float, self.val());
        }

        pub fn toComplex(self: autodiff.Var(N), comptime Complex: type) Complex {
            return numeric.cast(Complex, self.val());
        }
    };
}

pub fn Abs(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Abs: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn abs(x: anytype) autodiff.@"var".Abs(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Abs(X) = .{ .constant = undefined };

            numeric.absInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .abs,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.absInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Abs1(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Abs1: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn abs1(x: anytype) autodiff.@"var".Abs1(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Abs1(X) = .{ .constant = undefined };

            numeric.abs1Into(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .abs1,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.abs1Into(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Abs2(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Abs2: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn abs2(x: anytype) autodiff.@"var".Abs2(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Abs2(X) = .{ .constant = undefined };

            numeric.abs2Into(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .abs2,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.abs2Into(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Neg(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Neg: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn neg(x: anytype) autodiff.@"var".Neg(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Neg(X) = .{ .constant = undefined };

            numeric.negInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .neg,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.negInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Re(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Re: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn re(x: anytype) autodiff.@"var".Re(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Re(X) = .{ .constant = undefined };

            numeric.reInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .re,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.reInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Im(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Im: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn im(x: anytype) autodiff.@"var".Im(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Im(X) = .{ .constant = undefined };

            numeric.imInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .im,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.imInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Conj(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Conj: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn conj(x: anytype) autodiff.@"var".Conj(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Conj(X) = .{ .constant = undefined };

            numeric.conjInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .conj,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.conjInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Sign(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Sign: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn sign(x: anytype) autodiff.@"var".Sign(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Sign(X) = .{ .constant = undefined };

            numeric.signInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .sign,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.signInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Add: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Add: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn add(x: anytype, y: anytype) autodiff.@"var".Add(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Add(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.addInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .add,
        .left = idOf(x),
        .right = idOf(y),
        .val = undefined,
    };

    numeric.addInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Sub: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Sub: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn sub(x: anytype, y: anytype) autodiff.@"var".Sub(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Sub(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.subInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .sub,
        .left = idOf(x),
        .right = idOf(y),
        .val = undefined,
    };

    numeric.subInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Mul: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Mul: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn mul(x: anytype, y: anytype) autodiff.@"var".Mul(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Mul(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.mulInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .mul,
        .left = ensureTracked(tape, x),
        .right = ensureTracked(tape, y),
        .val = undefined,
    };

    numeric.mulInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z) or (!isVar(X) and !isVar(Y) and !isVar(Z)))
        @compileError("zsl.autodiff.@\"var\".Fma: at least one of X, Y or Z must be a var type, the others must be numeric or var types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    switch (comptime isVar(X)) {
        true => switch (comptime isVar(Y)) {
            true => switch (comptime isVar(Z)) {
                true => if (comptime X != Y or X != Z)
                    @compileError("zsl.autodiff.@\"var\".Fma: if X, Y and Z are all var types, they must be equal, got\n\tX = " ++
                        @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
                false => if (comptime X != Y)
                    @compileError("zsl.autodiff.@\"var\".Fma: if X and Y are all both types, they must be equal, got\n\tX = " ++
                        @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
            },
            false => switch (comptime isVar(Z)) {
                true => if (comptime X != Z)
                    @compileError("zsl.autodiff.@\"var\".Fma: if X and Z are all both types, they must be equal, got\n\tX = " ++
                        @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
                false => {},
            },
        },
        false => switch (comptime isVar(Y)) {
            true => switch (comptime isVar(Z)) {
                true => if (comptime Y != Z)
                    @compileError("zsl.autodiff.@\"var\".Fma: if Y and Z are all both types, they must be equal, got\n\tX = " ++
                        @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
                false => {},
            },
            false => switch (comptime isVar(Z)) {
                true => {},
                false => unreachable,
            },
        },
    }

    return if (comptime isVar(X)) X else if (comptime isVar(Y)) Y else Z;
}

pub fn fma(x: anytype, y: anytype, z: anytype) autodiff.@"var".Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const Z = @TypeOf(z);
    const R = autodiff.@"var".Fma(X, Y, Z);

    const x_val = valOf(x);
    const y_val = valOf(y);
    const z_val = valOf(z);

    if (!isTracked(x) and !isTracked(y) and !isTracked(z)) {
        var result: R = .{ .constant = undefined };
        numeric.fmaInto(&result.constant, x_val, y_val, z_val);
        return result;
    }

    const tape = getTape3(x, y, z);

    const idx = ensureTracked(tape, x);
    const idy = ensureTracked(tape, y);
    const idz = ensureTracked(tape, z);

    var mul_val: meta.Scalar(R) = undefined;
    numeric.mulInto(&mul_val, x_val, y_val);

    const id_mul = tape.pushAssumeCapacity(.{
        .op = .mul,
        .left = idx,
        .right = idy,
        .val = mul_val,
    });

    var fma_val: meta.Scalar(R) = undefined;
    numeric.fmaInto(&fma_val, x_val, y_val, z_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(.{
                .op = .add,
                .left = id_mul,
                .right = idz,
                .val = fma_val,
            }),
        },
    };
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Div: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Div: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn div(x: anytype, y: anytype) autodiff.@"var".Div(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Div(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.divInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .div,
        .left = ensureTracked(tape, x),
        .right = ensureTracked(tape, y),
        .val = undefined,
    };

    numeric.divInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn eq(x: anytype, y: @TypeOf(x)) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".eq: x and y must be vars, got\n\tx, y: " ++
            @typeName(X) ++ "\n");

    return numeric.eq(x.val(), y.val());
}

pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".ne: x and y must be vars, got\n\tx, y: " ++
            @typeName(X) ++ "\n");

    return numeric.ne(x.val(), y.val());
}

pub fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".lt: x and y must be vars, got\n\tx, y: " ++
            @typeName(X) ++ "\n");

    return numeric.lt(x.val(), y.val());
}

pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".le: x and y must be vars, got\n\tx, y: " ++
            @typeName(X) ++ "\n");

    return numeric.le(x.val(), y.val());
}

pub fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".gt: x and y must be vars, got\n\tx, y: " ++
            @typeName(X) ++ "\n");

    return numeric.gt(x.val(), y.val());
}

pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);

    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".ge: x and y must be vars, got\n\tx, y: " ++
            @typeName(X) ++ "\n");

    return numeric.ge(x.val(), y.val());
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Max: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Max: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn max(x: anytype, y: anytype) autodiff.@"var".Max(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Max(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.maxInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .max,
        .left = ensureTracked(tape, x),
        .right = ensureTracked(tape, y),
        .val = undefined,
    };

    numeric.maxInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Min: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Min: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn min(x: anytype, y: anytype) autodiff.@"var".Min(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Min(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.minInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .min,
        .left = ensureTracked(tape, x),
        .right = ensureTracked(tape, y),
        .val = undefined,
    };

    numeric.minInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Exp(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Exp: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn exp(x: anytype) autodiff.@"var".Exp(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Exp(X) = .{ .constant = undefined };

            numeric.expInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .exp,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.expInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Ln(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Ln: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn ln(x: anytype) autodiff.@"var".Ln(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Ln(X) = .{ .constant = undefined };

            numeric.lnInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .ln,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.lnInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Pow(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Pow: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Pow: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn pow(x: anytype, y: anytype) autodiff.@"var".Pow(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Pow(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.powInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .pow,
        .left = ensureTracked(tape, x),
        .right = ensureTracked(tape, y),
        .val = undefined,
    };

    numeric.powInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Sqrt(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Sqrt: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn sqrt(x: anytype) autodiff.@"var".Sqrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Sqrt(X) = .{ .constant = undefined };

            numeric.sqrtInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .sqrt,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.sqrtInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Cbrt(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Cbrt: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn cbrt(x: anytype) autodiff.@"var".Cbrt(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Cbrt(X) = .{ .constant = undefined };

            numeric.cbrtInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .cbrt,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.cbrtInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Hypot(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Hypot: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Hypot: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn hypot(x: anytype, y: anytype) autodiff.@"var".Hypot(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Hypot(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.hypotInto(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .hypot,
        .left = idOf(x),
        .right = idOf(y),
        .val = undefined,
    };

    numeric.hypotInto(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Sin(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Sin: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn sin(x: anytype) autodiff.@"var".Sin(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Sin(X) = .{ .constant = undefined };

            numeric.sinInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .sin,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.sinInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Cos(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Cos: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn cos(x: anytype) autodiff.@"var".Cos(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Cos(X) = .{ .constant = undefined };

            numeric.cosInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .cos,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.cosInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Tan(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Tan: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn tan(x: anytype) autodiff.@"var".Tan(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Tan(X) = .{ .constant = undefined };

            numeric.tanInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .tan,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.tanInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Asin(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Asin: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn asin(x: anytype) autodiff.@"var".Asin(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Asin(X) = .{ .constant = undefined };

            numeric.asinInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .asin,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.asinInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Acos(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Acos: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn acos(x: anytype) autodiff.@"var".Acos(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Acos(X) = .{ .constant = undefined };

            numeric.acosInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .acos,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.acosInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Atan(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Atan: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn atan(x: anytype) autodiff.@"var".Atan(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Atan(X) = .{ .constant = undefined };

            numeric.atanInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .atan,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.atanInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Atan2(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Atan2: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    comptime if (isVar(X) and isVar(Y) and X != Y)
        @compileError("zsl.autodiff.@\"var\".Atan2: if X and Y are both var types, they must be equal, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return if (comptime isVar(X)) X else Y;
}

pub fn atan2(x: anytype, y: anytype) autodiff.@"var".Atan2(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const R = autodiff.@"var".Atan2(X, Y);

    const x_val = valOf(x);
    const y_val = valOf(y);

    if (!isTracked(x) and !isTracked(y)) {
        var result: R = .{ .constant = undefined };

        numeric.atan2Into(&result.constant, x_val, y_val);

        return result;
    }

    const tape = getTape2(x, y);
    var node: autodiff.Tape(meta.Scalar(R)).Node = .{
        .op = .atan2,
        .left = ensureTracked(tape, x),
        .right = ensureTracked(tape, y),
        .val = undefined,
    };

    numeric.atan2Into(&node.val, x_val, y_val);

    return .{
        .tracked = .{
            .tape = tape,
            .id = tape.pushAssumeCapacity(node),
        },
    };
}

pub fn Sinh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Sinh: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn sinh(x: anytype) autodiff.@"var".Sinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Sinh(X) = .{ .constant = undefined };

            numeric.sinhInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .sinh,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.sinhInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Cosh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Cosh: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn cosh(x: anytype) autodiff.@"var".Cosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Cosh(X) = .{ .constant = undefined };

            numeric.coshInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .cosh,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.coshInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Tanh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Tanh: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn tanh(x: anytype) autodiff.@"var".Tanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Tanh(X) = .{ .constant = undefined };

            numeric.tanhInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .tanh,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.tanhInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Asinh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Asinh: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn asinh(x: anytype) autodiff.@"var".Asinh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Asinh(X) = .{ .constant = undefined };

            numeric.asinhInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .asinh,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.asinhInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Acosh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Acosh: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn acosh(x: anytype) autodiff.@"var".Acosh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Acosh(X) = .{ .constant = undefined };

            numeric.acoshInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .acosh,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.acoshInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

pub fn Atanh(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Atanh: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return X;
}

pub fn atanh(x: anytype) autodiff.@"var".Atanh(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    switch (x) {
        .constant => |cx| {
            var result: autodiff.@"var".Atanh(X) = .{ .constant = undefined };

            numeric.atanhInto(&result.constant, cx);

            return result;
        },
        .tracked => |tx| {
            var node: autodiff.Tape(meta.Scalar(X)).Node = .{
                .op = .atanh,
                .left = tx.id,
                .right = int.highest(usize),
                .val = undefined,
            };

            numeric.atanhInto(&node.val, x.val());

            return .{
                .tracked = .{
                    .tape = tx.tape,
                    .id = tx.tape.pushAssumeCapacity(node),
                },
            };
        },
    }
}

// Utils
inline fn valOf(x: anytype) if (isVar(@TypeOf(x))) meta.Scalar(@TypeOf(x)) else @TypeOf(x) {
    return if (comptime isVar(@TypeOf(x)))
        x.val()
    else
        x;
}

inline fn isTracked(x: anytype) bool {
    return if (comptime isVar(@TypeOf(x)))
        switch (x) {
            .constant => false,
            .tracked => true,
        }
    else
        false;
}

inline fn ensureTracked(tape: anytype, x: anytype) usize {
    return if (comptime isVar(@TypeOf(x)))
        switch (x) {
            .tracked => |t| t.id,
            .constant => |c| tape.pushAssumeCapacity(.{
                .op = .@"var",
                .left = int.highest(usize),
                .right = int.highest(usize),
                .val = c,
            }),
        }
    else
        tape.pushAssumeCapacity(.{
            .op = .@"var",
            .left = int.highest(usize),
            .right = int.highest(usize),
            .val = numeric.cast(meta.Child(@TypeOf(tape)).Numeric, x),
        });
}

inline fn getTape2(x: anytype, y: anytype) if (isVar(@TypeOf(x))) *autodiff.Tape(meta.Scalar(@TypeOf(x))) else *autodiff.Tape(meta.Scalar(@TypeOf(y))) {
    if (comptime isVar(@TypeOf(x))) {
        switch (x) {
            .tracked => |t| return t.tape,
            .constant => {},
        }
    }

    if (comptime isVar(@TypeOf(y))) {
        switch (y) {
            .tracked => |t| return t.tape,
            .constant => {},
        }
    }

    unreachable;
}

inline fn getTape3(x: anytype, y: anytype, z: anytype) if (isVar(@TypeOf(x))) *autodiff.Tape(meta.Scalar(@TypeOf(x))) else if (isVar(@TypeOf(y))) *autodiff.Tape(meta.Scalar(@TypeOf(y))) else *autodiff.Tape(meta.Scalar(@TypeOf(z))) {
    if (comptime isVar(@TypeOf(x))) {
        switch (x) {
            .tracked => |t| return t.tape,
            .constant => {},
        }
    }

    if (comptime isVar(@TypeOf(y))) {
        switch (y) {
            .tracked => |t| return t.tape,
            .constant => {},
        }
    }

    if (comptime isVar(@TypeOf(z))) {
        switch (z) {
            .tracked => |t| return t.tape,
            .constant => {},
        }
    }

    unreachable;
}

inline fn idOf(x: anytype) usize {
    return if (comptime isVar(@TypeOf(x)))
        switch (x) {
            .tracked => |t| t.id,
            .constant => int.highest(usize),
        }
    else
        int.highest(usize);
}
