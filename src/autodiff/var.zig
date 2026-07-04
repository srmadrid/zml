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
                .left = int.maxVal(usize),
                .right = int.maxVal(usize),
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
        // pub const Pow = autodiff.@"var".Pow;
        // pub const pow = autodiff.@"var".pow;
        // pub const Sqrt = autodiff.@"var".Sqrt;
        // pub const sqrt = autodiff.@"var".sqrt;
        // pub const Cbrt = autodiff.@"var".Cbrt;
        // pub const cbrt = autodiff.@"var".cbrt;
        // pub const Hypot = autodiff.@"var".Hypot;
        // pub const hypot = autodiff.@"var".hypot;

        // Trigonometric functions
        // pub const Sin = autodiff.@"var".Sin;
        // pub const sin = autodiff.@"var".sin;
        // pub const Cos = autodiff.@"var".Cos;
        // pub const cos = autodiff.@"var".cos;
        // pub const Tan = autodiff.@"var".Tan;
        // pub const tan = autodiff.@"var".tan;
        // pub const Asin = autodiff.@"var".Asin;
        // pub const asin = autodiff.@"var".asin;
        // pub const Acos = autodiff.@"var".Acos;
        // pub const acos = autodiff.@"var".acos;
        // pub const Atan = autodiff.@"var".Atan;
        // pub const atan = autodiff.@"var".atan;
        // pub const Atan2 = autodiff.@"var".Atan2;
        // pub const atan2 = autodiff.@"var".atan2;

        // Hyperbolic functions
        // pub const Sinh = autodiff.@"var".Sinh;
        // pub const sinh = autodiff.@"var".sinh;
        // pub const Cosh = autodiff.@"var".Cosh;
        // pub const cosh = autodiff.@"var".cosh;
        // pub const Tanh = autodiff.@"var".Tanh;
        // pub const tanh = autodiff.@"var".tanh;
        // pub const Asinh = autodiff.@"var".Asinh;
        // pub const asinh = autodiff.@"var".asinh;
        // pub const Acosh = autodiff.@"var".Acosh;
        // pub const acosh = autodiff.@"var".acosh;
        // pub const Atanh = autodiff.@"var".Atanh;
        // pub const atanh = autodiff.@"var".atanh;

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
                                if (!numeric.eq(node.val, 0)) {
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        numeric.div(numeric.mul(g, t.tape.nodes[node.left].val), node.val),
                                    );
                                }
                            },
                            .abs1 => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.sign(t.tape.nodes[node.left].val)),
                                );
                            },
                            .abs2 => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, numeric.mul(numeric.two(N), t.tape.nodes[node.left].val)),
                                );
                            },
                            .neg => {
                                numeric.subInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    g,
                                );
                            },
                            .re => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.re(g),
                                );
                            },
                            .im => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.im(g),
                                );
                            },
                            .conj => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.conj(g),
                                );
                            },
                            .sign => {},
                            .add => {
                                if (node.left != int.maxVal(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    );

                                if (node.right != int.maxVal(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                            },
                            .sub => {
                                if (node.left != int.maxVal(usize))
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    );

                                if (node.right != int.maxVal(usize))
                                    numeric.subInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                            },
                            .mul => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, t.tape.nodes[node.right].val),
                                );

                                numeric.addInto(
                                    &t.tape.nodes[node.right].grad,
                                    t.tape.nodes[node.right].grad,
                                    numeric.mul(g, t.tape.nodes[node.left].val),
                                );
                            },
                            .div => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(
                                        g,
                                        t.tape.nodes[node.right].val,
                                    ),
                                );

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
                                if (numeric.gt(t.tape.nodes[node.left].val, t.tape.nodes[node.right].val)) {
                                    if (node.left != int.maxVal(usize))
                                        numeric.addInto(
                                            &t.tape.nodes[node.left].grad,
                                            t.tape.nodes[node.left].grad,
                                            g,
                                        );
                                } else {
                                    if (node.right != int.maxVal(usize))
                                        numeric.addInto(
                                            &t.tape.nodes[node.right].grad,
                                            t.tape.nodes[node.right].grad,
                                            g,
                                        );
                                }
                            },
                            .min => {
                                if (numeric.lt(t.tape.nodes[node.left].val, t.tape.nodes[node.right].val)) {
                                    if (node.left != int.maxVal(usize))
                                        numeric.addInto(
                                            &t.tape.nodes[node.left].grad,
                                            t.tape.nodes[node.left].grad,
                                            g,
                                        );
                                } else {
                                    if (node.right != int.maxVal(usize))
                                        numeric.addInto(
                                            &t.tape.nodes[node.right].grad,
                                            t.tape.nodes[node.right].grad,
                                            g,
                                        );
                                }
                            },
                            .exp => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.mul(g, node.val),
                                );
                            },
                            .ln => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    numeric.div(g, t.tape.nodes[node.left].val),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .right = int.maxVal(usize),
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
                .left = int.maxVal(usize),
                .right = int.maxVal(usize),
                .val = c,
            }),
        }
    else
        tape.pushAssumeCapacity(.{
            .op = .@"var",
            .left = int.maxVal(usize),
            .right = int.maxVal(usize),
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
            .constant => int.maxVal(usize),
        }
    else
        int.maxVal(usize);
}
