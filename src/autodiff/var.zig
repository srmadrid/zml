const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

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
                .left = 0,
                .right = 0,
                .val = value,
            });

            return .{ .tracked = .{ .tape = tape, .id = id } };
        }

        // Basic operations
        // pub const Abs = autodiff.@"var".Abs;
        // pub const abs = autodiff.@"var".abs;
        // pub const Abs1 = autodiff.@"var".Abs1;
        // pub const abs1 = autodiff.@"var".abs1;
        // pub const Abs2 = autodiff.@"var".Abs2;
        // pub const abs2 = autodiff.@"var".abs2;
        // pub const Neg = autodiff.@"var".Neg;
        // pub const neg = autodiff.@"var".neg;
        // pub const Re = autodiff.@"var".Re;
        // pub const re = autodiff.@"var".re;
        // pub const Im = autodiff.@"var".Im;
        // pub const im = autodiff.@"var".im;
        // pub const Conj = autodiff.@"var".Conj;
        // pub const conj = autodiff.@"var".conj;
        // pub const Sign = autodiff.@"var".Sign;
        // pub const sign = autodiff.@"var".sign;

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

                    var i: usize = t.id;
                    while (true) : (i -= 1) {
                        const node = &t.tape.nodes[i];
                        const g = node.grad;

                        if (numeric.eq(g, 0))
                            if (i == 0) break else continue;

                        switch (node.op) {
                            .@"var" => {},
                            .add => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    g,
                                );

                                numeric.addInto(
                                    &t.tape.nodes[node.right].grad,
                                    t.tape.nodes[node.right].grad,
                                    g,
                                );
                            },
                            .sub => {
                                numeric.addInto(
                                    &t.tape.nodes[node.left].grad,
                                    t.tape.nodes[node.left].grad,
                                    g,
                                );

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
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    );
                                } else {
                                    numeric.addInto(
                                        &t.tape.nodes[node.right].grad,
                                        t.tape.nodes[node.right].grad,
                                        g,
                                    );
                                }
                            },
                            .min => {
                                if (numeric.lt(t.tape.nodes[node.left].val, t.tape.nodes[node.right].val)) {
                                    numeric.addInto(
                                        &t.tape.nodes[node.left].grad,
                                        t.tape.nodes[node.left].grad,
                                        g,
                                    );
                                } else {
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

                        if (i == 0)
                            break;
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

pub fn Add(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Add: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;

    return Var(numeric.Add(SX, SY));
}

pub fn add(x: anytype, y: anytype) autodiff.@"var".Add(@TypeOf(x), @TypeOf(y)) {
    if (comptime isVar(@TypeOf(x))) {
        if (comptime isVar(@TypeOf(y))) {
            switch (x) {
                .constant => |cx| switch (y) {
                    .constant => |cy| return .{ .constant = numeric.add(cx, cy) },
                    .tracked => |ty| return .{
                        .tracked = .{
                            .tape = ty.tape,
                            .id = ty.tape.pushAssumeCapacity(.{
                                .op = .add,
                                .left = ty.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cx,
                                }),
                                .right = ty.id,
                                .val = numeric.add(
                                    cx,
                                    y.val(),
                                ),
                            }),
                        },
                    },
                },
                .tracked => |tx| switch (y) {
                    .constant => |cy| return .{
                        .tracked = .{
                            .tape = tx.tape,
                            .id = tx.tape.pushAssumeCapacity(.{
                                .op = .add,
                                .left = tx.id,
                                .right = tx.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cy,
                                }),
                                .val = numeric.add(
                                    x.val(),
                                    cy,
                                ),
                            }),
                        },
                    },
                    .tracked => |ty| {
                        std.debug.assert(tx.tape == ty.tape);

                        return .{
                            .tracked = .{
                                .tape = tx.tape,
                                .id = tx.tape.pushAssumeCapacity(.{
                                    .op = .add,
                                    .left = tx.id,
                                    .right = ty.id,
                                    .val = numeric.add(
                                        x.val(),
                                        y.val(),
                                    ),
                                }),
                            },
                        };
                    },
                },
            }
        } else {
            switch (x) {
                .constant => |cx| return .{ .constant = numeric.add(cx, y) },
                .tracked => |tx| return .{
                    .tracked = .{
                        .tape = tx.tape,
                        .id = tx.tape.pushAssumeCapacity(.{
                            .op = .add,
                            .left = tx.id,
                            .right = tx.tape.pushAssumeCapacity(.{
                                .op = .@"var",
                                .left = 0,
                                .right = 0,
                                .val = y,
                            }),
                            .val = numeric.add(
                                x.val(),
                                y,
                            ),
                        }),
                    },
                },
            }
        }
    } else {
        switch (y) {
            .constant => |cy| return .{ .constant = numeric.add(x, cy) },
            .tracked => |ty| return .{
                .tracked = .{
                    .tape = ty.tape,
                    .id = ty.tape.pushAssumeCapacity(.{
                        .op = .add,
                        .left = ty.tape.pushAssumeCapacity(.{
                            .op = .@"var",
                            .left = 0,
                            .right = 0,
                            .val = x,
                        }),
                        .right = ty.id,
                        .val = numeric.add(
                            x,
                            y.val(),
                        ),
                    }),
                },
            },
        }
    }
}

pub fn Sub(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Sub: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;

    return Var(numeric.Sub(SX, SY));
}

pub fn sub(x: anytype, y: anytype) autodiff.@"var".Sub(@TypeOf(x), @TypeOf(y)) {
    if (comptime isVar(@TypeOf(x))) {
        if (comptime isVar(@TypeOf(y))) {
            switch (x) {
                .constant => |cx| switch (y) {
                    .constant => |cy| return .{ .constant = numeric.sub(cx, cy) },
                    .tracked => |ty| return .{
                        .tracked = .{
                            .tape = ty.tape,
                            .id = ty.tape.pushAssumeCapacity(.{
                                .op = .sub,
                                .left = ty.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cx,
                                }),
                                .right = ty.id,
                                .val = numeric.sub(
                                    cx,
                                    y.val(),
                                ),
                            }),
                        },
                    },
                },
                .tracked => |tx| switch (y) {
                    .constant => |cy| return .{
                        .tracked = .{
                            .tape = tx.tape,
                            .id = tx.tape.pushAssumeCapacity(.{
                                .op = .sub,
                                .left = tx.id,
                                .right = tx.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cy,
                                }),
                                .val = numeric.sub(
                                    x.val(),
                                    cy,
                                ),
                            }),
                        },
                    },
                    .tracked => |ty| {
                        std.debug.assert(tx.tape == ty.tape);

                        return .{
                            .tracked = .{
                                .tape = tx.tape,
                                .id = tx.tape.pushAssumeCapacity(.{
                                    .op = .sub,
                                    .left = tx.id,
                                    .right = ty.id,
                                    .val = numeric.sub(
                                        x.val(),
                                        y.val(),
                                    ),
                                }),
                            },
                        };
                    },
                },
            }
        } else {
            switch (x) {
                .constant => |cx| return .{ .constant = numeric.sub(cx, y) },
                .tracked => |tx| return .{
                    .tracked = .{
                        .tape = tx.tape,
                        .id = tx.tape.pushAssumeCapacity(.{
                            .op = .sub,
                            .left = tx.id,
                            .right = tx.tape.pushAssumeCapacity(.{
                                .op = .@"var",
                                .left = 0,
                                .right = 0,
                                .val = y,
                            }),
                            .val = numeric.sub(
                                x.val(),
                                y,
                            ),
                        }),
                    },
                },
            }
        }
    } else {
        switch (y) {
            .constant => |cy| return .{ .constant = numeric.sub(x, cy) },
            .tracked => |ty| return .{
                .tracked = .{
                    .tape = ty.tape,
                    .id = ty.tape.pushAssumeCapacity(.{
                        .op = .sub,
                        .left = ty.tape.pushAssumeCapacity(.{
                            .op = .@"var",
                            .left = 0,
                            .right = 0,
                            .val = x,
                        }),
                        .right = ty.id,
                        .val = numeric.sub(
                            x,
                            y.val(),
                        ),
                    }),
                },
            },
        }
    }
}

pub fn Mul(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Mul: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;

    return Var(numeric.Mul(SX, SY));
}

pub fn mul(x: anytype, y: anytype) autodiff.@"var".Mul(@TypeOf(x), @TypeOf(y)) {
    if (comptime isVar(@TypeOf(x))) {
        if (comptime isVar(@TypeOf(y))) {
            switch (x) {
                .constant => |cx| switch (y) {
                    .constant => |cy| return .{ .constant = numeric.mul(cx, cy) },
                    .tracked => |ty| return .{
                        .tracked = .{
                            .tape = ty.tape,
                            .id = ty.tape.pushAssumeCapacity(.{
                                .op = .mul,
                                .left = ty.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cx,
                                }),
                                .right = ty.id,
                                .val = numeric.mul(cx, y.val()),
                            }),
                        },
                    },
                },
                .tracked => |tx| switch (y) {
                    .constant => |cy| return .{
                        .tracked = .{
                            .tape = tx.tape,
                            .id = tx.tape.pushAssumeCapacity(.{
                                .op = .mul,
                                .left = tx.id,
                                .right = tx.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cy,
                                }),
                                .val = numeric.mul(x.val(), cy),
                            }),
                        },
                    },
                    .tracked => |ty| {
                        std.debug.assert(tx.tape == ty.tape);
                        return .{
                            .tracked = .{
                                .tape = tx.tape,
                                .id = tx.tape.pushAssumeCapacity(.{
                                    .op = .mul,
                                    .left = tx.id,
                                    .right = ty.id,
                                    .val = numeric.mul(x.val(), y.val()),
                                }),
                            },
                        };
                    },
                },
            }
        } else {
            switch (x) {
                .constant => |cx| return .{ .constant = numeric.mul(cx, y) },
                .tracked => |tx| return .{
                    .tracked = .{
                        .tape = tx.tape,
                        .id = tx.tape.pushAssumeCapacity(.{
                            .op = .mul,
                            .left = tx.id,
                            .right = tx.tape.pushAssumeCapacity(.{
                                .op = .@"var",
                                .left = 0,
                                .right = 0,
                                .val = y,
                            }),
                            .val = numeric.mul(x.val(), y),
                        }),
                    },
                },
            }
        }
    } else {
        switch (y) {
            .constant => |cy| return .{ .constant = numeric.mul(x, cy) },
            .tracked => |ty| return .{
                .tracked = .{
                    .tape = ty.tape,
                    .id = ty.tape.pushAssumeCapacity(.{
                        .op = .mul,
                        .left = ty.tape.pushAssumeCapacity(.{
                            .op = .@"var",
                            .left = 0,
                            .right = 0,
                            .val = x,
                        }),
                        .right = ty.id,
                        .val = numeric.mul(x, y.val()),
                    }),
                },
            },
        }
    }
}

pub fn Fma(comptime X: type, comptime Y: type, comptime Z: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or !meta.isNumeric(Z) or (!isVar(X) and !isVar(Y) and !isVar(Z)))
        @compileError("zsl.autodiff.@\"var\".Fma: at least one of X, Y or Z must be a var type, the others must be numeric or var types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n\tZ = " ++ @typeName(Z) ++ "\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;
    const SZ: type = if (isVar(Z)) meta.Scalar(Z) else Z;

    return Var(numeric.Fma(SX, SY, SZ));
}

pub fn fma(x: anytype, y: anytype, z: anytype) autodiff.@"var".Fma(@TypeOf(x), @TypeOf(y), @TypeOf(z)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);
    const Z = @TypeOf(z);

    const vx = if (comptime isVar(X)) x.val() else x;
    const vy = if (comptime isVar(Y)) y.val() else y;
    const vz = if (comptime isVar(Z)) z.val() else z;
    const vr = numeric.fma(vx, vy, vz);

    const x_tracked = (comptime isVar(X)) and x == .tracked;
    const y_tracked = (comptime isVar(Y)) and y == .tracked;
    const z_tracked = (comptime isVar(Z)) and z == .tracked;

    if (!x_tracked and !y_tracked and !z_tracked) {
        return .{ .constant = vr };
    }

    const tape = if (x_tracked) x.tracked.tape else if (y_tracked) y.tracked.tape else z.tracked.tape;

    if (x_tracked and y_tracked) std.debug.assert(x.tracked.tape == y.tracked.tape);
    if (y_tracked and z_tracked) std.debug.assert(y.tracked.tape == z.tracked.tape);
    if (x_tracked and z_tracked) std.debug.assert(x.tracked.tape == z.tracked.tape);

    const idx = if (x_tracked) x.tracked.id else tape.pushAssumeCapacity(.{ .op = .@"var", .left = 0, .right = 0, .val = vx });
    const idy = if (y_tracked) y.tracked.id else tape.pushAssumeCapacity(.{ .op = .@"var", .left = 0, .right = 0, .val = vy });
    const idz = if (z_tracked) z.tracked.id else tape.pushAssumeCapacity(.{ .op = .@"var", .left = 0, .right = 0, .val = vz });

    const idm = tape.pushAssumeCapacity(.{
        .op = .mul,
        .left = idx,
        .right = idy,
        .val = numeric.mul(vx, vy),
    });

    return .{ .tracked = .{
        .tape = tape,
        .id = tape.pushAssumeCapacity(.{
            .op = .add,
            .left = idm,
            .right = idz,
            .val = vr,
        }),
    } };
}

pub fn Div(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Div: at least one of X or Y must be a var type, the other must be a numeric or a var type, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;

    return Var(numeric.Div(SX, SY));
}

pub fn div(x: anytype, y: anytype) autodiff.@"var".Div(@TypeOf(x), @TypeOf(y)) {
    if (comptime isVar(@TypeOf(x))) {
        if (comptime isVar(@TypeOf(y))) {
            switch (x) {
                .constant => |cx| switch (y) {
                    .constant => |cy| return .{ .constant = numeric.div(cx, cy) },
                    .tracked => |ty| return .{
                        .tracked = .{
                            .tape = ty.tape,
                            .id = ty.tape.pushAssumeCapacity(.{
                                .op = .div,
                                .left = ty.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cx,
                                }),
                                .right = ty.id,
                                .val = numeric.div(cx, y.val()),
                            }),
                        },
                    },
                },
                .tracked => |tx| switch (y) {
                    .constant => |cy| return .{
                        .tracked = .{
                            .tape = tx.tape,
                            .id = tx.tape.pushAssumeCapacity(.{
                                .op = .div,
                                .left = tx.id,
                                .right = tx.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cy,
                                }),
                                .val = numeric.div(x.val(), cy),
                            }),
                        },
                    },
                    .tracked => |ty| {
                        std.debug.assert(tx.tape == ty.tape);
                        return .{
                            .tracked = .{
                                .tape = tx.tape,
                                .id = tx.tape.pushAssumeCapacity(.{
                                    .op = .div,
                                    .left = tx.id,
                                    .right = ty.id,
                                    .val = numeric.div(x.val(), y.val()),
                                }),
                            },
                        };
                    },
                },
            }
        } else {
            switch (x) {
                .constant => |cx| return .{ .constant = numeric.div(cx, y) },
                .tracked => |tx| return .{
                    .tracked = .{
                        .tape = tx.tape,
                        .id = tx.tape.pushAssumeCapacity(.{
                            .op = .div,
                            .left = tx.id,
                            .right = tx.tape.pushAssumeCapacity(.{
                                .op = .@"var",
                                .left = 0,
                                .right = 0,
                                .val = y,
                            }),
                            .val = numeric.div(x.val(), y),
                        }),
                    },
                },
            }
        }
    } else {
        switch (y) {
            .constant => |cy| return .{ .constant = numeric.div(x, cy) },
            .tracked => |ty| return .{
                .tracked = .{
                    .tape = ty.tape,
                    .id = ty.tape.pushAssumeCapacity(.{
                        .op = .div,
                        .left = ty.tape.pushAssumeCapacity(.{
                            .op = .@"var",
                            .left = 0,
                            .right = 0,
                            .val = x,
                        }),
                        .right = ty.id,
                        .val = numeric.div(x, y.val()),
                    }),
                },
            },
        }
    }
}

pub fn eq(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".eq: at least one of x or y must be a var, the other must be a numeric or a var, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.eq(
        if (comptime isVar(X)) x.val() else x,
        if (comptime isVar(Y)) y.val() else y,
    );
}

pub fn ne(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".ne: at least one of x or y must be a var, the other must be a numeric or a var, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.ne(
        if (comptime isVar(X)) x.val() else x,
        if (comptime isVar(Y)) y.val() else y,
    );
}

pub fn lt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".lt: at least one of x or y must be a var, the other must be a numeric or a var, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.lt(
        if (comptime isVar(X)) x.val() else x,
        if (comptime isVar(Y)) y.val() else y,
    );
}

pub fn le(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".le: at least one of x or y must be a var, the other must be a numeric or a var, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.le(
        if (comptime isVar(X)) x.val() else x,
        if (comptime isVar(Y)) y.val() else y,
    );
}

pub fn gt(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".gt: at least one of x or y must be a var, the other must be a numeric or a var, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.gt(
        if (comptime isVar(X)) x.val() else x,
        if (comptime isVar(Y)) y.val() else y,
    );
}

pub fn ge(x: anytype, y: anytype) bool {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".ge: at least one of x or y must be a var, the other must be a numeric or a var, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    return numeric.ge(
        if (comptime isVar(X)) x.val() else x,
        if (comptime isVar(Y)) y.val() else y,
    );
}

pub fn Max(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Max: at least one of X or Y must be a var type...\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;

    return Var(numeric.Max(SX, SY));
}

pub fn max(x: anytype, y: anytype) autodiff.@"var".Max(@TypeOf(x), @TypeOf(y)) {
    if (comptime isVar(@TypeOf(x))) {
        if (comptime isVar(@TypeOf(y))) {
            switch (x) {
                .constant => |cx| switch (y) {
                    .constant => |cy| return .{ .constant = numeric.max(cx, cy) },
                    .tracked => |ty| return .{
                        .tracked = .{
                            .tape = ty.tape,
                            .id = ty.tape.pushAssumeCapacity(.{
                                .op = .max,
                                .left = ty.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cx,
                                }),
                                .right = ty.id,
                                .val = numeric.max(cx, y.val()),
                            }),
                        },
                    },
                },
                .tracked => |tx| switch (y) {
                    .constant => |cy| return .{
                        .tracked = .{
                            .tape = tx.tape,
                            .id = tx.tape.pushAssumeCapacity(.{
                                .op = .max,
                                .left = tx.id,
                                .right = tx.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cy,
                                }),
                                .val = numeric.max(x.val(), cy),
                            }),
                        },
                    },
                    .tracked => |ty| {
                        std.debug.assert(tx.tape == ty.tape);
                        return .{
                            .tracked = .{
                                .tape = tx.tape,
                                .id = tx.tape.pushAssumeCapacity(.{
                                    .op = .max,
                                    .left = tx.id,
                                    .right = ty.id,
                                    .val = numeric.max(x.val(), y.val()),
                                }),
                            },
                        };
                    },
                },
            }
        } else {
            switch (x) {
                .constant => |cx| return .{ .constant = numeric.max(cx, y) },
                .tracked => |tx| return .{
                    .tracked = .{
                        .tape = tx.tape,
                        .id = tx.tape.pushAssumeCapacity(.{
                            .op = .max,
                            .left = tx.id,
                            .right = tx.tape.pushAssumeCapacity(.{
                                .op = .@"var",
                                .left = 0,
                                .right = 0,
                                .val = y,
                            }),
                            .val = numeric.max(x.val(), y),
                        }),
                    },
                },
            }
        }
    } else {
        switch (y) {
            .constant => |cy| return .{ .constant = numeric.max(x, cy) },
            .tracked => |ty| return .{
                .tracked = .{
                    .tape = ty.tape,
                    .id = ty.tape.pushAssumeCapacity(.{
                        .op = .max,
                        .left = ty.tape.pushAssumeCapacity(.{
                            .op = .@"var",
                            .left = 0,
                            .right = 0,
                            .val = x,
                        }),
                        .right = ty.id,
                        .val = numeric.max(x, y.val()),
                    }),
                },
            },
        }
    }
}

pub fn Min(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isNumeric(X) or !meta.isNumeric(Y) or (!isVar(X) and !isVar(Y)))
        @compileError("zsl.autodiff.@\"var\".Min: at least one of X or Y must be a var type...\n");

    const SX: type = if (isVar(X)) meta.Scalar(X) else X;
    const SY: type = if (isVar(Y)) meta.Scalar(Y) else Y;

    return Var(numeric.Min(SX, SY));
}

pub fn min(x: anytype, y: anytype) autodiff.@"var".Min(@TypeOf(x), @TypeOf(y)) {
    if (comptime isVar(@TypeOf(x))) {
        if (comptime isVar(@TypeOf(y))) {
            switch (x) {
                .constant => |cx| switch (y) {
                    .constant => |cy| return .{ .constant = numeric.min(cx, cy) },
                    .tracked => |ty| return .{
                        .tracked = .{
                            .tape = ty.tape,
                            .id = ty.tape.pushAssumeCapacity(.{
                                .op = .min,
                                .left = ty.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cx,
                                }),
                                .right = ty.id,
                                .val = numeric.min(cx, y.val()),
                            }),
                        },
                    },
                },
                .tracked => |tx| switch (y) {
                    .constant => |cy| return .{
                        .tracked = .{
                            .tape = tx.tape,
                            .id = tx.tape.pushAssumeCapacity(.{
                                .op = .min,
                                .left = tx.id,
                                .right = tx.tape.pushAssumeCapacity(.{
                                    .op = .@"var",
                                    .left = 0,
                                    .right = 0,
                                    .val = cy,
                                }),
                                .val = numeric.min(x.val(), cy),
                            }),
                        },
                    },
                    .tracked => |ty| {
                        std.debug.assert(tx.tape == ty.tape);
                        return .{
                            .tracked = .{
                                .tape = tx.tape,
                                .id = tx.tape.pushAssumeCapacity(.{
                                    .op = .min,
                                    .left = tx.id,
                                    .right = ty.id,
                                    .val = numeric.min(x.val(), y.val()),
                                }),
                            },
                        };
                    },
                },
            }
        } else {
            switch (x) {
                .constant => |cx| return .{ .constant = numeric.min(cx, y) },
                .tracked => |tx| return .{
                    .tracked = .{
                        .tape = tx.tape,
                        .id = tx.tape.pushAssumeCapacity(.{
                            .op = .min,
                            .left = tx.id,
                            .right = tx.tape.pushAssumeCapacity(.{
                                .op = .@"var",
                                .left = 0,
                                .right = 0,
                                .val = y,
                            }),
                            .val = numeric.min(x.val(), y),
                        }),
                    },
                },
            }
        }
    } else {
        switch (y) {
            .constant => |cy| return .{ .constant = numeric.min(x, cy) },
            .tracked => |ty| return .{
                .tracked = .{
                    .tape = ty.tape,
                    .id = ty.tape.pushAssumeCapacity(.{
                        .op = .min,
                        .left = ty.tape.pushAssumeCapacity(.{
                            .op = .@"var",
                            .left = 0,
                            .right = 0,
                            .val = x,
                        }),
                        .right = ty.id,
                        .val = numeric.min(x, y.val()),
                    }),
                },
            },
        }
    }
}

pub fn Exp(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Exp: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Var(numeric.Exp(meta.Scalar(X)));
}

pub fn exp(x: anytype) autodiff.@"var".Exp(@TypeOf(x)) {
    switch (x) {
        .constant => |cx| return .{ .constant = numeric.exp(cx) },
        .tracked => |tx| return .{
            .tracked = .{
                .tape = tx.tape,
                .id = tx.tape.pushAssumeCapacity(.{
                    .op = .exp,
                    .left = tx.id,
                    .right = 0,
                    .val = numeric.exp(x.val()),
                }),
            },
        },
    }
}

pub fn Ln(comptime X: type) type {
    comptime if (!meta.isNumeric(X) or !isVar(X))
        @compileError("zsl.autodiff.@\"var\".Ln: X must be a var type, got\n\tX = " ++ @typeName(X) ++ "\n");

    return Var(numeric.Ln(meta.Scalar(X)));
}

pub fn ln(x: anytype) autodiff.@"var".Ln(@TypeOf(x)) {
    switch (x) {
        .constant => |cx| return .{ .constant = numeric.ln(cx) },
        .tracked => |tx| return .{
            .tracked = .{
                .tape = tx.tape,
                .id = tx.tape.pushAssumeCapacity(.{
                    .op = .ln,
                    .left = tx.id,
                    .right = 0,
                    .val = numeric.ln(x.val()),
                }),
            },
        },
    }
}
