const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const EnsureArray = types.EnsureArray;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");

const dense = @import("dense.zig");
const strided = @import("strided.zig");
const sparse = @import("sparse.zig");

pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    opts: struct {
        order: ?types.Order = null,
        // eventually will add axis, axes and keepdims opts, to apply2 as well
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), ReturnType1(op, Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X))
        @compileError("apply1: x must be an array, got " ++ @typeName(X));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 1 and @typeInfo(@TypeOf(op)).@"fn".params.len != 2))
        @compileError("apply1: op must be a function of one argument, or a function of two arguments with the second argument being a context, got " ++ @typeName(@TypeOf(op)));

    switch (comptime types.arrayType(@TypeOf(x))) {
        .dense => return dense.apply1(allocator, x, op, .{ .order = opts.order }, ctx),
        .strided => return strided.apply1(allocator, x, op, .{ .order = opts.order }, ctx),
        .sparse => @compileError("apply1 not implemented for sparse arrays yet"),
        .numeric => unreachable,
    }
}

pub fn apply1_(
    o: anytype,
    x: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("apply1_: o must be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O))
        @compileError("apply1_: o must be an array, got " ++ @typeName(O));

    comptime if (@typeInfo(@TypeOf(op_)) != .@"fn" or (@typeInfo(@TypeOf(op_)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op_)).@"fn".params.len != 3))
        @compileError("apply1_: op_ must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op_)));

    if (comptime !types.isArray(X)) {
        switch (comptime types.arrayType(O)) {
            .dense => return dense.apply1_(o, x, op_, ctx),
            .strided => return strided.apply1_(o, x, op_, ctx),
            .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(X)) {
                .dense => return dense.apply1_(o, x, op_, ctx), // dense dense apply1_
                .strided => return strided.apply1_(o, x, op_, ctx), // dense strided apply1_
                .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(X)) {
                .dense => return strided.apply1_(o, x, op_, ctx), // strided dense apply1_
                .strided => return strided.apply1_(o, x, op_, ctx), // strided strided apply1_
                .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    }
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isArray(Y))
        @compileError("apply2: at least one of x or y must be an array, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isArray(X)) {
        switch (comptime types.arrayType(Y)) {
            .dense => return dense.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
            .strided => return strided.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
            .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isArray(Y)) {
        switch (comptime types.arrayType(X)) {
            .dense => return dense.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
            .strided => return strided.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
            .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.arrayType(X)) {
            .dense => switch (comptime types.arrayType(Y)) {
                .dense => return dense.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
                .strided => return strided.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
                .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(Y)) {
                .dense => return strided.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
                .strided => return strided.apply2(allocator, x, y, op, .{ .order = opts.order }, ctx),
                .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2 not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    }
}

pub fn apply2_(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("apply2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O))
        @compileError("apply2_: o must be an array, got " ++ @typeName(O));

    comptime if (@typeInfo(@TypeOf(op_)) != .@"fn" or (@typeInfo(@TypeOf(op_)).@"fn".params.len != 3 and @typeInfo(@TypeOf(op_)).@"fn".params.len != 4))
        @compileError("apply2_: op must be a function of three arguments, or a function of four arguments with the fourth argument being a context, got " ++ @typeName(@TypeOf(op_)));

    if (comptime !types.isArray(X)) {
        if (comptime !types.isArray(Y)) {
            switch (comptime types.arrayType(O)) {
                .dense => return dense.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            }
        }

        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(Y)) {
                .dense => return dense.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(Y)) {
                .dense => return strided.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isArray(Y)) {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(X)) {
                .dense => return dense.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(X)) {
                .dense => return strided.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(X)) {
                .dense => switch (comptime types.arrayType(Y)) {
                    .dense => return dense.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .strided => switch (comptime types.arrayType(Y)) {
                    .dense => return strided.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(X)) {
                .dense => switch (comptime types.arrayType(Y)) {
                    .dense => return strided.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .strided => switch (comptime types.arrayType(Y)) {
                    .dense => return strided.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    }
}

// Basic operations
pub inline fn abs(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Scalar(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{
                .allocator = .{ .type = std.mem.Allocator, .required = true },
                .copy = .{ .type = bool, .required = false },
            },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.abs,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn abs_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .copy = .{ .type = bool, .required = false },
                    },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            types.validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.abs_, ctx);
}

pub inline fn abs2(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Scalar(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.abs2,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn abs2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.abs2_, ctx);
}

// Exponential functions
pub inline fn exp(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn exp_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.exp_, ctx);
}

pub inline fn exp10(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp10,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn exp10_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.exp10_, ctx);
}

pub inline fn exp2(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp2,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn exp2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.exp2_, ctx);
}

pub inline fn exp10m1(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp10m1,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn exp10m1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.exp10m1_, ctx);
}

pub inline fn exp2m1(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp2m1,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn exp2m1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.exp2m1_, ctx);
}

pub inline fn expm1(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.expm1,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn expm1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.expm1_, ctx);
}

pub inline fn log(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn log_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.log_, ctx);
}

pub inline fn log10(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log10,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn log10_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.log10_, ctx);
}

pub inline fn log2(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log2,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn log2_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.log2_, ctx);
}

pub inline fn log10p1(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log10p1,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn log10p1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.log10p1_, ctx);
}

pub inline fn log2p1(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log2p1,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn log2p1_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.log2p1_, ctx);
}

pub inline fn log1p(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log1p,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn log1p_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.log1p_, ctx);
}

// Power functions
pub inline fn pow(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.pow,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn pow_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            types.validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.pow_, ctx);
}

pub inline fn sqrt(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sqrt,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn sqrt_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.sqrt_, ctx);
}

pub inline fn cbrt(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cbrt,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn cbrt_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.cbrt_, ctx);
}

pub inline fn hypot(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.hypot,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn hypot_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            types.validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.hypot_, ctx);
}

// Trigonometric functions
pub inline fn sin(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sin,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn sin_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.sin_, ctx);
}

pub inline fn cos(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cos,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn cos_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.cos_, ctx);
}

pub inline fn tan(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.tan,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn tan_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.tan_, ctx);
}

pub inline fn asin(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.asin,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn asin_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.asin_, ctx);
}

pub inline fn acos(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.acos,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn acos_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.acos_, ctx);
}

pub inline fn atan(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.atan,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn atan_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.atan_, ctx);
}

pub inline fn atan2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.atan2,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn atan2_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            types.validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.atan2_, ctx);
}

pub inline fn sinpi(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sinpi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn sinpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.sinpi_, ctx);
}

pub inline fn cospi(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cospi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn cospi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.cospi_, ctx);
}

pub inline fn tanpi(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.tanpi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn tanpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.tanpi_, ctx);
}

pub inline fn asinpi(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.asinpi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn asinpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.asinpi_, ctx);
}

pub inline fn acospi(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.acospi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn acospi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.acospi_, ctx);
}

pub inline fn atanpi(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.atanpi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn atanpi_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.atanpi_, ctx);
}

pub inline fn atan2pi(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.atan2pi,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn atan2pi_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            types.validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.atan2pi_, ctx);
}

// Hyperbolic functions
pub inline fn sinh(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sinh,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn sinh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.sinh_, ctx);
}

pub inline fn cosh(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cosh,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn cosh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.cosh_, ctx);
}

pub inline fn tanh(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.tanh,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn tanh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.tanh_, ctx);
}

pub inline fn asinh(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.asinh,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn asinh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.asinh_, ctx);
}

pub inline fn acosh(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.acosh,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn acosh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.acosh_, ctx);
}

pub inline fn atanh(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.atanh,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn atanh_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.atanh_, ctx);
}

// Error and gamma functions
pub inline fn erf(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.erf,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn erf_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.erf_, ctx);
}

pub inline fn erfc(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.erfc,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn erfc_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.erfc_, ctx);
}

pub inline fn gamma(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.gamma,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn gamma_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.gamma_, ctx);
}

pub inline fn lgamma(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.lgamma,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn lgamma_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.lgamma_, ctx);
}

// Nearest integer operations
pub inline fn ceil(
    allocator: std.mem.Allocator,
    x: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Numeric(@TypeOf(x))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.ceil,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn ceil_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.ceil_, ctx);
}

pub inline fn add(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            types.validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.add,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn add_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        }
    };

    return apply2_(o, x, y, ops.add_, ctx);
}

pub inline fn sub(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            types.validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.sub,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn sub_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        }
    };

    return apply2_(o, x, y, ops.sub_, ctx);
}

pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            types.validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.add,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn mul_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        }
    };

    return apply2_(o, x, y, ops.mul_, ctx);
}

pub inline fn div(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.div,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn div_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(C)) {
            // Both types are arbitrary precision
            if (O == C) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            // Only the output is arbitrary precision, so we need the output's
            // allocator to perform the casting
            types.validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            types.validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.div_, ctx);
}

pub inline fn eq(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    return apply2(allocator, x, y, ops.eq, .{ .order = opts.order }, .{});
}

pub inline fn eq_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.eq_, ctx);
}

pub inline fn ne(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    return apply2(allocator, x, y, ops.ne, .{ .order = opts.order }, .{});
}

pub inline fn ne_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.ne_, ctx);
}

pub inline fn lt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    return apply2(allocator, x, y, ops.lt, .{ .order = opts.order }, .{});
}

pub inline fn lt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.lt_, ctx);
}

pub inline fn le(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    return apply2(allocator, x, y, ops.le, .{ .order = opts.order }, .{});
}

pub inline fn le_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.le_, ctx);
}

pub inline fn gt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    return apply2(allocator, x, y, ops.gt, .{ .order = opts.order }, .{});
}

pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.gt_, ctx);
}

pub inline fn ge(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    return apply2(allocator, x, y, ops.ge, .{ .order = opts.order }, .{});
}

pub inline fn ge_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.ge_, ctx);
}

pub inline fn max(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.max,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn max_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.max_, ctx);
}

pub inline fn min(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    opts: struct {
        order: ?types.Order = null,
    },
    ctx: anytype,
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.min,
        .{ .order = opts.order },
        ctx,
    );
}

pub inline fn min_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.min_, ctx);
}
