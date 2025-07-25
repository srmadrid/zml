const std = @import("std");

// Maybe in place versions (_) also require allocator if the coerced type of the inputs is arbitrary precision

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;
const validateContext = types.validateContext;
const getFieldOrDefault = types.getFieldOrDefault;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");
const Array = array.Array;

const dense = @import("dense.zig");
const strided = @import("strided.zig");

pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    options: struct {
        order: ?array.Order = null,
        // eventually will add axis, axes and keepdims options, to apply2 as well
    },
    ctx: anytype,
) !Array(ReturnType1(op, Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isSlice(X))
        @compileError("apply1: x must be an array or slice, got " ++ @typeName(X));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 1 and @typeInfo(@TypeOf(op)).@"fn".params.len != 2))
        @compileError("apply1: op must be a function of one argument, or a function of two arguments with the second argument being a context, got " ++ @typeName(@TypeOf(op)));

    switch (x.flags.storage) {
        .dense => return dense.apply1(Numeric(X), allocator, x, op, options.order orelse x.flags.order, ctx),
        .strided => return strided.apply1(Numeric(X), allocator, x, op, options.order orelse x.flags.order, ctx),
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

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1_: o must be an array or slice, got " ++ @typeName(O));

    comptime if (@typeInfo(@TypeOf(op_)) != .@"fn" or (@typeInfo(@TypeOf(op_)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op_)).@"fn".params.len != 3))
        @compileError("apply1_: op_ must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op_)));

    if (comptime !types.isArray(X) and !types.isSlice(X)) {
        switch (o.flags.storage) {
            .dense => return dense.apply1_(Numeric(O), o, Numeric(X), x, op_, ctx),
            .strided => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, ctx),
        }
    } else {
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => return dense.apply1_(Numeric(O), o, Numeric(X), x, op_, ctx),
                .strided => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, ctx),
            },
            .strided => switch (x.flags.storage) {
                .dense => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, ctx),
                .strided => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, ctx),
            },
        }
    }
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isSlice(X) and
        !types.isArray(Y) and !types.isSlice(Y))
        @compileError("apply2: at least one of x or y must be an array or slice, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op)) != .@"fn" or (@typeInfo(@TypeOf(op)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op)).@"fn".params.len != 3))
        @compileError("apply2: op must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op)));

    if (comptime !types.isArray(X) and !types.isSlice(X)) {
        // x is a scalar, only consider y's storage
        switch (y.flags.storage) {
            .dense => return dense.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse y.flags.order, ctx),
            .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse y.flags.order, ctx),
        }
    } else if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
        // y is a scalar, only consider x's storage
        switch (x.flags.storage) {
            .dense => return dense.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse x.flags.order, ctx),
            .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse x.flags.order, ctx),
        }
    } else {
        switch (x.flags.storage) {
            .dense => switch (y.flags.storage) {
                .dense => return dense.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse x.flags.order.resolve2(y.flags.order), ctx),
                .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse x.flags.order.resolve2(y.flags.order), ctx),
            },
            .strided => switch (y.flags.storage) {
                .dense => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse x.flags.order.resolve2(y.flags.order), ctx),
                .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, options.order orelse x.flags.order.resolve2(y.flags.order), ctx),
            },
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
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply2_: o must be an array or slice, got " ++ @typeName(O));

    comptime if (@typeInfo(@TypeOf(op_)) != .@"fn" or (@typeInfo(@TypeOf(op_)).@"fn".params.len != 3 and @typeInfo(@TypeOf(op_)).@"fn".params.len != 4))
        @compileError("apply2_: op must be a function of three arguments, or a function of four arguments with the fourth argument being a context, got " ++ @typeName(@TypeOf(op_)));

    if (comptime !types.isArray(X) and !types.isSlice(X)) {
        if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
            // x and y are scalars, only consider o's storage
            switch (o.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
            }
        }

        // x is a scalar, only consider o and y's storage
        switch (o.flags.storage) {
            .dense => switch (y.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
            },
            .strided => switch (y.flags.storage) {
                .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
            },
        }
    } else if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
        // y is a scalar, only consider o and x's storage
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
            },
            .strided => switch (x.flags.storage) {
                .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
            },
        }
    } else {
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => switch (y.flags.storage) {
                    .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                },
                .strided => switch (y.flags.storage) {
                    .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                },
            },
            .strided => switch (x.flags.storage) {
                .dense => switch (y.flags.storage) {
                    .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                },
                .strided => switch (y.flags.storage) {
                    .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, ctx),
                },
            },
        }
    }
}

// Basic operations
pub inline fn abs(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Scalar(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{
                .allocator = .{ .type = std.mem.Allocator, .required = true },
                .copy = .{ .type = bool, .required = false },
            },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.abs,
        .{ .order = options.order },
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
                validateContext(
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
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.abs_, ctx);
}

pub inline fn abs2(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Scalar(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.abs2,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.abs2_, ctx);
}

// Exponential functions
pub inline fn exp(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.exp_, ctx);
}

pub inline fn exp10(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp10,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.exp10_, ctx);
}

pub inline fn exp2(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp2,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.exp2_, ctx);
}

pub inline fn exp10m1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp10m1,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.exp10m1_, ctx);
}

pub inline fn exp2m1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.exp2m1,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.exp2m1_, ctx);
}

pub inline fn expm1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.expm1,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.expm1_, ctx);
}

pub inline fn log(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.log_, ctx);
}

pub inline fn log10(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log10,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.log10_, ctx);
}

pub inline fn log2(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log2,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.log2_, ctx);
}

pub inline fn log10p1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log10p1,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.log10p1_, ctx);
}

pub inline fn log2p1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log2p1,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.log2p1_, ctx);
}

pub inline fn log1p(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.log1p,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.log1p_, ctx);
}

// Power functions
pub inline fn pow(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.pow,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.pow_, ctx);
}

pub inline fn sqrt(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sqrt,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.sqrt_, ctx);
}

pub inline fn cbrt(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cbrt,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.cbrt_, ctx);
}

pub inline fn hypot(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.hypot,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.hypot_, ctx);
}

// Trigonometric functions
pub inline fn sin(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sin,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.sin_, ctx);
}

pub inline fn cos(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cos,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.cos_, ctx);
}

pub inline fn tan(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.tan,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.tan_, ctx);
}

pub inline fn asin(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.asin,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.asin_, ctx);
}

pub inline fn acos(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.acos,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.acos_, ctx);
}

pub inline fn atan(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.atan,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.atan_, ctx);
}

pub inline fn atan2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.atan2,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.atan2_, ctx);
}

pub inline fn sinpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sinpi,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.sinpi_, ctx);
}

pub inline fn cospi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cospi,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.cospi_, ctx);
}

pub inline fn tanpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.tanpi,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.tanpi_, ctx);
}

pub inline fn asinpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.asinpi,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.asinpi_, ctx);
}

pub inline fn acospi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.acospi,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.acospi_, ctx);
}

pub inline fn atanpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.atanpi,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.atanpi_, ctx);
}

pub inline fn atan2pi(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.atan2pi,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.atan2pi_, ctx);
}

// Hyperbolic functions
pub inline fn sinh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.sinh,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.sinh_, ctx);
}

pub inline fn cosh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.cosh,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.cosh_, ctx);
}

pub inline fn tanh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.tanh,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.tanh_, ctx);
}

pub inline fn asinh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.asinh,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.asinh_, ctx);
}

pub inline fn acosh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.acosh,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.acosh_, ctx);
}

pub inline fn atanh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.atanh,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.atanh_, ctx);
}

// Error and gamma functions
pub inline fn erf(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.erf,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.erf_, ctx);
}

pub inline fn erfc(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.erfc,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.erfc_, ctx);
}

pub inline fn gamma(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.gamma,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.gamma_, ctx);
}

pub inline fn lgamma(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.lgamma,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.lgamma_, ctx);
}

// Nearest integer operations
pub inline fn ceil(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !@TypeOf(x) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.ceil,
        .{ .order = options.order },
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

    comptime if (types.isArbitraryPrecision(O)) {
        if (types.isArbitraryPrecision(X)) {
            // Both types are arbitrary precision
            if (O == X) {
                // Equal types: output can be used for the operations, needing
                // only the output's allocator
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(X)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply1_(o, x, ops.ceil_, ctx);
}

pub inline fn add(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.add,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        }
    };

    return apply2_(o, x, y, ops.add_, ctx);
}

pub inline fn sub(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.sub,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        }
    };

    return apply2_(o, x, y, ops.sub_, ctx);
}

pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        if (types.numericType(C) == .int) {
            validateContext(
                @TypeOf(ctx),
                .{ .mode = .{ .type = int.Mode, .required = false } },
            );
        } else {
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2(
        allocator,
        x,
        y,
        ops.add,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
                validateContext(
                    @TypeOf(ctx),
                    .{
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                        .mode = .{ .type = int.Mode, .required = false },
                    },
                );
            } else {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            if (types.numericType(C) == .int) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .mode = .{ .type = int.Mode, .required = false } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            }
        }
    };

    return apply2_(o, x, y, ops.mul_, ctx);
}

pub inline fn div(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.div,
        .{ .order = options.order },
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
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                // Different types: internal allocator is required to perform
                // the operation at `x`'s precision, and then cast the result to
                // the output with the output's allocator
                validateContext(
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
            validateContext(
                @TypeOf(ctx),
                .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
            );
        }
    } else {
        if (types.isArbitraryPrecision(C)) {
            // Only the input is arbitrary precision, so we need the internal
            // allocator to perform the operation at `x`'s precision
            validateContext(
                @TypeOf(ctx),
                .{
                    .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        } else {
            // Both types are fixed precision, no special requirements
            validateContext(@TypeOf(ctx), .{});
        }
    };

    return apply2_(o, x, y, ops.div_, ctx);
}

pub inline fn eq(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.eq, .{ .order = options.order }, .{});
}

pub inline fn eq_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.eq_, ctx);
}

pub inline fn ne(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.ne, .{ .order = options.order }, .{});
}

pub inline fn ne_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.ne_, ctx);
}

pub inline fn lt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.lt, .{ .order = options.order }, .{});
}

pub inline fn lt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.lt_, ctx);
}

pub inline fn le(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.le, .{ .order = options.order }, .{});
}

pub inline fn le_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.le_, ctx);
}

pub inline fn gt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.gt, .{ .order = options.order }, .{});
}

pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.gt_, ctx);
}

pub inline fn ge(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.ge, .{ .order = options.order }, .{});
}

pub inline fn ge_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));

    comptime if (types.isArbitraryPrecision(O)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.ge_, ctx);
}

pub inline fn max(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.max,
        .{ .order = options.order },
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
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.max_, ctx);
}

pub inline fn min(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        order: ?array.Order = null,
    },
    ctx: anytype,
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2(
        allocator,
        x,
        y,
        ops.min,
        .{ .order = options.order },
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
        validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        validateContext(@TypeOf(ctx), .{});
    };

    return apply2_(o, x, y, ops.min_, ctx);
}
