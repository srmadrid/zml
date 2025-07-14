const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const EnsureFloat = types.EnsureFloat;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const scast = types.scast;
const cast = types.cast;
const needsAllocator = types.needsAllocator;

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
        writeable: bool = true,
        // Add order
    },
) !Array(ReturnType1(op, Numeric(@TypeOf(x)))) {
    const X: type = @TypeOf(x);

    comptime if (!types.isArray(X) and !types.isSlice(X))
        @compileError("apply1: x must be an array or slice, got " ++ @typeName(X));

    switch (x.flags.storage) {
        .dense => return dense.apply1(allocator, Numeric(X), x, op, options.writeable),
        .strided => return strided.apply1(allocator, Numeric(X), x, op, options.writeable),
    }
}

pub fn apply1_(
    o: anytype,
    x: anytype,
    comptime op_: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("apply1_: o must be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1_: o must be an array or slice, got " ++ @typeName(O));

    if (o.flags.writeable == false)
        return error.ArrayNotWriteable;

    switch (o.flags.storage) {
        .dense => switch (x.flags.storage) {
            .dense => return dense.apply1_(Numeric(O), o, Numeric(X), x, op_, options.allocator),
            .strided => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, options.allocator),
        },
        .strided => switch (x.flags.storage) {
            .dense => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, options.allocator),
            .strided => return strided.apply1_(Numeric(O), o, Numeric(X), x, op_, options.allocator),
        },
    }
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    options: struct {
        order: ?array.Order = null,
        writeable: bool = true,
    },
) !Array(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isArray(X) and !types.isSlice(X) and
        !types.isArray(Y) and !types.isSlice(Y))
        @compileError("apply2: at least one of x or y must be an array or slice, got " ++
            @typeName(X) ++ " and " ++ @typeName(Y));

    if (comptime !types.isArray(X) and !types.isSlice(X)) {
        // x is a scalar, only consider y's storage
        switch (y.flags.storage) {
            .dense => return dense.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
            .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
        }
    } else if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
        // y is a scalar, only consider x's storage
        switch (x.flags.storage) {
            .dense => return dense.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
            .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
        }
    } else {
        switch (x.flags.storage) {
            .dense => switch (y.flags.storage) {
                .dense => return dense.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
                .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
            },
            .strided => switch (y.flags.storage) {
                .dense => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
                .strided => return strided.apply2(allocator, Numeric(X), x, Numeric(Y), y, op, .{ .order = options.order, .writeable = options.writeable }),
            },
        }
    }
}

pub fn apply2_(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1: o must be an array or slice, got " ++ @typeName(O));

    if (comptime !types.isArray(X) and !types.isSlice(X)) {
        if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
            // x and y are scalars, only consider o's storage
            switch (o.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
            }
        }

        // x is a scalar, only consider o and y's storage
        switch (o.flags.storage) {
            .dense => switch (y.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
            },
            .strided => switch (y.flags.storage) {
                .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
            },
        }
    } else if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
        // y is a scalar, only consider o and x's storage
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
            },
            .strided => switch (x.flags.storage) {
                .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
            },
        }
    } else {
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => switch (y.flags.storage) {
                    .dense => return dense.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                },
                .strided => switch (y.flags.storage) {
                    .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                },
            },
            .strided => switch (x.flags.storage) {
                .dense => switch (y.flags.storage) {
                    .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                },
                .strided => switch (y.flags.storage) {
                    .dense => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
                    .strided => return strided.apply2_(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_, options.allocator),
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
        writeable: bool = true,
    },
) !Array(Scalar(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.abs, .{ .writeable = options.writeable });
}

pub inline fn abs_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.abs_, .{ .allocator = options.allocator });
}

// Exponential functions
pub inline fn exp(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.exp, .{ .writeable = options.writeable });
}

pub inline fn exp_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.exp_, .{ .allocator = options.allocator });
}

pub inline fn exp10(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.exp10, .{ .writeable = options.writeable });
}

pub inline fn exp10_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.exp10_, .{ .allocator = options.allocator });
}

pub inline fn exp2(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.exp2, .{ .writeable = options.writeable });
}

pub inline fn exp2_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.exp2_, .{ .allocator = options.allocator });
}

pub inline fn exp10m1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.exp10m1, .{ .writeable = options.writeable });
}

pub inline fn exp10m1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.exp10m1_, .{ .allocator = options.allocator });
}

pub inline fn exp2m1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.exp2m1, .{ .writeable = options.writeable });
}

pub inline fn exp2m1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.exp2m1_, .{ .allocator = options.allocator });
}

pub inline fn expm1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.expm1, .{ .writeable = options.writeable });
}

pub inline fn expm1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.expm1_, .{ .allocator = options.allocator });
}

pub inline fn log(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.log, .{ .writeable = options.writeable });
}

pub inline fn log_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.log_, .{ .allocator = options.allocator });
}

pub inline fn log10(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.log10, .{ .writeable = options.writeable });
}

pub inline fn log10_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.log10_, .{ .allocator = options.allocator });
}

pub inline fn log2(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.log2, .{ .writeable = options.writeable });
}

pub inline fn log2_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.log2_, .{ .allocator = options.allocator });
}

pub inline fn log10p1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.log10p1, .{ .writeable = options.writeable });
}

pub inline fn log10p1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.log10p1_, .{ .allocator = options.allocator });
}

pub inline fn log2p1(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.log2p1, .{ .writeable = options.writeable });
}

pub inline fn log2p1_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.log2p1_, .{ .allocator = options.allocator });
}

pub inline fn log1p(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.log1p, .{ .writeable = options.writeable });
}

pub inline fn log1p_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.log1p_, .{ .allocator = options.allocator });
}

// Power functions
pub inline fn pow(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply2(allocator, x, y, ops.pow, .{ .writeable = options.writeable });
}

pub inline fn pow_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.pow_, .{ .allocator = options.allocator });
}

pub inline fn sqrt(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.sqrt, .{ .writeable = options.writeable });
}

pub inline fn sqrt_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.sqrt_, .{ .allocator = options.allocator });
}

pub inline fn cbrt(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.cbrt, .{ .writeable = options.writeable });
}

pub inline fn cbrt_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.cbrt_, .{ .allocator = options.allocator });
}

pub inline fn hypot(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply2(allocator, x, y, ops.hypot, .{ .writeable = options.writeable });
}

pub inline fn hypot_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.hypot_, .{ .allocator = options.allocator });
}

// Trigonometric functions
pub inline fn sin(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.sin, .{ .writeable = options.writeable });
}

pub inline fn sin_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.sin_, .{ .allocator = options.allocator });
}

pub inline fn cos(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.cos, .{ .writeable = options.writeable });
}

pub inline fn cos_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.cos_, .{ .allocator = options.allocator });
}

pub inline fn tan(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.tan, .{ .writeable = options.writeable });
}

pub inline fn tan_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.tan_, .{ .allocator = options.allocator });
}

pub inline fn asin(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.asin, .{ .writeable = options.writeable });
}

pub inline fn asin_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.asin_, .{ .allocator = options.allocator });
}

pub inline fn acos(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.acos, .{ .writeable = options.writeable });
}

pub inline fn acos_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.acos_, .{ .allocator = options.allocator });
}

pub inline fn atan(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.atan, .{ .writeable = options.writeable });
}

pub inline fn atan_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.atan_, .{ .allocator = options.allocator });
}

pub inline fn atan2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply2(allocator, x, y, ops.atan2, .{ .writeable = options.writeable });
}

pub inline fn atan2_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.atan2_, .{ .allocator = options.allocator });
}

pub inline fn sinpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.sinpi, .{ .writeable = options.writeable });
}

pub inline fn sinpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.sinpi_, .{ .allocator = options.allocator });
}

pub inline fn cospi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.cospi, .{ .writeable = options.writeable });
}

pub inline fn cospi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.cospi_, .{ .allocator = options.allocator });
}

pub inline fn tanpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.tanpi, .{ .writeable = options.writeable });
}

pub inline fn tanpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.tanpi_, .{ .allocator = options.allocator });
}

pub inline fn asinpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.asinpi, .{ .writeable = options.writeable });
}

pub inline fn asinpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.asinpi_, .{ .allocator = options.allocator });
}

pub inline fn acospi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.acospi, .{ .writeable = options.writeable });
}

pub inline fn acospi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.acospi_, .{ .allocator = options.allocator });
}

pub inline fn atanpi(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.atanpi, .{ .writeable = options.writeable });
}

pub inline fn atanpi_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.atanpi_, .{ .allocator = options.allocator });
}

pub inline fn atan2pi(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply2(allocator, x, y, ops.atan2pi, .{ .writeable = options.writeable });
}

pub inline fn atan2pi_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.atan2pi_, .{ .allocator = options.allocator });
}

// Hyperbolic functions
pub inline fn sinh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.sinh, .{ .writeable = options.writeable });
}

pub inline fn sinh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.sinh_, .{ .allocator = options.allocator });
}

pub inline fn cosh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.cosh, .{ .writeable = options.writeable });
}

pub inline fn cosh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.cosh_, .{ .allocator = options.allocator });
}

pub inline fn tanh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.tanh, .{ .writeable = options.writeable });
}

pub inline fn tanh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.tanh_, .{ .allocator = options.allocator });
}

pub inline fn asinh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.asinh, .{ .writeable = options.writeable });
}

pub inline fn asinh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.asinh_, .{ .allocator = options.allocator });
}

pub inline fn acosh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.acosh, .{ .writeable = options.writeable });
}

pub inline fn acosh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.acosh_, .{ .allocator = options.allocator });
}

pub inline fn atanh(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.atanh, .{ .writeable = options.writeable });
}

pub inline fn atanh_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.atanh_, .{ .allocator = options.allocator });
}

// Error and gamma functions
pub inline fn erf(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.erf, .{ .writeable = options.writeable });
}

pub inline fn erf_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.erf_, .{ .allocator = options.allocator });
}

pub inline fn erfc(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.erfc, .{ .writeable = options.writeable });
}

pub inline fn erfc_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.erfc_, .{ .allocator = options.allocator });
}

pub inline fn gamma(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.gamma, .{ .writeable = options.writeable });
}

pub inline fn gamma_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.gamma_, .{ .allocator = options.allocator });
}

pub inline fn lgamma(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(EnsureFloat(Numeric(@TypeOf(x)))) {
    return apply1(allocator, x, ops.lgamma, .{ .writeable = options.writeable });
}

pub inline fn lgamma_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.lgamma_, .{ .allocator = options.allocator });
}

// Nearest integer operations
pub inline fn ceil(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !@TypeOf(x) {
    return apply1(allocator, x, ops.ceil, .{ .writeable = options.writeable });
}

pub inline fn ceil_(
    o: anytype,
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, x, ops.ceil_, .{ .allocator = options.allocator });
}

inline fn addDefault(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.add(x, y, .{
        .mode = int.Mode.default,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

inline fn addWrap(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.add(x, y, .{
        .mode = int.Mode.wrap,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

inline fn addSaturate(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.add(x, y, .{
        .mode = int.Mode.saturate,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

pub inline fn add(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
        mode: int.Mode = .default,
    },
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    switch (options.mode) {
        .default => return apply2(allocator, x, y, addDefault, .{ .writeable = options.writeable }),
        .wrap => return apply2(allocator, x, y, addWrap, .{ .writeable = options.writeable }),
        .saturate => return apply2(allocator, x, y, addSaturate, .{ .writeable = options.writeable }),
    }
}

inline fn addDefault_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.add_(o, x, y, .{
        .mode = int.Mode.default,
        .allocator = options.allocator,
    });
}

inline fn addWrap_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.add_(o, x, y, .{
        .mode = int.Mode.wrap,
        .allocator = options.allocator,
    });
}

inline fn addSaturate_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.add_(o, x, y, .{
        .mode = int.Mode.saturate,
        .allocator = options.allocator,
    });
}

pub inline fn add_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        mode: int.Mode = .default,
    },
) !void {
    switch (options.mode) {
        .default => return apply2_(o, x, y, addDefault_, .{ .allocator = options.allocator }),
        .wrap => return apply2_(o, x, y, addWrap_, .{ .allocator = options.allocator }),
        .saturate => return apply2_(o, x, y, addSaturate_, .{ .allocator = options.allocator }),
    }
}

inline fn subDefault(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.sub(x, y, .{
        .mode = int.Mode.default,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

inline fn subWrap(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.sub(x, y, .{
        .mode = int.Mode.wrap,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

inline fn subSaturate(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.sub(x, y, .{
        .mode = int.Mode.saturate,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

pub inline fn sub(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
        mode: int.Mode = .default,
    },
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    switch (options.mode) {
        .default => return apply2(allocator, x, y, subDefault, .{ .writeable = options.writeable }),
        .wrap => return apply2(allocator, x, y, subWrap, .{ .writeable = options.writeable }),
        .saturate => return apply2(allocator, x, y, subSaturate, .{ .writeable = options.writeable }),
    }
}

inline fn subDefault_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.sub_(o, x, y, .{
        .mode = int.Mode.default,
        .allocator = options.allocator,
    });
}

inline fn subWrap_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.sub_(o, x, y, .{
        .mode = int.Mode.wrap,
        .allocator = options.allocator,
    });
}

inline fn subSaturate_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.sub_(o, x, y, .{
        .mode = int.Mode.saturate,
        .allocator = options.allocator,
    });
}

pub inline fn sub_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        mode: int.Mode = .default,
    },
) !void {
    switch (options.mode) {
        .default => return apply2_(o, x, y, subDefault_, .{ .allocator = options.allocator }),
        .wrap => return apply2_(o, x, y, subWrap_, .{ .allocator = options.allocator }),
        .saturate => return apply2_(o, x, y, subSaturate_, .{ .allocator = options.allocator }),
    }
}

inline fn mulDefault(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.mul(x, y, .{
        .mode = int.Mode.default,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

inline fn mulWrap(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.mul(x, y, .{
        .mode = int.Mode.wrap,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

inline fn mulSaturate(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        writeable: bool = true,
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    return ops.mul(x, y, .{
        .mode = int.Mode.saturate,
        .allocator = options.allocator,
        .writeable = options.writeable,
    });
}

pub inline fn mul(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
        mode: int.Mode = .default,
    },
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    switch (options.mode) {
        .default => return apply2(allocator, x, y, mulDefault, .{ .writeable = options.writeable }),
        .wrap => return apply2(allocator, x, y, mulWrap, .{ .writeable = options.writeable }),
        .saturate => return apply2(allocator, x, y, mulSaturate, .{ .writeable = options.writeable }),
    }
}

inline fn mulDefault_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.mul_(o, x, y, .{
        .mode = int.Mode.default,
        .allocator = options.allocator,
    });
}

inline fn mulWrap_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.mul_(o, x, y, .{
        .mode = int.Mode.wrap,
        .allocator = options.allocator,
    });
}

inline fn mulSaturate_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return ops.mul_(o, x, y, .{
        .mode = int.Mode.saturate,
        .allocator = options.allocator,
    });
}

pub inline fn mul_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
        mode: int.Mode = .default,
    },
) !void {
    switch (options.mode) {
        .default => return apply2_(o, x, y, mulDefault_, .{ .allocator = options.allocator }),
        .wrap => return apply2_(o, x, y, mulWrap_, .{ .allocator = options.allocator }),
        .saturate => return apply2_(o, x, y, mulSaturate_, .{ .allocator = options.allocator }),
    }
}

pub inline fn div(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    return apply2(allocator, x, y, ops.div, .{ .writeable = options.writeable });
}

pub inline fn div_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.div_, .{ .allocator = options.allocator });
}

pub inline fn eq(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.eq, .{ .writeable = options.writeable });
}

pub inline fn eq_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.eq_, .{ .allocator = options.allocator });
}

pub inline fn ne(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.ne, .{ .writeable = options.writeable });
}

pub inline fn ne_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.ne_, .{ .allocator = options.allocator });
}

pub inline fn lt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.lt, .{ .writeable = options.writeable });
}

pub inline fn lt_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.lt_, .{ .allocator = options.allocator });
}

pub inline fn le(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.le, .{ .writeable = options.writeable });
}

pub inline fn le_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.le_, .{ .allocator = options.allocator });
}

pub inline fn gt(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.gt, .{ .writeable = options.writeable });
}

pub inline fn gt_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.gt_, .{ .allocator = options.allocator });
}

pub inline fn ge(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(bool) {
    return apply2(allocator, x, y, ops.ge, .{ .writeable = options.writeable });
}

pub inline fn ge_(
    o: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply2_(o, x, y, ops.ge_, .{ .allocator = options.allocator });
}
