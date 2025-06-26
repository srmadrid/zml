const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
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

// Only for in-pleace operations
pub fn apply1_(
    o: anytype,
    comptime op_: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1: o must be an array or slice, got " ++ @typeName(O));

    if (o.flags.writeable == false)
        return error.ArrayNotWriteable;

    switch (o.flags.storage) {
        .dense => return dense.apply1_(Numeric(O), o, op_, options.allocator),
        .strided => return strided.apply1_(Numeric(O), o, op_, options.allocator),
    }
}

pub fn apply1_to(
    o: anytype,
    x: anytype,
    comptime op_to: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("apply1_to: o must be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1_to: o must be an array or slice, got " ++ @typeName(O));

    if (o.flags.writeable == false)
        return error.ArrayNotWriteable;

    switch (o.flags.storage) {
        .dense => switch (x.flags.storage) {
            .dense => return dense.apply1_to(Numeric(O), o, Numeric(X), x, op_to, options.allocator),
            .strided => return strided.apply1_to(Numeric(O), o, Numeric(X), x, op_to, options.allocator),
        },
        .strided => switch (x.flags.storage) {
            .dense => return strided.apply1_to(Numeric(O), o, Numeric(X), x, op_to, options.allocator),
            .strided => return strided.apply1_to(Numeric(O), o, Numeric(X), x, op_to, options.allocator),
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
    y: anytype,
    comptime op_: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    comptime var O: type = @TypeOf(o);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O) and !types.isSlice(O))
        @compileError("apply1: o must be an array or slice, got " ++ @typeName(O));

    if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
        // y is a scalar, only consider o's storage
        switch (o.flags.storage) {
            .dense => return dense.apply2_(Numeric(O), o, Numeric(Y), y, op_, options.allocator),
            .strided => return strided.apply2_(Numeric(O), o, Numeric(Y), y, op_, options.allocator),
        }
    } else {
        switch (o.flags.storage) {
            .dense => switch (y.flags.storage) {
                .dense => return dense.apply2_(Numeric(O), o, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(Y), y, op_, options.allocator),
            },
            .strided => switch (y.flags.storage) {
                .dense => return strided.apply2_(Numeric(O), o, Numeric(Y), y, op_, options.allocator),
                .strided => return strided.apply2_(Numeric(O), o, Numeric(Y), y, op_, options.allocator),
            },
        }
    }
}

pub fn apply2_to(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime op_to: anytype,
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
                .dense => return dense.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
            }
        }

        // x is a scalar, only consider o and y's storage
        switch (o.flags.storage) {
            .dense => switch (y.flags.storage) {
                .dense => return dense.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
            },
            .strided => switch (y.flags.storage) {
                .dense => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
            },
        }
    } else if (comptime !types.isArray(Y) and !types.isSlice(Y)) {
        // y is a scalar, only consider o and x's storage
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => return dense.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
            },
            .strided => switch (x.flags.storage) {
                .dense => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
            },
        }
    } else {
        switch (o.flags.storage) {
            .dense => switch (x.flags.storage) {
                .dense => switch (y.flags.storage) {
                    .dense => return dense.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                    .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                },
                .strided => switch (y.flags.storage) {
                    .dense => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                    .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                },
            },
            .strided => switch (x.flags.storage) {
                .dense => switch (y.flags.storage) {
                    .dense => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                    .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                },
                .strided => switch (y.flags.storage) {
                    .dense => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                    .strided => return strided.apply2_to(Numeric(O), o, Numeric(X), x, Numeric(Y), y, op_to, options.allocator),
                },
            },
        }
    }
}

pub inline fn abs(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(Scalar(Scalar(@TypeOf(x)))) {
    return apply1(allocator, x, ops.abs, .{ .writeable = options.writeable });
}

pub inline fn abs_(
    o: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, ops.abs_, .{ .allocator = options.allocator });
}

pub inline fn ceil(
    allocator: std.mem.Allocator,
    x: anytype,
    options: struct {
        writeable: bool = true,
    },
) !Array(Scalar(@TypeOf(x))) {
    return apply1(allocator, x, ops.ceil, .{ .writeable = options.writeable });
}

pub inline fn ceil_(
    o: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) !void {
    return apply1_(o, ops.ceil_, .{ .allocator = options.allocator });
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
    return apply2(allocator, x, y, ops.add, .{ .writeable = options.writeable });
}
