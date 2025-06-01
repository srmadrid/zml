const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;
const CoerceToArray = types.CoerceToArray;
const Child = types.Child;
const needsAllocator = types.needsAllocator;

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");
const array = @import("array.zig");

pub inline fn add(
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.add not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.add not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.add(x, y, .{ .mode = options.mode }),
                .float => return float.add(x, y),
                .cfloat => return cfloat.add(x, y),
                .integer => {
                    if (options.allocator == null)
                        @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " requires an allocator");

                    // return integer.add(x, y, .{ .allocator = options.allocator, .mode = options.mode });
                    return .empty;
                },
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.add(x, y, .{ .mode = options.mode }),
                .float => return float.add(x, y),
                .cfloat => return cfloat.add(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.add(x, y),
                .cfloat => return cfloat.add(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.add(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn add_(
    out: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !void {
    comptime var O: type = @TypeOf(out);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.add_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y) or
        types.isArray(O) or types.isSlice(O))
        @compileError("zml.add_ not implemented for arrays or slices yet.");
    // return array.add(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.add_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.add_(out, x, y, .{ .mode = options.mode }),
                .float => return float.add_(out, x, y),
                .cfloat => return cfloat.add_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.add_(out, x, y, .{ .mode = options.mode }),
                .float => return float.add_(out, x, y),
                .cfloat => return cfloat.add_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.add_(out, x, y),
                .cfloat => return cfloat.add_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.add_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.add_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn sub(
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.sub not implemented for arrays or slices yet.");
    // return array.sub(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.sub not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.sub(x, y, .{ .mode = options.mode }),
                .float => return float.sub(x, y),
                .cfloat => return cfloat.sub(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.sub(x, y, .{ .mode = options.mode }),
                .float => return float.sub(x, y),
                .cfloat => return cfloat.sub(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.sub(x, y),
                .cfloat => return cfloat.sub(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.sub(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn sub_(
    out: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !void {
    comptime var O: type = @TypeOf(out);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.sub_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y) or
        types.isArray(O) or types.isSlice(O))
        @compileError("zml.sub_ not implemented for arrays or slices yet.");
    // return array.sub(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.sub_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.sub_(out, x, y, .{ .mode = options.mode }),
                .float => return float.sub_(out, x, y),
                .cfloat => return cfloat.sub_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.sub_(out, x, y, .{ .mode = options.mode }),
                .float => return float.sub_(out, x, y),
                .cfloat => return cfloat.sub_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.sub_(out, x, y),
                .cfloat => return cfloat.sub_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.sub_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.sub_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn mul(
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.mul not implemented for arrays or slices yet.");
    // return array.mul(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.mul not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.mul(x, y, .{ .mode = options.mode }),
                .float => return float.mul(x, y),
                .cfloat => return cfloat.mul(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.mul(x, y, .{ .mode = options.mode }),
                .float => return float.mul(x, y),
                .cfloat => return cfloat.mul(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.mul(x, y),
                .cfloat => return cfloat.mul(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.mul(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn mul_(
    out: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        mode: int.Mode = .default,
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !void {
    comptime var O: type = @TypeOf(out);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y) or
        types.isArray(O) or types.isSlice(O))
        @compileError("zml.mul_ not implemented for arrays or slices yet.");
    // return array.mul(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.mul_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.mul_(out, x, y, .{ .mode = options.mode }),
                .float => return float.mul_(out, x, y),
                .cfloat => return cfloat.mul_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.mul_(out, x, y, .{ .mode = options.mode }),
                .float => return float.mul_(out, x, y),
                .cfloat => return cfloat.mul_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.mul_(out, x, y),
                .cfloat => return cfloat.mul_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.mul_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.mul_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn div(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.div not implemented for arrays or slices yet.");
    // return array.div(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.div not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.div(x, y),
                .float => return float.div(x, y),
                .cfloat => return cfloat.div(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.div(x, y),
                .float => return float.div(x, y),
                .cfloat => return cfloat.div(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.div(x, y),
                .cfloat => return cfloat.div(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.div(x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn div_(
    out: anytype,
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !void {
    comptime var O: type = @TypeOf(out);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    _ = options.allocator;
    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.div_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y) or
        types.isArray(O) or types.isSlice(O))
        @compileError("zml.div_ not implemented for arrays or slices yet.");
    // return array.div(x, y, options.allocator);

    switch (types.numericType(X)) {
        .bool => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.div_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.div_(out, x, y),
                .float => return float.div_(out, x, y),
                .cfloat => return cfloat.div_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(Y)) {
                .bool, .int => return int.div_(out, x, y),
                .float => return float.div_(out, x, y),
                .cfloat => return cfloat.div_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool, .int, .float => return float.div_(out, x, y),
                .cfloat => return cfloat.div_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(Y)) {
                .bool, .int, .float, .cfloat => return cfloat.div_(out, x, y),
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(Y)) {
                .unsupported => unreachable,
                else => @compileError("zml.div_ between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn eq(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.zml.eq not implemented for arrays or slices yet.");
    // return array.eq(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.eq not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.eq(x, y),
                .float => return float.eq(x, y),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.eq(x, y),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.eq not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.eq between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn ne(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.ne not implemented for arrays or slices yet.");
    // return array.ne(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.ne(x, y),
                .float => return float.ne(x, y),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.ne(x, y),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ne not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ne between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn lt(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.lt not implemented for arrays or slices yet.");
    // return array.lt(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.lt not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.lt(x, y),
                .float => return float.lt(x, y),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.lt(x, y),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.lt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.lt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn le(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.le not implemented for arrays or slices yet.");
    // return array.le(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.le not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.le(x, y),
                .float => return float.le(x, y),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.le(x, y),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.le not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.le between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn gt(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.gt not implemented for arrays or slices yet.");
    // return array.gt(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.gt not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.gt(x, y),
                .float => return float.gt(x, y),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.gt(x, y),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.gt not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.gt between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn ge(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !CoerceToArray(Coerce(@TypeOf(x), @TypeOf(y)), bool) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.ge not implemented for arrays or slices yet.");
    // return array.ge(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ge not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int => return int.ge(x, y),
                .float => return float.ge(x, y),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return float.ge(x, y),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.ge not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.ge between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn max(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.max not implemented for arrays or slices yet.");
    // return array.max(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.max not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @max(x, y),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @max(x, y),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.max not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.max between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn min(
    x: anytype,
    y: anytype,
    options: struct {
        allocator: if (needsAllocator(Coerce(@TypeOf(x), @TypeOf(y)))) std.mem.Allocator else void = {},
    },
) !Coerce(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (types.isArray(X) or types.isSlice(X) or
        types.isArray(Y) or types.isSlice(Y))
        @compileError("zml.min not implemented for arrays or slices yet.");
    // return array.min(x, y, options.allocator);

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.min not defined for bools"),
        .int => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @min(x, y),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .int, .float => return @min(x, y),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .integer => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
        .expression => {
            switch (types.numericType(Y)) {
                .bool => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .cfloat => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .complex => @compileError("zml.min not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                .unsupported => unreachable,
                else => @compileError("zml.min between " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

// Include:
// - sub
// - mul
// - div
// - mod
// - pow
// - neg
// - abs
// - sqrt
// - sin
// - cos
// - tan
// - asin
// - acos
// - atan
// - sinh
// - cosh
// - tanh
// - asinh
// - acosh
// - atanh
// - sec
// - csc
// - cot
// - asec
// - acsc
// - acot
// - sech
// - csch
// - coth
// - asech
// - acsch
// - acoth
// - log
// - log10
// - log2
// - log1p
// - exp
// - any other basic math functions

pub inline fn abs(
    x: anytype,
    options: struct {
        allocator: if (needsAllocator(@TypeOf(x))) std.mem.Allocator else void = {},
    },
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (types.isArray(X) or types.isSlice(X))
        @compileError("zml.abs not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs not defined for " ++ @typeName(X)),
        .int => return int.abs,
        .float => return float.abs(x),
        .cfloat => return cfloat.abs(x),
        .unsupported => unreachable,
        else => @compileError("zml.abs not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn abs_(
    out: anytype,
    x: anytype,
    options: struct {
        allocator: if (needsAllocator(@TypeOf(x))) std.mem.Allocator else void = {},
    },
) !void {
    comptime var O: type = @TypeOf(out);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.abs_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X))
        @compileError("zml.abs_ not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.abs_ not defined for " ++ @typeName(X)),
        .int => out.* = try cast(O, int.abs(x), .{}),
        .float => out.* = try cast(O, float.abs(x), .{}),
        .cfloat => out.* = try cast(O, cfloat.abs(x), .{}),
        .unsupported => unreachable,
        else => @compileError("zml.abs_ not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn ceil(
    x: anytype,
    options: struct {
        allocator: if (needsAllocator(@TypeOf(x))) std.mem.Allocator else void = {},
    },
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (types.isArray(X) or types.isSlice(X))
        @compileError("zml.ceil not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        .int => return x,
        .float => return float.ceil(x),
        .cfloat => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        .complex => @compileError("zml.ceil not defined for " ++ @typeName(X)),
        .unsupported => unreachable,
        else => @compileError("zml.ceil not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn ceil_(
    out: anytype,
    x: anytype,
    options: struct {
        allocator: if (needsAllocator(@TypeOf(x))) std.mem.Allocator else void = {},
    },
) !void {
    comptime var O: type = @TypeOf(out);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ceil_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (types.isArray(O) or types.isSlice(O) or
        types.isArray(X) or types.isSlice(X))
        @compileError("zml.ceil_ not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool => @compileError("zml.ceil_ not defined for " ++ @typeName(X)),
        .int => out.* = try cast(O, x, .{}),
        .float => out.* = try cast(O, float.ceil(x), .{}),
        .cfloat => @compileError("zml.ceil_ not defined for " ++ @typeName(X)),
        .complex => @compileError("zml.ceil_ not defined for " ++ @typeName(X)),
        .unsupported => unreachable,
        else => @compileError("zml.ceil_ not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn copy(
    x: anytype,
    options: struct {
        allocator: if (needsAllocator(@TypeOf(x))) std.mem.Allocator else void = {},
    },
) @TypeOf(x) {
    const X: type = @TypeOf(x);

    comptime if (types.isArray(X) or types.isSlice(X))
        @compileError("zml.copy not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool, .int, .float, .cfloat => return x,
        .unsupported => unreachable,
        else => @compileError("zml.copy not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn deinit(
    x: anytype,
    options: struct {
        allocator: if (needsAllocator(Child(@TypeOf(x)))) std.mem.Allocator else void = {},
    },
) void {
    const X: type = Child(@TypeOf(x));

    comptime if (types.isArray(X) or types.isSlice(X))
        @compileError("zml.deinit not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(X)) {
        .bool, .int, .float, .cfloat => {},
        .unsupported => unreachable,
        else => @compileError("zml.deinit not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

test {}
