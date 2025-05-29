const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;
const Scalar = types.Scalar;
const Coerce = types.Coerce;
const CoerceToArray = types.CoerceToArray;

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");
const array = @import("array.zig");

pub inline fn add(
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) !Coerce(@TypeOf(left), @TypeOf(right)) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);

    if (types.isArray(L) or types.isSlice(L) or
        types.isArray(R) or types.isSlice(R))
        @compileError("zml.add not implemented for arrays or slices yet.");
    // return array.add(left, right, options.allocator);

    switch (types.numericType(L)) {
        .bool => {
            switch (types.numericType(R)) {
                .bool => @compileError("add not defined for two bools"),
                .int => return int.add(left, right, .{ .mode = options.mode }),
                .float => return float.add(left, right),
                .cfloat => return cfloat.add(left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(R)) {
                .bool, .int => return int.add(left, right, .{ .mode = options.mode }),
                .float => return float.add(left, right),
                .cfloat => return cfloat.add(left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(R)) {
                .bool, .int, .float => return float.add(left, right),
                .cfloat => return cfloat.add(left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .bool, .int, .float, .cfloat => return cfloat.add(left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn add_(
    out: anytype,
    left: anytype,
    right: anytype,
    options: struct {
        comptime mode: int.Mode = .default,
        allocator: ?std.mem.Allocator = null,
    },
) void {
    const O: type = types.Child(@TypeOf(out));
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.add_ requires the output to be a pointer to a mutable type, got " ++ @typeName(O));

    if (types.isArray(L) or types.isSlice(L) or
        types.isArray(R) or types.isSlice(R))
        @compileError("zml.add not implemented for arrays or slices yet.");
    // return array.add(left, right, options.allocator);

    switch (types.numericType(L)) {
        .bool => {
            switch (types.numericType(R)) {
                .bool => @compileError("add not defined for two bools"),
                .int => return int.add_(out, left, right, .{ .mode = options.mode }),
                .float => return float.add_(out, left, right),
                .cfloat => return cfloat.add_(out, left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .int => {
            switch (types.numericType(R)) {
                .bool, .int => return int.add_(out, left, right, .{ .mode = options.mode }),
                .float => return float.add_(out, left, right),
                .cfloat => return cfloat.add_(out, left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(R)) {
                .bool, .int, .float => return float.add_(out, left, right),
                .cfloat => return cfloat.add_(out, left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .cfloat => {
            switch (types.numericType(R)) {
                .bool, .int, .float, .cfloat => return cfloat.add_(out, left, right),
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .integer => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .complex => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .expression => {
            switch (types.numericType(R)) {
                .unsupported => unreachable,
                else => @compileError("add between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .unsupported => unreachable,
    }
}

pub inline fn eq(
    left: anytype,
    right: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) CoerceToArray(Coerce(@TypeOf(left), @TypeOf(right)), bool) {
    const L: type = @TypeOf(left);
    const R: type = @TypeOf(right);

    if (types.isArray(L) or types.isSlice(L) or
        types.isArray(R) or types.isSlice(R))
        @compileError("zml.cmp not implemented for arrays or slices yet.");
    // return array.cmp(left, right, options.allocator);

    _ = options.allocator;
    switch (types.numericType(L)) {
        .bool => @compileError("cmp not defined for bools"),
        .int => {
            switch (types.numericType(R)) {
                .bool => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .int => return int.eq(left, right),
                .float => return float.eq(left, right),
                .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .unsupported => unreachable,
                else => @compileError("cmp between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .float => {
            switch (types.numericType(R)) {
                .bool => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .int, .float => return float.eq(left, right),
                .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .unsupported => unreachable,
                else => @compileError("cmp between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
        .integer => {
            switch (types.numericType(R)) {
                .bool => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .unsupported => unreachable,
                else => @compileError("cmp between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .rational => {
            switch (types.numericType(R)) {
                .bool => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .unsupported => unreachable,
                else => @compileError("cmp between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .real => {
            switch (types.numericType(R)) {
                .bool => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .unsupported => unreachable,
                else => @compileError("cmp between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
            }
        },
        .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
        .expression => {
            switch (types.numericType(R)) {
                .bool => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .cfloat => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .complex => @compileError("cmp not defined for " ++ @typeName(L) ++ " and " ++ @typeName(R)),
                .unsupported => unreachable,
                else => @compileError("cmp between " ++ @typeName(L) ++ " and " ++ @typeName(R) ++ " not implemented yet"),
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

pub fn abs(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator,
    },
) @TypeOf(x) {
    _ = options.allocator;
    switch (types.numericType(@TypeOf(x))) {
        .bool => @compileError("abs does not support bool arguments"),
        .int => return int.abs,
        .float => return float.abs(x),
        .cfloat => return cfloat.abs(x),
        .integer => @compileError("abs not implemented for integers yet"),
        .rational => @compileError("abs not implemented for rationals yet"),
        .real => @compileError("abs not implemented for reals yet"),
        .complex => @compileError("abs not implemented for complexes yet"),
        .expression => @compileError("abs not implemented for expressions yet"),
        .unsupported => unreachable,
    }
}

pub inline fn deinit(
    x: anytype,
    options: struct {
        allocator: ?std.mem.Allocator = null,
    },
) void {
    if (types.isArray(@TypeOf(x)) or types.isSlice(@TypeOf(x)))
        @compileError("zml.deinit not implemented for arrays or slices yet.");

    _ = options.allocator;
    switch (types.numericType(@TypeOf(x))) {
        .bool, .int, .float, .cfloat => {},
        .integer => @compileError("deinit not implemented for integers yet"),
        .rational => @compileError("deinit not implemented for rationals yet"),
        .real => @compileError("deinit not implemented for reals yet"),
        .complex => @compileError("deinit not implemented for complexes yet"),
        .expression => @compileError("deinit not implemented for expressions yet"),
        .unsupported => unreachable,
    }
}

test {}
