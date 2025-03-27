const std = @import("std");
const types = @import("types.zig");

/// Adds two elements of any supported numeric type and stores the result in the
/// output variable. `left` and `right` must be coercible by `Coerce`, and `out`
/// must be a pointer of type `Coerce(@TypeOf(left), @TypeOf(right))`. The
/// optional allocator is needed only if unmanaged types are used.
pub inline fn add(
    out: anytype,
    left: anytype,
    right: anytype,
    options: struct {
        allocator_out: ?std.mem.Allocator = null,
        allocator_left: ?std.mem.Allocator = null,
        allocator_right: ?std.mem.Allocator = null,
    },
) void {
    // This implementations seems like the way to go?
    const T = @TypeOf(out);
    const K = @TypeOf(left);
    const V = @TypeOf(right);

    if (!types.canCoerce(K, T) or !types.canCoerce(V, T))
        if (types.Coerce(types.Coerce(K, V), T) != T)
            @compileError("Cannot coerce " ++ @typeName(K) ++ " and " ++ @typeName(V) ++ " to " ++ @typeName(T));

    // Casting may be expensive especially to arbitrary precision, invoke instead
    // builtin functions to add fixed to arbitrary precision types without
    // converting to arbitrary precision first. For now this is fine since
    // arbitrary precision types are not supported yet.
    const o = out;
    const l = types.cast(T, left, .{ .allocator = options.allocator_left });
    const r = types.cast(T, right, .{ .allocator = options.allocator_right });

    switch (types.numericType(T)) {
        .bool => @compileError("Cannot add boolean types"),
        .int => o.* = l + r,
        .float => o.* = l + r,
        .cfloat => o.* = T.add(l, r),
        .integer => @compileError("Not implemented yet"),
        .rational => @compileError("Not implemented yet"),
        .real => @compileError("Not implemented yet"),
        .complex => @compileError("Not implemented yet"),
        .expression => @compileError("Not implemented yet"),
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
