const std = @import("std");
const types = @import("types.zig");

/// Adds two elements of any supported numeric type and stores the result in the
/// output variable. `left` and `right` must be coercible by `Coerce`, and `out`
/// must be a pointer of type `Coerce(@TypeOf(left), @TypeOf(right))`. The
/// optional allocator is needed only if unmanaged types are used.
pub inline fn add(allocator_out: anytype, out: anytype, allocator_left: anytype, left: anytype, allocator_right: anytype, right: anytype) void {
    if (@TypeOf(allocator_out) != ?std.mem.Allocator or @TypeOf(allocator_left) != ?std.mem.Allocator or @TypeOf(allocator_right) != ?std.mem.Allocator)
        @compileError("Allocators must be optional allocators.");

    const T = @TypeOf(out);
    const K = @TypeOf(left);
    const V = @TypeOf(right);

    const tnumeric = types.numericType(T);
    const knumeric = types.numericType(K);
    const vnumeric = types.numericType(V);
    _ = tnumeric;
    _ = knumeric;
    _ = vnumeric;

    if (T != types.Coerce(K, V))
        @compileError("Output type must be: Coerce(@TypeOf(left), @TypeOf(right)).");
}
