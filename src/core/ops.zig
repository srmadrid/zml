const std = @import("std");
const types = @import("types.zig");

/// Adds two elements of any supported numeric type and stores the result in the
/// output variable. `left` and `right` must be coercible by `Coerce`, and `out`
/// must be a pointer of type `Coerce(@TypeOf(left), @TypeOf(right))`. The
/// optional allocator is needed only if unmanaged types are used.
pub inline fn add(allocator_out: anytype, out: anytype, allocator_left: anytype, left: anytype, allocator_right: anytype, right: anytype) void {
    if (@TypeOf(allocator_out) != ?std.mem.Allocator or @TypeOf(allocator_left) != ?std.mem.Allocator or @TypeOf(allocator_right) != ?std.mem.Allocator) @compileError("Allocators must be optional allocators.");

    const T = @TypeOf(out);
    const K = @TypeOf(left);
    const V = @TypeOf(right);

    if (!types.isScalar(T) or !types.isScalar(K) or !types.isScalar(V)) @compileError("Only scalar types are supported.");

    if (T != types.Coerce(K, V)) @compileError("Cannot coerce " ++ @typeName(K) ++ " and " ++ @typeName(V) ++ " to " ++ @typeName(T) ++ ", instead use " ++ @typeName(types.Coerce(K, V)));
}

/// Subtracts two elements of any supported type and stores the result in the
/// output variable. `left` and `right` must be of the same type, and `out` must
/// be a pointer of that same type.
//pub inline fn _sub(out: anytype, left: anytype, right: anytype) void {}

/// Multiplies two elements of any supported type and stores the result in the
/// output variable. `left` and `right` must be of the same type, and `out` must
/// be a pointer of that same type.
//pub inline fn _mul(out: anytype, left: anytype, right: anytype) void {}

/// Divides two elements of any supported type and stores the result in the
/// output variable. `left` and `right` must be of the same type, and `out` must
/// be a pointer of that same type.
//pub inline fn _div(out: anytype, left: anytype, right: anytype) void {}
fn a() void {
    return;
}
