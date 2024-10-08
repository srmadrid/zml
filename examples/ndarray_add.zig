const std = @import("std");
const zml = @import("zml");

pub fn main() !void {
    const a: std.mem.Allocator = std.heap.page_allocator;

    var B: zml.NDArray(f64) = try zml.NDArray(f64).initFlags(a, &[_]usize{ 1, 4 }, zml.ndarray.Flags{ .RowMajorContiguous = false, .ColumnMajorContiguous = true });
    defer B.deinit();
    for (0..B.size) |i| {
        B.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("B =\n", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{B.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n", .{});
    }

    var C: zml.NDArray(f64) = try zml.NDArray(f64).initFlags(a, &[_]usize{ 3, 1 }, zml.ndarray.Flags{ .RowMajorContiguous = true, .ColumnMajorContiguous = false });
    defer C.deinit();
    C.setAll(10);
    std.debug.print("\nC =\n", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{C.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n", .{});
    }

    var D: zml.NDArray(f64) = try zml.NDArray(f64).initFlags(a, &[_]usize{ 3, 4 }, zml.ndarray.Flags{ .RowMajorContiguous = false, .ColumnMajorContiguous = true });
    defer D.deinit();
    try D.mult(B, C);
    // try zml.NDArray(u64).add(&D, B, C);
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{D.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n\n", .{});
    }
}
