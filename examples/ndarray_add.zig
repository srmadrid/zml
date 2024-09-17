const std = @import("std");
const cml = @import("camel");

pub fn main() !void {
    const a: std.mem.Allocator = std.heap.page_allocator;

    var B: cml.NDArray(f64) = try cml.NDArray(f64).initFlags(a, &[_]usize{ 1, 4 }, cml.ndarray.Flags{ .RowMajorContiguous = false, .ColumnMajorContiguous = true });
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

    var C: cml.NDArray(f64) = try cml.NDArray(f64).initFlags(a, &[_]usize{ 3, 1 }, cml.ndarray.Flags{ .RowMajorContiguous = true, .ColumnMajorContiguous = false });
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

    var D: cml.NDArray(f64) = try cml.NDArray(f64).initFlags(a, &[_]usize{ 3, 4 }, cml.ndarray.Flags{ .RowMajorContiguous = false, .ColumnMajorContiguous = true });
    defer D.deinit();
    try D.mult(B, C);
    // try cml.NDArray(u64).add(&D, B, C);
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{D.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n\n", .{});
    }
}
