const std = @import("std");
const NDArray = @import("../../ndarray.zig").NDArray;

pub fn getrf(comptime T: type, a: *NDArray(T), ipiv: *NDArray(usize)) !void {
    return;
}
