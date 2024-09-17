const std = @import("std");
const NDArray = @import("../../ndarray.zig").NDArray;
const NDArrayError = @import("../../ndarray.zig").NDArrayError;

pub fn getrf(comptime T: type, a: *NDArray(T), ipiv: *NDArray(usize)) !void {
    return;
}
