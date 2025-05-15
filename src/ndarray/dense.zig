const ndarray = @import("../ndarray.zig");
const NDArray = ndarray.NDArray;

pub inline fn index(comptime T: type, array: *const NDArray(T), position: []const usize) usize {
    var idx: usize = 0;
    for (0..array.ndim) |i| {
        idx += position[i] * array.metadata.dense.strides[i];
    }

    return idx;
}

pub inline fn checkPosition(comptime T: type, array: *const NDArray(T), position: []const usize) !void {
    if (position.len > array.ndim) {
        return ndarray.Error.DimensionMismatch;
    }

    for (0..position.len) |i| {
        if (position[i] >= array.shape[i]) {
            return ndarray.Error.PositionOutOfBounds;
        }
    }
}
