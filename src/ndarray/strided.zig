const types = @import("../types.zig");
const cast = types.cast;
const int = @import("../int.zig");
const ndarray = @import("../ndarray.zig");
const NDArray = ndarray.NDArray;

pub inline fn index(comptime T: type, array: *const NDArray(T), position: []const usize) usize {
    var idx: usize = array.metadata.strided.offset;
    for (0..array.ndim) |i| {
        const stride: isize = array.metadata.strided.strides[i];
        if (stride < 0) {
            idx -= position[i] * cast(usize, int.abs(stride), .{});
        } else {
            idx += position[i] * cast(usize, stride, .{});
        }
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
