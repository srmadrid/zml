const int = @import("../int.zig");

const numeric = @import("../numeric.zig");

const array = @import("../array.zig");

pub fn broadcastShapes(shapes: []const []const usize) !struct { ndim: usize, shape: [array.max_dimensions]usize } {
    if (shapes.len == 0)
        return array.Error.ZeroDimension;

    var ndim: usize = 0;
    for (shapes) |shape| {
        if (shape.len == 0)
            return array.Error.ZeroDimension;

        if (shape.len > array.max_dimensions)
            return array.Error.TooManyDimensions;

        ndim = int.max(ndim, shape.len);
    }

    var result: [array.max_dimensions]usize = undefined;

    var out_dim: usize = ndim;
    while (out_dim > 0) {
        out_dim -= 1;

        var max_dim: usize = 1;

        for (shapes) |shape| {
            const trailing_dims_missing = ndim - shape.len;
            if (out_dim >= trailing_dims_missing) {
                const shape_idx = out_dim - trailing_dims_missing;
                const current_dim = shape[shape_idx];

                if (current_dim == 0)
                    return array.Error.ZeroDimension;

                if (current_dim > max_dim) {
                    if (max_dim != 1)
                        return array.Error.NotBroadcastable;

                    max_dim = current_dim;
                } else if (current_dim != 1 and current_dim != max_dim) {
                    return array.Error.NotBroadcastable;
                }
            }
        }

        result[out_dim] = max_dim;
    }

    return .{
        .ndim = ndim,
        .shape = result,
    };
}

/// Determines the optimal memory layout order (`.c` or `.f`) for an output
/// array so that simultaneous iteration with an input array achieves optimal
/// cache locality.
pub fn orderFromStrides(strides: []const isize) array.Order {
    if (strides.len <= 1)
        return array.Order.default;

    var c_votes: usize = 0;
    var f_votes: usize = 0;
    var prev_stride: ?usize = null;

    for (strides) |stride| {
        const abs_stride = int.abs(stride);

        if (abs_stride == 0)
            continue;

        if (prev_stride) |prev| {
            if (prev > abs_stride) {
                c_votes += 1;
            } else if (prev < abs_stride) {
                f_votes += 1;
            }
        }

        prev_stride = abs_stride;
    }

    if (c_votes > f_votes)
        return .c;

    if (f_votes > c_votes)
        return .f;

    return array.Order.default;
}

/// Checks if the given axes form a valid permutation of `[0, ..., ndim - 1]`.
pub fn isPermutation(axes: []const usize) bool {
    if (axes.len > array.max_dimensions)
        return false;

    var mask: u64 = 0;
    for (axes) |axis| {
        if (axis >= axes.len) return false;

        mask |= (@as(u64, 1) << @intCast(axis));
    }

    const target: u64 = @truncate((@as(u64, 1) << @intCast(axes.len)) - 1);

    return mask == target;
}

pub fn checkIndex(shape: []const usize, index: []const usize) !void {
    if (index.len > shape.len)
        return array.Error.DimensionMismatch;

    var i: usize = 0;
    while (i < index.len) : (i += 1) {
        if (index[i] >= shape[i]) {
            return array.Error.PositionOutOfBounds;
        }
    }
}

inline fn getPtr(ptr: anytype, offset: isize) @TypeOf(ptr) {
    return if (offset >= 0)
        ptr + numeric.cast(usize, offset)
    else
        ptr - numeric.cast(usize, -offset);
}
