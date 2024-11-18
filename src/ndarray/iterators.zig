const std = @import("std");
const zml = @import("../zml.zig");
const core = @import("../core/core.zig");

const ndarray = @import("ndarray.zig");
const NDArray = ndarray.NDArray;
const Flags = ndarray.Flags;
const Error = ndarray.Error;
const MaxDimensions = ndarray.MaxDimensions;

/// Max number of arrays for a multiIteratior.
const MaxArrays: usize = 8;

/// Iterator for one `NDArray`.
///
/// For now, iterators are associated to an array, but if, for example, we just
/// want to use a transpose once, we have to create the transpose view, and
/// then an iterator for that view. If iterators were associated with a shape
/// and a stride, you would only need an iterator for the transpose case, but
/// that iterator would be (considerably) heavier.
///
/// Maybe, make a standard next that doesnt take iteration order and instead
/// calls `nextOrder` with the best order for the array?
pub fn Iterator(comptime T: type) type {
    return struct {
        const Self: type = @This();
        /// The dimensions of the iterator.
        ndim: usize,
        /// The shape of the iterator.
        shape: [MaxDimensions]usize,
        /// The strides of the iterator.
        strides: [MaxDimensions]usize,
        /// The flags of the iterator.
        flags: Flags,
        /// Current position.
        position: [MaxDimensions]usize,
        /// Current index.
        index: usize,

        /// Initializes an iterator for a given array.
        ///
        /// **Description**:
        ///
        /// Initializes an iterator for a given array.
        ///
        /// **Input Parameters**:
        /// - `array`: the array for the Iterator.
        ///
        /// **Return Values**:
        /// - `Self`: the initialized iterator.
        pub fn init(array: NDArray(T)) Self {
            const ndim: usize = array.shape.len;
            var shape: [MaxDimensions]usize = [_]usize{0} ** MaxDimensions;
            var strides: [MaxDimensions]usize = [_]usize{0} ** MaxDimensions;
            for (0..ndim) |i| {
                shape[i] = array.shape[i];
                strides[i] = array.strides[i];
            }

            return Self{
                .ndim = ndim,
                .shape = shape,
                .strides = strides,
                .flags = array.flags,
                .position = [_]usize{0} ** MaxDimensions,
                .index = 0,
            };
        }

        /// Iterates to the next element.
        ///
        /// **Description**:
        ///
        /// Iterates to the next element in the most efficient memory order, and
        /// returns it if it exists, or null otherwise. Also updates `position`
        /// and `index`.
        ///
        /// Due to the nature of the implementation, if the end of the iterator
        /// is reached, apart from `null` being returned, `position` and `index`
        /// are reset to 0, so there is no need to call `reset`.
        ///
        /// **Input Parameters**:
        /// - `self`: the iterator.
        ///
        /// **Return Values**:
        /// - `usize`: the index of the next item.
        /// - `null`: reached end.
        pub fn next(self: *Self) ?usize {
            return self.nextOrder(self.flags.RowMajorContiguous);
        }

        /// Iterates to the next element.
        ///
        /// **Description**:
        ///
        /// Iterates to the next element in the chosen memory order, and returns
        /// it if it exists, or null otherwise. Also updates `position` and
        /// `index`.
        ///
        /// Due to the nature of the implementation, if the end of the iterator
        /// is reached, apart from `null` being returned, `position` and `index`
        /// are reset to 0, so there is no need to call `reset`.
        ///
        /// **Input Parameters**:
        /// - `self`: the iterator.
        /// - `rowMajorOrder`: what order to iterate with. If `true`, iterates
        /// in row major order (right to left, `(...,0,0)->(...,0,1)->...->
        /// (...,0,n)->(...,1,0))`, which is very efficient for arrays stored in
        /// row major order; if `false`, iterates in column major order (left to
        /// right, `(0,0,...)->(1,0,...)->...->(n,0,...)->(0,1,...)`), which is
        /// very efficient for arrays stored in column major order.
        ///
        /// **Return Values**:
        /// - `usize`: the index of the next item.
        /// - `null`: reached end.
        pub fn nextOrder(self: *Self, rowMajorOrder: bool) ?usize {
            var carry: usize = 1;
            var change: i64 = undefined;
            var prev: i64 = undefined;
            var index: i64 = undefined;
            var stride: i64 = undefined;
            var mustBreak: bool = false;
            if (rowMajorOrder) {
                for (0..self.ndim) |i| {
                    prev = @intCast(self.position[self.ndim - i - 1]);
                    self.position[self.ndim - i - 1] += 1;
                    if (self.position[self.ndim - i - 1] >= self.shape[self.ndim - i - 1]) {
                        self.position[self.ndim - i - 1] = 0;
                    } else {
                        carry = 0;
                        mustBreak = true;
                    }
                    // change = @intCast(self.position[self.ndim - i - 1]) - prev;
                    change = @intCast(self.position[self.ndim - i - 1]);
                    change -= prev;
                    // self.index += change * self.strides[self.ndim - i - 1];
                    index = @intCast(self.index);
                    // index += change * self.strides[self.ndim - i - 1];
                    stride = @intCast(self.strides[self.ndim - i - 1]);
                    index += change * stride;
                    self.index = @intCast(index);
                    if (mustBreak) {
                        break;
                    }
                }
            } else {
                for (0..self.ndim) |i| {
                    prev = @intCast(self.position[i]);
                    self.position[i] += 1;
                    if (self.position[i] >= self.shape[i]) {
                        self.position[i] = 0;
                    } else {
                        carry = 0;
                        mustBreak = true;
                    }
                    // change = @intCast(self.position[i]) - prev;
                    change = @intCast(self.position[i]);
                    change -= prev;
                    // self.index += change * self.strides[i];
                    index = @intCast(self.index);
                    // index += change * self.strides[i];
                    stride = @intCast(self.strides[i]);
                    index += change * stride;
                    self.index = @intCast(index);
                    if (mustBreak) {
                        break;
                    }
                }
            }

            if (carry == 1) {
                return null;
            }

            return self.index;
        }
    };
}

/// Iterator for multiple `NDArray`.
pub fn MultiIterator(comptime T: type) type {
    return struct {
        const Self: type = @This();
        /// The number of arrays.
        narray: usize,
        /// The number of dimetions of the broadcast.
        ndim: usize,
        /// Subiterators for the arrays.
        iterators: [MaxArrays]Iterator(T),
        /// Broadcasted shape.
        shape: [MaxDimensions]usize,
        /// Broadcasted strides.
        strides: [MaxDimensions]usize,
        /// Broadcasted size.
        size: usize,
        /// Flags of the (theoretical) broadcasted array.
        flags: Flags,
        /// Current position.
        position: [MaxDimensions]usize,
        /// Current index.
        index: usize,

        /// Initializes a multi iterator for the given arrays.
        ///
        /// **Description**:
        ///
        /// Initializes a multi iterator for the given arrays and flags.
        ///
        /// **Input Parameters**:
        /// - `array`: the array for the Iterator.
        /// - `flags`: the flags for the iterator.
        ///
        /// **Return Values**:
        /// - `Self`: the initialized iterator.
        pub fn init(arrays: []const NDArray(T), flags: Flags) !Self {
            if (arrays.len == 0) {
                return IteratorError.NoArrays;
            }

            const narray: usize = arrays.len;
            if (narray > MaxArrays) {
                return IteratorError.TooManyArrays;
            }

            var ndim: usize = 0;
            var iterators: [MaxArrays]Iterator(T) = undefined;
            for (0..narray) |i| {
                iterators[i] = Iterator(T).init(arrays[i]);
                if (arrays[i].shape.len > ndim) {
                    ndim = arrays[i].shape.len;
                }
            }

            var shape: [MaxDimensions]usize = [_]usize{0} ** MaxDimensions;
            for (0..narray) |i| {
                //std.debug.print("Shape = [  ", .{});
                //for (0..ndim) |help| {
                //    std.debug.print("{}  ", .{shape[help]});
                //}
                //std.debug.print("]\n", .{});
                //std.debug.print("Array shape = [  ", .{});
                //for (0..arrays[i].shape.len) |help| {
                //    std.debug.print("{}  ", .{arrays[i].shape[help]});
                //}
                //std.debug.print("]\n", .{});
                for (0..arrays[i].shape.len) |j| {
                    if (shape[ndim - j - 1] != arrays[i].shape[arrays[i].shape.len - j - 1]) {
                        //std.debug.print("Disagreement = {} != {}\n", .{ shape[ndim - j - 1], arrays[i].shape[arrays[i].shape.len - j - 1] });
                        if (shape[ndim - j - 1] == 1 or shape[ndim - j - 1] == 0) {
                            shape[ndim - j - 1] = arrays[i].shape[arrays[i].shape.len - j - 1];
                        } else if (arrays[i].shape[arrays[i].shape.len - j - 1] != 1) {
                            return IteratorError.NotBroadcastable;
                        }
                    }
                }
            }

            var size: usize = 1;
            var strides: [MaxDimensions]usize = [_]usize{0} ** MaxDimensions;
            if (flags.order == .RowMajor) {
                for (0..ndim) |i| {
                    strides[ndim - i - 1] = size;
                    size *= shape[ndim - i - 1];
                }
            } else {
                for (0..ndim) |i| {
                    strides[i] = size;
                    size *= shape[i];
                }
            }

            return Self{
                .narray = narray,
                .ndim = ndim,
                .iterators = iterators,
                .shape = shape,
                .strides = strides,
                .size = size,
                .flags = flags,
                .position = [_]usize{0} ** MaxDimensions,
                .index = 0,
            };
        }

        /// Initializes a multi iterator for the given arrays.
        ///
        /// **Description**:
        ///
        /// Initializes a multi iterator for the given arrays. The first array
        /// must have the correct broadcasted shape.
        ///
        /// **Input Parameters**:
        /// - `array`: the array for the Iterator.
        /// - `flags`: the flags for the iterator.
        ///
        /// **Return Values**:
        /// - `Self`: the initialized iterator.
        pub fn initCheck(first: NDArray(T), rest: []const NDArray(T)) !Self {
            // Not implemented yet, this is done to be able to compile.
            return init([_]NDArray(T){first} ++ rest, first.flags);
        }

        /// Iterates to the next element.
        ///
        /// **Description**:
        ///
        /// Iterates to the next element in the most efficient memory order, and
        /// returns it if it exists, or null otherwise. Also updates `position`
        /// and `index`.
        ///
        /// Due to the nature of the implementation, if the end of the iterator
        /// is reached, apart from `null` being returned, `position` and `index`
        /// are reset to 0, so there is no need to call `reset`.
        ///
        /// **Input Parameters**:
        /// - `self`: the iterator.
        ///
        /// **Return Values**:
        /// - `T`: the next item.
        /// - `null`: reached end.
        pub fn next(self: *Self) ?usize {
            var trueCount: usize = 0;
            for (0..self.narray) |i| {
                if (self.iterators[i].flags.order == .RowMajor) {
                    trueCount += 1;
                }
            }
            const order: bool = self.narray <= (trueCount * 2);

            return self.nextOrder(order);
        }

        /// Iterates to the next element.
        ///
        /// **Description**:
        ///
        /// Iterates to the next element in the chosen memory order, and returns
        /// it if it exists, or null otherwise. Also updates `position` and
        /// `index`.
        ///
        /// Due to the nature of the implementation, if the end of the iterator
        /// is reached, apart from `null` being returned, `position` and `index`
        /// are reset to 0, so there is no need to call `reset`.
        ///
        /// **Input Parameters**:
        /// - `self`: the iterator.
        /// - `rowMajorOrder`: what order to iterate with. If `true`, iterates
        /// in row major order (right to left, `(...,0,0)->(...,0,1)->...->
        /// (...,0,n)->(...,1,0))`, which is very efficient for arrays stored in
        /// row major order; if `false`, iterates in column major order (left to
        /// right, `(0,0,...)->(1,0,...)->...->(n,0,...)->(0,1,...)`), which is
        /// very efficient for arrays stored in column major order.
        ///
        /// **Return Values**:
        /// - `T`: the next item.
        /// - `null`: reached end.
        pub fn nextOrder(self: *Self, rowMajorOrder: bool) ?usize {
            var carry: usize = 1;
            var change: i64 = undefined;
            var prev: i64 = undefined;
            var index: i64 = undefined;
            var stride: i64 = undefined;
            var mustBreak: bool = false;
            if (rowMajorOrder) {
                for (0..self.ndim) |i| {
                    prev = @intCast(self.position[self.ndim - i - 1]);
                    self.position[self.ndim - i - 1] += 1;
                    if (self.position[self.ndim - i - 1] >= self.shape[self.ndim - i - 1]) {
                        self.position[self.ndim - i - 1] = 0;
                    } else {
                        carry = 0;
                        mustBreak = true;
                    }
                    // change = @intCast(self.position[self.ndim - i - 1]) - prev;
                    change = @intCast(self.position[self.ndim - i - 1]);
                    change -= prev;
                    // self.index += change * self.strides[self.ndim - i - 1];
                    index = @intCast(self.index);
                    // index += change * self.strides[self.ndim - i - 1];
                    stride = @intCast(self.strides[self.ndim - i - 1]);
                    index += change * stride;
                    self.index = @intCast(index);
                    if (mustBreak) {
                        break;
                    }
                }
            } else {
                for (0..self.ndim) |i| {
                    prev = @intCast(self.position[i]);
                    self.position[i] += 1;
                    if (self.position[i] >= self.shape[i]) {
                        self.position[i] = 0;
                    } else {
                        carry = 0;
                        mustBreak = true;
                    }
                    // change = @intCast(self.position[i]) - prev;
                    change = @intCast(self.position[i]);
                    change -= prev;
                    // self.index += change * self.strides[i];
                    index = @intCast(self.index);
                    // index += change * self.strides[i];
                    stride = @intCast(self.strides[i]);
                    index += change * stride;
                    self.index = @intCast(index);
                    if (mustBreak) {
                        break;
                    }
                }
            }

            for (0..self.narray) |i| {
                for (0..self.iterators[i].ndim) |j| {
                    //std.debug.print("ndim:{},i:[{}],", .{ self.iterators[i].ndim, i });
                    prev = @intCast(self.iterators[i].position[self.iterators[i].ndim - j - 1]);
                    if (self.iterators[i].shape[self.iterators[i].ndim - j - 1] == 1) {
                        continue;
                    }
                    self.iterators[i].position[self.iterators[i].ndim - j - 1] = self.position[self.ndim - j - 1];
                    // change = @intCast(self.position[length - i - 1]) - prev;
                    change = @intCast(self.iterators[i].position[self.iterators[i].ndim - j - 1]);
                    change -= prev;
                    // self.index += change * self.array.strides[length - i - 1];
                    index = @intCast(self.iterators[i].index);
                    // index += change * self.array.strides[length - i - 1];
                    stride = @intCast(self.iterators[i].strides[self.iterators[i].ndim - j - 1]);
                    index += change * stride;
                    self.iterators[i].index = @intCast(index);
                }
                //std.debug.print("\n", .{});
            }

            if (carry == 1) {
                return null;
            }

            return self.index;
        }
    };
}

/// Errors that can occur when workin with array iterators.
pub const IteratorError = error{
    /// Too many arrays.
    TooManyArrays,
    /// No input arrays.
    NoArrays,
    /// Incompatible shapes for broadcasting.
    NotBroadcastable,
};
