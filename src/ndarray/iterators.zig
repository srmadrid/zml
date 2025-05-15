const std = @import("std");
const types = @import("../types.zig");
const cast = types.cast;

const ndarray = @import("../ndarray.zig");
const NDArray = ndarray.NDArray;

const maxArrays: usize = 8;

pub fn Iterator(comptime T: type) type {
    return struct {
        array: *const NDArray(T),
        position: [ndarray.maxDimensions]usize,
        index: usize,

        pub fn init(array: *const NDArray(T)) Iterator(T) {
            switch (array.flags.storage) {
                .dense => {
                    return Iterator(T){
                        .array = array,
                        .position = .{0} ** ndarray.maxDimensions,
                        .index = 0,
                    };
                },
                .strided => {
                    return Iterator(T){
                        .array = array,
                        .position = .{0} ** ndarray.maxDimensions,
                        .index = array.metadata.strided.offset,
                    };
                },
                .csr => {
                    // Check if the elemnt in (0, 0) is default or a value. If
                    // it is a value, then index = 1 (or wherever it is),
                    // otherwise index = 0 (the position for the default value).
                    std.debug.print("CSR storage not implemented yet", .{});
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                },
            }

            return Iterator(T){
                .array = array,
                .position = .{0} ** ndarray.maxDimensions,
                .index = 0,
            };
        }

        pub fn next(self: *Iterator(T)) ?usize {
            switch (self.array.flags.storage) {
                .dense, .strided => {
                    switch (self.array.flags.order) {
                        .rowMajor => {
                            return self.nextAO(self.array.ndim - 1, .rightToLeft);
                        },
                        .columnMajor => {
                            return self.nextAO(0, .leftToRight);
                        },
                        .other => unreachable,
                    }
                },
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                    return null;
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                    return null;
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                    return null;
                },
            }
        }

        pub fn nextAO(self: *Iterator(T), axis: usize, order: IterationOrder) ?usize {
            switch (self.array.flags.storage) {
                .dense => {
                    var carry: bool = true;
                    var mustBreak: bool = false;
                    var change: isize = undefined;
                    var prev: isize = undefined;
                    var index: isize = undefined;
                    var stride: isize = undefined;

                    const ax = if (order == .rightToLeft) self.array.ndim - axis - 1 else axis;

                    for (ax..self.array.ndim) |i| {
                        const idx = if (order == .rightToLeft) self.array.ndim - i - 1 else i;
                        prev = @intCast(self.position[idx]);
                        self.position[idx] += 1;

                        if (self.position[idx] >= self.array.shape[idx]) {
                            self.position[idx] = 0;
                        } else {
                            carry = false;
                            mustBreak = true;
                        }

                        change = cast(isize, self.position[idx], .{});
                        change -= prev;
                        index = @intCast(self.index);
                        stride = @intCast(self.array.metadata.dense.strides[idx]);
                        index += change * stride;
                        self.index = @intCast(index);

                        if (mustBreak) {
                            break;
                        }
                    }

                    if (carry) {
                        return null;
                    }

                    return self.index;
                },
                .strided => {
                    return null;
                },
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                    return null;
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                    return null;
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                    return null;
                },
            }
        }
    };
}
/// Iterator for multiple `NDArray`.
pub fn MultiIterator(comptime T: type) type {
    return struct {
        /// The number of arrays.
        narray: usize,
        /// The number of dimetions of the broadcast.
        ndim: usize,
        /// Subiterators for the arrays.
        iterators: [maxArrays]Iterator(T),
        /// Broadcasted shape.
        shape: [ndarray.maxDimensions]usize,
        /// Broadcasted strides.
        strides: [ndarray.maxDimensions]isize,
        /// Offset of the first element of the iterator.
        offset: usize,
        /// Broadcasted size.
        size: usize,
        /// Flags of the (theoretical) broadcasted array.
        flags: ndarray.Flags,
        /// Current position.
        position: [ndarray.maxDimensions]usize,
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
        /// - `MultiIterator(T)`: the initialized iterator.
        pub fn init(arrays: []const NDArray(T), flags: ndarray.Flags) !MultiIterator(T) {
            if (arrays.len == 0) {
                return Error.NoArrays;
            }

            const narray: usize = arrays.len;
            if (narray > maxArrays) {
                return Error.TooManyArrays;
            }

            var ndim: usize = 0;
            var iterators: [maxArrays]Iterator(T) = undefined;
            for (0..narray) |i| {
                iterators[i] = Iterator(T).init(arrays[i]);
                if (arrays[i].ndim > ndim) {
                    ndim = arrays[i].ndim;
                }
            }

            var shape: [ndarray.maxDimensions]usize = .{0} ** ndarray.maxDimensions;
            for (0..narray) |i| {
                for (0..arrays[i].ndim) |j| {
                    if (shape[ndim - j - 1] != arrays[i].shape[arrays[i].ndim - j - 1]) {
                        //std.debug.print("Disagreement = {} != {}\n", .{ shape[ndim - j - 1], arrays[i].shape[arrays[i].shape.len - j - 1] });
                        if (shape[ndim - j - 1] == 1 or shape[ndim - j - 1] == 0) {
                            shape[ndim - j - 1] = arrays[i].shape[arrays[i].ndim - j - 1];
                        } else if (arrays[i].shape[arrays[i].ndim - j - 1] != 1) {
                            return Error.NotBroadcastable;
                        }
                    }
                }
            }

            var size: usize = 1;
            var strides: [ndarray.maxDimensions]isize = .{0} ** ndarray.maxDimensions;

            for (0..ndim) |i| {
                const idx = if (flags.order == .RowMajor) ndim - i - 1 else i;
                strides[idx] = @intCast(size);
                size *= shape[idx];
            }

            return MultiIterator(T){
                .narray = narray,
                .ndim = ndim,
                .iterators = iterators,
                .shape = shape,
                .strides = strides,
                .offset = 0,
                .size = size,
                .flags = flags,
                .position = [_]usize{0} ** ndarray.maxDimensions,
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
        /// - `MultiIterator(T)`: the initialized iterator.
        pub fn initCheck(first: NDArray(T), rest: []const NDArray(T)) !MultiIterator(T) {
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
        pub fn next(self: *MultiIterator(T)) ?usize {
            var rowCount: usize = 0;
            for (0..self.narray) |i| {
                if (self.iterators[i].flags.order == .RowMajor) {
                    rowCount += 1;
                }
            }
            const order: ndarray.Order = if (self.narray <= (rowCount * 2)) .RowMajor else .ColumnMajor;

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
        pub fn nextOrder(self: *MultiIterator(T), order: ndarray.Order) ?usize {
            var carry: usize = 1;
            var change: i64 = undefined;
            var prev: i64 = undefined;
            var index: i64 = undefined;
            var stride: i64 = undefined;
            var mustBreak: bool = false;
            if (order == .RowMajor) {
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
pub const Error = error{
    /// Too many arrays.
    TooManyArrays,
    /// No input arrays.
    NoArrays,
    /// Incompatible shapes for broadcasting.
    NotBroadcastable,
};

pub const IterationOrder = enum {
    leftToRight,
    rightToLeft,
};
