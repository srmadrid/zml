const std = @import("std");
const types = @import("../types.zig");
const scast = types.scast;
const IterationOrder = types.IterationOrder;

const array = @import("../array.zig");

const max_arrays: usize = 8;

pub fn Iterator(comptime T: type) type {
    return struct {
        arr: array.AnyArray(T),
        position: [array.max_dimensions]usize,
        index: usize,

        pub fn init(arr: array.AnyArray(T)) Iterator(T) {
            switch (arr) {
                .dense => {
                    return Iterator(T){
                        .arr = arr,
                        .position = .{0} ** array.max_dimensions,
                        .index = 0,
                    };
                },
                .strided => {
                    return Iterator(T){
                        .arr = arr,
                        .position = .{0} ** array.max_dimensions,
                        .index = arr.strided.offset,
                    };
                },
                .sparse => {
                    return Iterator(T){
                        .arr = arr,
                        .position = .{0} ** array.max_dimensions,
                        .index = 0,
                    };
                },
            }
        }

        pub fn next(self: *Iterator(T)) ?usize {
            switch (self.arr) {
                .dense, .strided => |arr| {
                    switch (arr.flags.order) {
                        .row_major => {
                            return self.nextAO(arr.ndim - 1, .right_to_left);
                        },
                        .col_major => {
                            return self.nextAO(0, .left_to_right);
                        },
                    }
                },
                .sparse => {
                    return 0;
                },
            }
        }

        pub fn nextAO(self: *Iterator(T), axis: usize, order: IterationOrder) ?usize {
            switch (self.arr) {
                .dense => |arr| {
                    var carry: bool = true;
                    var mustBreak: bool = false;
                    var change: isize = undefined;
                    var prev: isize = undefined;
                    var index: isize = undefined;
                    var stride: isize = undefined;

                    const ax = if (order == .right_to_left) arr.ndim - axis - 1 else axis;

                    for (ax..arr.ndim) |i| {
                        const idx = if (order == .right_to_left) arr.ndim - i - 1 else i;
                        prev = scast(isize, self.position[idx]);
                        self.position[idx] += 1;

                        if (self.position[idx] >= arr.shape[idx]) {
                            self.position[idx] = 0;
                        } else {
                            carry = false;
                            mustBreak = true;
                        }

                        change = scast(isize, self.position[idx]) - prev;
                        index = scast(isize, self.index);
                        stride = scast(isize, arr.strides[idx]);
                        index += change * stride;
                        self.index = scast(usize, index);

                        if (mustBreak) {
                            break;
                        }
                    }

                    if (carry) {
                        return null;
                    }

                    return self.index;
                },
                .strided => |arr| {
                    var carry: bool = true;
                    var mustBreak: bool = false;
                    var change: isize = undefined;
                    var prev: isize = undefined;
                    var index: isize = undefined;
                    var stride: isize = undefined;

                    const ax = if (order == .right_to_left) arr.ndim - axis - 1 else axis;

                    for (ax..arr.ndim) |i| {
                        const idx = if (order == .right_to_left) arr.ndim - i - 1 else i;
                        prev = scast(isize, self.position[idx]);
                        self.position[idx] += 1;

                        if (self.position[idx] >= arr.shape[idx]) {
                            self.position[idx] = 0;
                        } else {
                            carry = false;
                            mustBreak = true;
                        }

                        change = scast(isize, self.position[idx]) - prev;
                        index = scast(isize, self.index);
                        stride = arr.strides[idx];
                        index += change * stride;
                        self.index = scast(usize, index);

                        if (mustBreak) {
                            break;
                        }
                    }

                    if (carry) {
                        return null;
                    }

                    return self.index;
                },
                .sparse => {
                    // Sparse arrays do not have a defined order, so we return null.
                    return null;
                },
            }
        }
    };
}

// /// Iterator for multiple `Array`.
// pub fn MultiIterator(comptime T: type) type {
//     return struct {
//         /// The number of arrs.
//         narr: usize,
//         /// The number of dimetions of the broadscast.
//         ndim: usize,
//         /// Subiterators for the arrs.
//         iterators: [max_arrays]Iterator(T),
//         /// Broadscasted shape.
//         shape: [array.max_dimensions]usize,
//         /// Broadscasted strides.
//         strides: [array.max_dimensions]isize,
//         /// Offset of the first element of the iterator.
//         offset: usize,
//         /// Broadscasted size.
//         size: usize,
//         /// Flags of the (theoretical) broadscasted arr.
//         flags: array.Flags,
//         /// Current position.
//         position: [array.max_dimensions]usize,
//         /// Current index.
//         index: usize,

//         /// Initializes a multi iterator for the given arrs.
//         ///
//         /// **Description**:
//         ///
//         /// Initializes a multi iterator for the given arrs and flags.
//         ///
//         /// **Input Parameters**:
//         /// - `arr`: the arr for the Iterator.
//         /// - `flags`: the flags for the iterator.
//         ///
//         /// **Return Values**:
//         /// - `MultiIterator(T)`: the initialized iterator.
//         pub fn init(arrs: []const Array(T), flags: array.Flags) !MultiIterator(T) {
//             if (arrs.len == 0) {
//                 return Error.NoArrays;
//             }

//             const narr: usize = arrs.len;
//             if (narr > max_arrays) {
//                 return Error.TooManyArrays;
//             }

//             var ndim: usize = 0;
//             var iterators: [max_arrays]Iterator(T) = undefined;
//             for (0..narr) |i| {
//                 iterators[i] = Iterator(T).init(arrs[i]);
//                 if (arrs[i].ndim > ndim) {
//                     ndim = arrs[i].ndim;
//                 }
//             }

//             var shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
//             for (0..narr) |i| {
//                 for (0..arrs[i].ndim) |j| {
//                     if (shape[ndim - j - 1] != arrs[i].shape[arrs[i].ndim - j - 1]) {
//                         //std.debug.print("Disagreement = {} != {}\n", .{ shape[ndim - j - 1], arrs[i].shape[arrs[i].shape.len - j - 1] });
//                         if (shape[ndim - j - 1] == 1 or shape[ndim - j - 1] == 0) {
//                             shape[ndim - j - 1] = arrs[i].shape[arrs[i].ndim - j - 1];
//                         } else if (arrs[i].shape[arrs[i].ndim - j - 1] != 1) {
//                             return Error.NotBroadscastable;
//                         }
//                     }
//                 }
//             }

//             var size: usize = 1;
//             var strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;

//             for (0..ndim) |i| {
//                 const idx = if (flags.order == .RowMajor) ndim - i - 1 else i;
//                 strides[idx] = @intCast(size);
//                 size *= shape[idx];
//             }

//             return MultiIterator(T){
//                 .narr = narr,
//                 .ndim = ndim,
//                 .iterators = iterators,
//                 .shape = shape,
//                 .strides = strides,
//                 .offset = 0,
//                 .size = size,
//                 .flags = flags,
//                 .position = [_]usize{0} ** array.max_dimensions,
//                 .index = 0,
//             };
//         }

//         /// Initializes a multi iterator for the given arrs.
//         ///
//         /// **Description**:
//         ///
//         /// Initializes a multi iterator for the given arrs. The first arr
//         /// must have the correct broadscasted shape.
//         ///
//         /// **Input Parameters**:
//         /// - `arr`: the arr for the Iterator.
//         /// - `flags`: the flags for the iterator.
//         ///
//         /// **Return Values**:
//         /// - `MultiIterator(T)`: the initialized iterator.
//         pub fn initCheck(first: Array(T), rest: []const Array(T)) !MultiIterator(T) {
//             // Not implemented yet, this is done to be able to compile.
//             return init([_]Array(T){first} ++ rest, first.flags);
//         }

//         /// Iterates to the next element.
//         ///
//         /// **Description**:
//         ///
//         /// Iterates to the next element in the most efficient memory order, and
//         /// returns it if it exists, or null otherwise. Also updates `position`
//         /// and `index`.
//         ///
//         /// Due to the nature of the implementation, if the end of the iterator
//         /// is reached, apart from `null` being returned, `position` and `index`
//         /// are reset to 0, so there is no need to call `reset`.
//         ///
//         /// **Input Parameters**:
//         /// - `self`: the iterator.
//         ///
//         /// **Return Values**:
//         /// - `T`: the next item.
//         /// - `null`: reached end.
//         pub fn next(self: *MultiIterator(T)) ?usize {
//             var rowCount: usize = 0;
//             for (0..self.narr) |i| {
//                 if (self.iterators[i].flags.order == .RowMajor) {
//                     rowCount += 1;
//                 }
//             }
//             const order: array.Order = if (self.narr <= (rowCount * 2)) .RowMajor else .ColumnMajor;

//             return self.nextOrder(order);
//         }

//         /// Iterates to the next element.
//         ///
//         /// **Description**:
//         ///
//         /// Iterates to the next element in the chosen memory order, and returns
//         /// it if it exists, or null otherwise. Also updates `position` and
//         /// `index`.
//         ///
//         /// Due to the nature of the implementation, if the end of the iterator
//         /// is reached, apart from `null` being returned, `position` and `index`
//         /// are reset to 0, so there is no need to call `reset`.
//         ///
//         /// **Input Parameters**:
//         /// - `self`: the iterator.
//         /// - `row_majorOrder`: what order to iterate with. If `true`, iterates
//         /// in row major order (right to left, `(...,0,0)->(...,0,1)->...->
//         /// (...,0,n)->(...,1,0))`, which is very efficient for arrs stored in
//         /// row major order; if `false`, iterates in column major order (left to
//         /// right, `(0,0,...)->(1,0,...)->...->(n,0,...)->(0,1,...)`), which is
//         /// very efficient for arrs stored in column major order.
//         ///
//         /// **Return Values**:
//         /// - `T`: the next item.
//         /// - `null`: reached end.
//         pub fn nextOrder(self: *MultiIterator(T), order: array.Order) ?usize {
//             var carry: usize = 1;
//             var change: i64 = undefined;
//             var prev: i64 = undefined;
//             var index: i64 = undefined;
//             var stride: i64 = undefined;
//             var mustBreak: bool = false;
//             if (order == .RowMajor) {
//                 for (0..self.ndim) |i| {
//                     prev = @intCast(self.position[self.ndim - i - 1]);
//                     self.position[self.ndim - i - 1] += 1;
//                     if (self.position[self.ndim - i - 1] >= self.shape[self.ndim - i - 1]) {
//                         self.position[self.ndim - i - 1] = 0;
//                     } else {
//                         carry = 0;
//                         mustBreak = true;
//                     }
//                     // change = @intCast(self.position[self.ndim - i - 1]) - prev;
//                     change = @intCast(self.position[self.ndim - i - 1]);
//                     change -= prev;
//                     // self.index += change * self.strides[self.ndim - i - 1];
//                     index = @intCast(self.index);
//                     // index += change * self.strides[self.ndim - i - 1];
//                     stride = @intCast(self.strides[self.ndim - i - 1]);
//                     index += change * stride;
//                     self.index = @intCast(index);
//                     if (mustBreak) {
//                         break;
//                     }
//                 }
//             } else {
//                 for (0..self.ndim) |i| {
//                     prev = @intCast(self.position[i]);
//                     self.position[i] += 1;
//                     if (self.position[i] >= self.shape[i]) {
//                         self.position[i] = 0;
//                     } else {
//                         carry = 0;
//                         mustBreak = true;
//                     }
//                     // change = @intCast(self.position[i]) - prev;
//                     change = @intCast(self.position[i]);
//                     change -= prev;
//                     // self.index += change * self.strides[i];
//                     index = @intCast(self.index);
//                     // index += change * self.strides[i];
//                     stride = @intCast(self.strides[i]);
//                     index += change * stride;
//                     self.index = @intCast(index);
//                     if (mustBreak) {
//                         break;
//                     }
//                 }
//             }

//             for (0..self.narr) |i| {
//                 for (0..self.iterators[i].ndim) |j| {
//                     //std.debug.print("ndim:{},i:[{}],", .{ self.iterators[i].ndim, i });
//                     prev = @intCast(self.iterators[i].position[self.iterators[i].ndim - j - 1]);
//                     if (self.iterators[i].shape[self.iterators[i].ndim - j - 1] == 1) {
//                         continue;
//                     }
//                     self.iterators[i].position[self.iterators[i].ndim - j - 1] = self.position[self.ndim - j - 1];
//                     // change = @intCast(self.position[length - i - 1]) - prev;
//                     change = @intCast(self.iterators[i].position[self.iterators[i].ndim - j - 1]);
//                     change -= prev;
//                     // self.index += change * self.arr.strides[length - i - 1];
//                     index = @intCast(self.iterators[i].index);
//                     // index += change * self.arr.strides[length - i - 1];
//                     stride = @intCast(self.iterators[i].strides[self.iterators[i].ndim - j - 1]);
//                     index += change * stride;
//                     self.iterators[i].index = @intCast(index);
//                 }
//                 //std.debug.print("\n", .{});
//             }

//             if (carry == 1) {
//                 return null;
//             }

//             return self.index;
//         }
//     };
// }

/// Errors that can occur when workin with arr iterators.
pub const Error = error{
    /// Too many arrs.
    TooManyArrays,
    /// No input arrs.
    NoArrays,
    /// Incompatible shapes for broadscasting.
    NotBroadscastable,
};
