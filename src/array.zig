const std = @import("std");
const types = @import("types.zig");
const scast = types.scast;
const cast = types.cast;

const ops = @import("ops.zig");

const int = @import("int.zig");

const dense = @import("array/dense.zig");
const strided = @import("array/strided.zig");

pub const Iterator = @import("array/iterators.zig").Iterator;
//pub const MultiIterator = @import("array/iterators.zig").MultiIterator;
pub const IterationOrder = @import("array/iterators.zig").IterationOrder;

pub const maxDimensions = 8;

pub fn Array(comptime T: type) type {
    // Catch any attempt to create with unsupported type.
    _ = types.numericType(T);

    return struct {
        data: []T,
        ndim: usize,
        shape: [maxDimensions]usize,
        size: usize,
        base: ?*const Array(T),
        flags: Flags,
        metadata: Metadata,

        pub const empty: Array(T) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** maxDimensions,
            .size = 0,
            .base = null,
            .flags = .{},
            .metadata = .{ .dense = .{
                .strides = .{0} ** maxDimensions,
            } },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            shape: []const usize,
            options: struct {
                order: Order = .rowMajor,
                storage: Storage = .dense,
            },
        ) !Array(T) {
            if (shape.len > maxDimensions) {
                return Error.TooManyDimensions;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return Error.ZeroDimension;
                }
            }

            switch (options.storage) {
                .dense => return dense.init(allocator, T, shape, options.order),
                .strided => return Error.InvalidFlags,
            }
        }

        pub fn full(
            allocator: std.mem.Allocator,
            shape: []const usize,
            value: anytype,
            options: struct {
                order: Order = .rowMajor,
                storage: Storage = .dense,
            },
        ) !Array(T) {
            if (shape.len > maxDimensions) {
                return Error.TooManyDimensions;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return Error.ZeroDimension;
                }
            }

            switch (options.storage) {
                .dense => return dense.full(allocator, T, shape, value, options.order),
                .strided => return Error.InvalidFlags,
            }
        }

        pub fn arange(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            step: anytype,
            options: struct {
                writeable: bool = true,
            },
        ) !Array(T) {
            comptime if (types.isComplex(T))
                @compileError("array.arange does not support " ++ @typeName(T));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.arange: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.arange: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(step)))
                @compileError("array.arange: step cannot be complex, got " ++ @typeName(@TypeOf(step)));

            return dense.arange(allocator, T, start, stop, step, options.writeable);
        }

        pub fn linspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            num: usize,
            options: struct {
                writeable: bool = true,
                endpoint: bool = true,
            },
        ) !Array(T) {
            comptime if (types.isComplex(T))
                @compileError("array.linspace does not support " ++ @typeName(T));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.linspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.linspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (num == 0) {
                return Error.ZeroDimension;
            }

            return dense.linspace(allocator, T, start, stop, num, options.writeable, options.endpoint);
        }

        /// Cleans up the array by deinitializing its elements. If the array holds
        /// a fixed precision type this is does not do anything.
        pub fn cleanup(self: *Array(T), allocator: ?std.mem.Allocator) void {
            if (types.isArbitraryPrecision(T)) {
                switch (self.flags.storage) {
                    .dense => dense.cleanup(T, allocator, self.data),
                    .strided => {}, // strided.cleanup(T, self),
                }
            }
        }

        /// Deinitializes the array, freeing its data if it owns it.
        ///
        /// If the array holds an arbitrary precision type, it will not free the
        /// elements; call `cleanup` first to do that.
        pub fn deinit(self: *Array(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.ownsData) {
                allocator.?.free(self.data);
            }

            self.* = undefined;
        }

        pub fn set(self: *Array(T), position: []const usize, value: T) !void {
            if (!self.flags.writeable) {
                return Error.ArrayNotWriteable;
            }

            switch (self.flags.storage) {
                .dense => return dense.set(T, self, position, value),
                .strided => return strided.set(T, self, position, value),
            }
        }

        pub fn get(self: Array(T), position: []const usize) !*T {
            switch (self.flags.storage) {
                .dense => return dense.get(T, self, position),
                .strided => return strided.get(T, self, position),
            }
        }

        pub fn reshape(self: *const Array(T), shape: []const usize) !Array(T) {
            if (shape.len > maxDimensions) {
                return Error.TooManyDimensions;
            }

            if (shape.len == 0) {
                return Error.ZeroDimension;
            }

            switch (self.flags.storage) {
                .dense => return dense.reshape(T, self, shape),
                .strided => return strided.reshape(T, self, shape),
            }
        }

        pub fn ravel(self: *const Array(T)) !Array(T) {
            if (self.ndim == 0) {
                return Error.ZeroDimension;
            }

            switch (self.flags.storage) {
                .dense => return dense.ravel(T, self),
                .strided => return error.NeedDense, // ravel can only be done on dense arrays without copying
            }
        }

        pub fn transpose(
            self: *const Array(T),
            axes: ?[]const usize,
        ) !Array(T) {
            const axess: []const usize =
                axes orelse
                trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axess.len == 0) {
                return Error.ZeroDimension;
            }

            if (axess.len > self.ndim) {
                return Error.TooManyDimensions;
            }

            if (!isPermutation(self.ndim, axess)) {
                return Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]
            }

            switch (self.flags.storage) {
                .dense => return dense.transpose(T, self, axess),
                .strided => return strided.transpose(T, self, axess),
            }
        }

        pub fn slice(self: *const Array(T), ranges: []const Range) !Array(T) {
            if (ranges.len == 0 or ranges.len > self.ndim) {
                return error.DimensionMismatch;
            }

            switch (self.flags.storage) {
                .dense => return dense.slice(T, self, ranges),
                .strided => return strided.slice(T, self, ranges),
            }
        }

        pub fn broadcast(self: *const Array(T), shape: []const usize) !Array(T) {
            if (shape.len > maxDimensions or shape.len < self.ndim) {
                return Error.TooManyDimensions;
            }

            var i: isize = scast(isize, shape.len - 1);
            const diff: isize = scast(isize, shape.len - self.ndim);
            while (i >= 0) : (i -= 1) {
                if (i - diff >= 0 and
                    self.shape[scast(usize, i - diff)] != 1 and
                    self.shape[scast(usize, i - diff)] != shape[scast(usize, i)])
                {
                    return Error.NotBroadcastable;
                }

                if (shape[scast(usize, i)] == 0) {
                    return Error.ZeroDimension;
                }
            }

            switch (self.flags.storage) {
                .dense => return dense.broadcast(T, self, shape),
                .strided => return strided.broadcast(T, self, shape),
            }
        }
    };
}

const arrops = @import("array/ops.zig");
pub const apply1 = arrops.apply1;
pub const apply1_ = arrops.apply1_;
pub const apply1_to = arrops.apply1_to;
pub const apply2 = arrops.apply2;
pub const apply2_ = arrops.apply2_;
pub const apply2_to = arrops.apply2_to;

pub const add = arrops.add;

pub const abs = arrops.abs;
pub const abs_ = arrops.abs_;
pub const ceil = arrops.ceil;
pub const ceil_ = arrops.ceil_;

pub const Broadcast = struct {
    ndim: usize,
    shape: [maxDimensions]usize,
};

pub fn broadcastShapes(
    shapes: []const []const usize,
) !Broadcast {
    if (shapes.len == 0) {
        return Error.ZeroDimension;
    }

    var ndim: usize = 0;
    for (shapes) |shape| {
        if (shape.len == 0) {
            return Error.ZeroDimension;
        }

        if (shape.len > maxDimensions) {
            return Error.TooManyDimensions;
        }

        if (shape.len > ndim) {
            ndim = shape.len;
        }
    }

    var result: [maxDimensions]usize = .{0} ** maxDimensions;
    var i: isize = scast(isize, ndim - 1);
    while (i >= 0) : (i -= 1) {
        var max_dim: usize = 1;
        for (shapes) |shape| {
            const diff: isize = scast(isize, ndim - shape.len);
            if (i - diff >= 0) {
                if (shape[scast(usize, i - diff)] == 0) {
                    return Error.ZeroDimension;
                }

                if (shape[scast(usize, i - diff)] > max_dim) {
                    if (max_dim != 1) {
                        return Error.NotBroadcastable;
                    }

                    max_dim = shape[scast(usize, i - diff)];
                }

                if (shape[scast(usize, i - diff)] != 1 and
                    shape[scast(usize, i - diff)] != max_dim)
                {
                    return Error.NotBroadcastable;
                }
            }
        }

        result[scast(usize, i)] = max_dim;
    }

    return .{
        .ndim = ndim,
        .shape = result,
    };
}

/// Checks if the given axes form a valid permutation of `[0, ..., ndim - 1]`.
pub fn isPermutation(
    ndim: usize,
    axes: []const usize,
) bool {
    if (ndim != axes.len) {
        return false; // axes must match the shape length
    }

    if (ndim == 0 or ndim > maxDimensions) {
        return false; // empty or too many dimensions is not a valid permutation
    }

    var seen: [maxDimensions]bool = .{false} ** maxDimensions;
    for (axes) |axis| {
        if (axis >= ndim) {
            return false; // axis out of bounds
        }

        if (seen[axis]) {
            return false; // duplicate axis
        }

        seen[axis] = true;
    }

    return true; // is a permutation
}

pub fn trivialPermutation(ndim: usize) [maxDimensions]usize {
    var result: [maxDimensions]usize = .{0} ** maxDimensions;

    if (ndim == 0 or ndim > maxDimensions)
        return result; // empty or too many dimensions, return trivial permutation

    for (result[0..ndim], 0..) |*axis, i| {
        axis.* = i;
    }

    return result;
}

pub fn trivialReversePermutation(ndim: usize) [maxDimensions]usize {
    var result: [maxDimensions]usize = .{0} ** maxDimensions;

    if (ndim == 0 or ndim > maxDimensions)
        return result; // empty or too many dimensions, return empty permutation

    for (result[0..ndim], 0..) |*axis, i| {
        axis.* = ndim - i - 1;
    }

    return result;
}

pub const Error = error{
    ArrayNotWriteable,
    TooManyDimensions,
    InvalidFlags,
    InvalidAxes,
    ZeroDimension,
    NotImplemented,
    NotBroadcastable,
    DimensionMismatch,
    PositionOutOfBounds,
    InvalidRange,
    RangeOutOfBounds,
    ZeroStep,
    NeedDense,
};

pub const Flags = packed struct {
    order: Order = .rowMajor,
    storage: Storage = .dense,
    ownsData: bool = true,
    writeable: bool = true,
};

pub const Order = enum(u1) {
    rowMajor,
    columnMajor,
};

pub const Storage = enum(u1) {
    dense,
    strided,
    //csr,
    //csc,
    //coo,
};

pub const Metadata = union(Storage) {
    dense: Dense,
    strided: Strided,
    //csr: CSR,
    //csc: CSC,
    //coo: COO,

    pub const Dense = struct {
        strides: [maxDimensions]usize,
    };

    pub const Strided = struct {
        strides: [maxDimensions]isize,
        offset: usize,
    };

    pub const CSR = struct {
        // # of nonzero elements
        // data holds the data, with the first element (index 0) being the default value, and the rest being the nonzero elements, so row[i] and column[i] refer to data[i+1]
    };

    pub const CSC = struct {};

    pub const COO = struct {};
};

pub const Range = struct {
    start: usize,
    stop: usize,
    step: isize,

    pub const all: Range = .{ .start = 0, .stop = int.max(usize), .step = 1 };

    pub const all_reverse: Range = .{ .start = int.max(usize), .stop = int.max(usize), .step = -1 };

    pub fn init(start: ?usize, stop: ?usize, step: ?isize) !Range {
        const range: Range = .{
            .start = start orelse int.max(usize),
            .stop = stop orelse int.max(usize),
            .step = step orelse 1,
        };

        if (step == 0) {
            return Error.ZeroStep;
        }

        if (((range.step > 0 and range.start >= range.stop) or
            (range.step < 0 and range.start <= range.stop)) and
            (range.start != int.max(usize) and range.stop != int.max(usize)))
        {
            return Error.RangeOutOfBounds;
        }

        return range;
    }

    pub fn single(index: usize) Range {
        return Range{ .start = index, .stop = index + 1, .step = 1 };
    }

    pub fn len(self: Range) usize {
        if (self.start == self.stop) {
            return 0;
        }

        if (self.step > 0) {
            return (self.stop - self.start + scast(usize, self.step) - 1) / scast(usize, self.step);
        }

        return (self.start - self.stop + scast(usize, int.abs(self.step)) - 1) / scast(usize, int.abs(self.step));
    }
};
