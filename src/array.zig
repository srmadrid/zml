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

            return dense.linspace(allocator, T, start, stop, num, options.writeable);
        }

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

        pub fn get(self: *const Array(T), position: []const usize) !*T {
            switch (self.flags.storage) {
                .dense => return dense.get(T, self, position),
                .strided => return strided.get(T, self, position),
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
    };
}

pub const Error = error{
    ArrayNotWriteable,
    TooManyDimensions,
    InvalidFlags,
    ZeroDimension,
    NotImplemented,
    DimensionMismatch,
    PositionOutOfBounds,
    InvalidRange,
    RangeOutOfBounds,
    ZeroStep,
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

        if ((range.step > 0 and range.start >= range.stop) or
            (range.step < 0 and range.start <= range.stop))
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
