const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;
const int = @import("int.zig");

const dense = @import("ndarray/dense.zig");
const strided = @import("ndarray/strided.zig");

pub const Iterator = @import("ndarray/iterators.zig").Iterator;
//pub const MultiIterator = @import("ndarray/iterators.zig").MultiIterator;

pub const maxDimensions = 8;

pub fn NDArray(comptime T: type) type {
    // Catch any attempt to create with unsupported type.
    _ = types.numericType(T);

    return struct {
        data: []T,
        ndim: usize,
        shape: [maxDimensions]usize,
        size: usize,
        base: ?*const NDArray(T),
        flags: Flags,
        metadata: Metadata,

        pub const empty: NDArray(T) = .{
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

        pub fn init(allocator: std.mem.Allocator, shape: []const usize, flags: Flags) !NDArray(T) {
            if (shape.len > maxDimensions) {
                return Error.TooManyDimensions;
            }

            if (!flags.ownsData) {
                return Error.InvalidFlags;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return Error.ZeroDimension;
                }
            }

            switch (flags.storage) {
                .dense => {
                    if (flags.order == .other) {
                        return Error.InvalidFlags;
                    }

                    var size: usize = 1;
                    var shapes: [maxDimensions]usize = .{0} ** maxDimensions;
                    var strides: [maxDimensions]usize = .{0} ** maxDimensions;
                    if (shape.len > 0) {
                        for (0..shape.len) |i| {
                            const idx: usize = if (flags.order == .rowMajor) shape.len - i - 1 else i;

                            strides[idx] = size;
                            size *= shape[idx];

                            shapes[i] = shape[i];
                        }
                    }

                    return NDArray(T){
                        .data = try allocator.alloc(T, size),
                        .ndim = shape.len,
                        .shape = shapes,
                        .size = size,
                        .base = null,
                        .flags = flags,
                        .metadata = .{ .dense = .{
                            .strides = strides,
                        } },
                    };
                },
                .strided => return Error.InvalidFlags,
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
            }
        }

        pub fn deinit(self: *NDArray(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.ownsData) {
                allocator.?.free(self.data);
            }

            self.ndim = 0;
            self.shape = .{0} ** maxDimensions;
            self.size = 0;
            self.base = null;

            switch (self.flags.storage) {
                .dense => {
                    self.metadata.dense.strides = .{0} ** maxDimensions;
                },
                .strided => {
                    self.metadata.strided.strides = .{0} ** maxDimensions;
                    self.metadata.strided.offset = 0;
                },
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                },
            }
        }

        pub fn set(self: *NDArray(T), position: []const usize, value: T) !void {
            _ = types.numericType(@TypeOf(value));

            switch (self.flags.storage) {
                .dense => {
                    if (self.size == 1) {
                        self.data[0] = value;
                        return;
                    }

                    try dense.checkPosition(T, self, position);

                    self.data[dense.index(T, self, position)] = value;
                },
                .strided => {
                    if (self.size == 1) {
                        self.data[self.metadata.strided.offset] = value;
                        return;
                    }

                    try strided.checkPosition(T, self, position);

                    self.data[strided.index(T, self, position)] = value;
                },
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
            }
        }

        pub fn get(self: *const NDArray(T), position: []const usize) !*T {
            switch (self.flags.storage) {
                .dense => {
                    if (self.size == 1) {
                        return &self.data[0];
                    }

                    try dense.checkPosition(T, self, position);

                    return &self.data[dense.index(T, self, position)];
                },
                .strided => {
                    if (self.size == 1) {
                        return &self.data[self.metadata.strided.offset];
                    }

                    try strided.checkPosition(T, self, position);

                    return &self.data[strided.index(T, self, position)];
                },
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
            }
        }

        pub fn slice(self: *const NDArray(T), slices: []const Slice) !NDArray(T) {
            switch (self.flags.storage) {
                .dense, .strided => {
                    if (slices.len == 0 or slices.len > self.ndim) {
                        return error.DimensionMismatch;
                    }

                    var ndim: usize = self.ndim;
                    var size: usize = 1;
                    var shapes: [maxDimensions]usize = .{0} ** maxDimensions;
                    var strides: [maxDimensions]isize = .{0} ** maxDimensions;
                    var offset: usize = if (self.flags.storage == .dense) 0 else self.metadata.strided.offset;

                    var i: usize = 0;
                    var j: usize = 0;
                    while (i < self.ndim) {
                        const stride: isize = if (self.flags.storage == .dense)
                            cast(isize, self.metadata.dense.strides[i], .{})
                        else
                            self.metadata.strided.strides[i];

                        if (i >= slices.len) {
                            shapes[j] = self.shape[i];
                            strides[j] = stride;
                            size *= self.shape[i];
                            j += 1;
                            i += 1;
                            continue;
                        }

                        if (slices[i].step > 0) {
                            if (slices[i].start >= self.shape[i] or slices[i].stop > self.shape[i]) {
                                return Error.SliceOutOfBounds;
                            }
                        } else if (slices[i].step == 0) {
                            if (slices[i].start >= self.shape[i]) {
                                return Error.SliceOutOfBounds;
                            }
                        } else {
                            if (slices[i].stop >= self.shape[i] or slices[i].start > self.shape[i]) {
                                return Error.SliceOutOfBounds;
                            }
                        }

                        const len: usize = slices[i].len();
                        if (len == 0 or len == 1) {
                            ndim -= 1;
                        } else {
                            shapes[j] = len;
                            strides[j] = stride * slices[i].step;
                            size *= len;
                            j += 1;
                        }

                        if (stride < 0) {
                            offset -= slices[i].start * cast(usize, int.abs(stride), .{});
                        } else {
                            offset += slices[i].start * cast(usize, stride, .{});
                        }

                        i += 1;
                    }

                    return NDArray(T){
                        .data = self.data,
                        .ndim = ndim,
                        .shape = shapes,
                        .size = size,
                        .base = if (self.flags.ownsData) self else self.base,
                        .flags = .{
                            .order = .other,
                            .storage = .strided,
                            .ownsData = false,
                            .writeable = self.flags.writeable,
                        },
                        .metadata = .{ .strided = .{
                            .strides = strides,
                            .offset = offset,
                        } },
                    };
                },
                .csr => {
                    std.debug.print("CSR storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .csc => {
                    std.debug.print("CSC storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
                .coo => {
                    std.debug.print("COO storage not implemented yet", .{});
                    return Error.NotImplemented;
                },
            }
        }
    };
}

pub const Error = error{
    TooManyDimensions,
    InvalidFlags,
    ZeroDimension,
    NotImplemented,
    DimensionMismatch,
    PositionOutOfBounds,
    SliceOutOfBounds,
};

pub const Flags = packed struct {
    order: Order = .rowMajor,
    storage: Storage = .dense,
    ownsData: bool = true,
    writeable: bool = true,
};

pub const Order = enum(u2) {
    rowMajor,
    columnMajor,
    other,
};

pub const Storage = enum(u3) {
    dense,
    strided,
    csr,
    csc,
    coo,
};

pub const Metadata = union(Storage) {
    dense: Dense,
    strided: Strided,
    csr: CSR,
    csc: CSC,
    coo: COO,

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

pub const Slice = struct {
    start: usize,
    stop: usize,
    step: isize,

    pub fn init(start: usize, stop: usize, step: isize) !Slice {
        if (step == 0 and start != stop) {
            return Error.ZeroDimension;
        }

        if (step > 0) {
            if (start >= stop) {
                return Error.SliceOutOfBounds;
            }
        } else if (step < 0) {
            if (start <= stop) {
                return Error.SliceOutOfBounds;
            }
        }

        return Slice{ .start = start, .stop = stop, .step = step };
    }

    pub fn len(self: Slice) usize {
        if (self.step == 0 or self.start == self.stop) {
            return 1;
        }

        if (self.step > 0) {
            return (self.stop - self.start + cast(usize, self.step, .{}) - 1) / cast(usize, self.step, .{});
        }

        return (self.start - self.stop + cast(usize, int.abs(self.step), .{}) - 1) / cast(usize, int.abs(self.step), .{});
    }
};
