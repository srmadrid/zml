const std = @import("std");
const types = @import("types.zig");
const cast = types.cast;

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

        pub fn set(
            self: *NDArray(T),
            position: []const usize,
            value: anytype,
            options: struct {
                allocator: ?std.mem.Allocator = null,
            },
        ) !void {
            _ = types.numericType(@TypeOf(value));

            switch (self.flags.storage) {
                .dense => {
                    if (self.size == 1) {
                        self.data[0] = value;
                        return;
                    }

                    try dense.checkPosition(self, position);

                    self.data[dense.index(self, position)] = cast(T, value, .{ .allocator = options.allocator });
                },
                .strided => {
                    if (self.size == 1) {
                        self.data[self.metadata.strided.offset] = value;
                        return;
                    }

                    try strided.checkPosition(position);

                    self.data[strided.index(self, position)] = cast(T, value, .{ .allocator = options.allocator });
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

        pub fn get(self: *NDArray(T), position: []const usize) !*T {
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

                    try dense.checkPosition(T, self, position);

                    return &self.data[dense.index(T, self, position)];
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
