const std = @import("std");

const opts = @import("options");

const types = @import("types.zig");
const scast = types.scast;
const cast = types.cast;
const validateContext = types.validateContext;
const getFieldOrDefault = types.getFieldOrDefault;

const ops = @import("ops.zig");

const int = @import("int.zig");
const float = @import("float.zig");

const dense = @import("array/dense.zig");
const strided = @import("array/strided.zig");

pub const Iterator = @import("array/iterators.zig").Iterator;
//pub const MultiIterator = @import("array/iterators.zig").MultiIterator;
pub const IterationOrder = @import("array/iterators.zig").IterationOrder;

pub const max_dimensions = opts.max_dimensions;

pub fn Array(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Array requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: []T,
        ndim: usize,
        shape: [max_dimensions]usize,
        size: usize,
        base: ?*const Array(T),
        flags: Flags,
        metadata: Metadata,

        pub const empty: Array(T) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** max_dimensions,
            .size = 0,
            .base = null,
            .flags = .{},
            .metadata = .{ .dense = .{
                .strides = .{0} ** max_dimensions,
            } },
        };

        pub fn init(
            allocator: std.mem.Allocator,
            shape: []const usize,
            options: struct {
                order: Order = .col_major,
                storage: Storage = .dense,
            },
        ) !Array(T) {
            if (shape.len > max_dimensions) {
                return Error.TooManyDimensions;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return Error.ZeroDimension;
                }
            }

            switch (options.storage) {
                .dense => return dense.init(T, allocator, shape, options.order),
                .strided => return Error.InvalidFlags,
            }
        }

        pub fn full(
            allocator: std.mem.Allocator,
            shape: []const usize,
            value: anytype,
            options: struct {
                order: Order = .col_major,
                storage: Storage = .dense,
            },
            ctx: anytype,
        ) !Array(T) {
            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            if (shape.len > max_dimensions) {
                return Error.TooManyDimensions;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return Error.ZeroDimension;
                }
            }

            switch (options.storage) {
                .dense => return dense.full(T, allocator, shape, value, options.order, ctx),
                .strided => return Error.InvalidFlags,
            }
        }

        pub fn arange(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            step: anytype,
            ctx: anytype,
        ) !Array(T) {
            comptime if (types.isComplex(T))
                @compileError("array.arange does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.arange: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.arange: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.arange: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            comptime if (!types.isNumeric(@TypeOf(step)))
                @compileError("array.arange: step must be numeric, got " ++ @typeName(@TypeOf(step)));

            comptime if (types.isComplex(@TypeOf(step)))
                @compileError("array.arange: step cannot be complex, got " ++ @typeName(@TypeOf(step)));

            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            return dense.arange(T, allocator, start, stop, step, ctx);
        }

        pub fn linspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            num: usize,
            endpoint: bool,
            retstep: ?*T,
            ctx: anytype,
        ) !Array(T) {
            comptime if (types.isComplex(T))
                @compileError("array.linspace does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.linspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.linspace: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.linspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (num == 0) {
                return Error.ZeroDimension;
            }

            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            return dense.linspace(T, allocator, start, stop, num, endpoint, retstep, ctx);
        }

        pub fn logspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            num: usize,
            base: anytype,
            endpoint: bool,
            ctx: anytype,
        ) !Array(T) {
            comptime if (types.isComplex(T))
                @compileError("array.logspace does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.logspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.logspace: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.logspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (num == 0) {
                return Error.ZeroDimension;
            }

            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            return dense.logspace(T, allocator, start, stop, num, base, endpoint, ctx);
        }

        /// Cleans up the array by deinitializing its elements. If the array holds
        /// a fixed precision type this is does not do anything.
        pub fn cleanup(self: *Array(T), ctx: anytype) void {
            comptime if (types.isArbitraryPrecision(T)) {
                validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                validateContext(@TypeOf(ctx), .{});
            };

            if (comptime types.isArbitraryPrecision(T)) {
                switch (self.flags.storage) {
                    .dense => dense.cleanup(T, self.data, ctx),
                    .strided => strided.cleanup(T, self, self.size, self.flags.order.toIterationOrder(), ctx),
                }
            }
        }

        /// Deinitializes the array, freeing its data if it owns it.
        ///
        /// If the array holds an arbitrary precision type, it will not free the
        /// elements; call `cleanup` first to do that.
        pub fn deinit(self: *Array(T), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data);
            }

            self.* = undefined;
        }

        pub fn set(self: *Array(T), position: []const usize, value: T) !void {
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
            if (shape.len > max_dimensions) {
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
            const axes_: []const usize =
                axes orelse
                trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axes_.len == 0) {
                return Error.ZeroDimension;
            }

            if (axes_.len > self.ndim) {
                return Error.TooManyDimensions;
            }

            if (!isPermutation(self.ndim, axes_)) {
                return Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]
            }

            switch (self.flags.storage) {
                .dense => return dense.transpose(T, self, axes_),
                .strided => return strided.transpose(T, self, axes_),
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
            if (shape.len > max_dimensions)
                return Error.TooManyDimensions;

            if (shape.len < self.ndim) {
                return Error.TooLittleDimensions;
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
pub const apply2 = arrops.apply2;
pub const apply2_ = arrops.apply2_;

pub const add = arrops.add;
pub const add_ = arrops.add_;
pub const sub = arrops.sub;
pub const sub_ = arrops.sub_;
pub const mul = arrops.mul;
pub const mul_ = arrops.mul_;
pub const div = arrops.div;
pub const div_ = arrops.div_;

pub const eq = arrops.eq;
pub const eq_ = arrops.eq_;
pub const ne = arrops.ne;
pub const ne_ = arrops.ne_;
pub const lt = arrops.lt;
pub const lt_ = arrops.lt_;
pub const le = arrops.le;
pub const le_ = arrops.le_;
pub const gt = arrops.gt;
pub const gt_ = arrops.gt_;
pub const ge = arrops.ge;
pub const ge_ = arrops.ge_;

pub const max = arrops.max;
pub const max_ = arrops.max_;
pub const min = arrops.min;
pub const min_ = arrops.min_;

// Basic operations
pub const abs = arrops.abs;
pub const abs_ = arrops.abs_;
pub const abs2 = arrops.abs2;
pub const abs2_ = arrops.abs2_;

// Exponential functions
pub const exp = arrops.exp;
pub const exp_ = arrops.exp_;
pub const exp10 = arrops.exp10;
pub const exp10_ = arrops.exp10_;
pub const exp2 = arrops.exp2;
pub const exp2_ = arrops.exp2_;
pub const exp10m1 = arrops.exp10m1;
pub const exp10m1_ = arrops.exp10m1_;
pub const exp2m1 = arrops.exp2m1;
pub const exp2m1_ = arrops.exp2m1_;
pub const expm1 = arrops.expm1;
pub const expm1_ = arrops.expm1_;
pub const log = arrops.log;
pub const log_ = arrops.log_;
pub const log10 = arrops.log10;
pub const log10_ = arrops.log10_;
pub const log2 = arrops.log2;
pub const log2_ = arrops.log2_;
pub const log10p1 = arrops.log10p1;
pub const log10p1_ = arrops.log10p1_;
pub const log2p1 = arrops.log2p1;
pub const log2p1_ = arrops.log2p1_;
pub const log1p = arrops.log1p;
pub const log1p_ = arrops.log1p_;

// Power functions
pub const pow = arrops.pow;
pub const pow_ = arrops.pow_;
pub const sqrt = arrops.sqrt;
pub const sqrt_ = arrops.sqrt_;
pub const cbrt = arrops.cbrt;
pub const cbrt_ = arrops.cbrt_;
pub const hypot = arrops.hypot;
pub const hypot_ = arrops.hypot_;

// Trigonometric functions
pub const sin = arrops.sin;
pub const sin_ = arrops.sin_;
pub const cos = arrops.cos;
pub const cos_ = arrops.cos_;
pub const tan = arrops.tan;
pub const tan_ = arrops.tan_;
pub const asin = arrops.asin;
pub const asin_ = arrops.asin_;
pub const acos = arrops.acos;
pub const acos_ = arrops.acos_;
pub const atan = arrops.atan;
pub const atan_ = arrops.atan_;
pub const atan2 = arrops.atan2;
pub const atan2_ = arrops.atan2_;
pub const sinpi = arrops.sinpi;
pub const sinpi_ = arrops.sinpi_;
pub const cospi = arrops.cospi;
pub const cospi_ = arrops.cospi_;
pub const tanpi = arrops.tanpi;
pub const tanpi_ = arrops.tanpi_;
pub const asinpi = arrops.asinpi;
pub const asinpi_ = arrops.asinpi_;
pub const acospi = arrops.acospi;
pub const acospi_ = arrops.acospi_;
pub const atanpi = arrops.atanpi;
pub const atanpi_ = arrops.atanpi_;
pub const atan2pi = arrops.atan2pi;
pub const atan2pi_ = arrops.atan2pi_;

// Hyperbolic functions
pub const sinh = arrops.sinh;
pub const sinh_ = arrops.sinh_;
pub const cosh = arrops.cosh;
pub const cosh_ = arrops.cosh_;
pub const tanh = arrops.tanh;
pub const tanh_ = arrops.tanh_;
pub const asinh = arrops.asinh;
pub const asinh_ = arrops.asinh_;
pub const acosh = arrops.acosh;
pub const acosh_ = arrops.acosh_;
pub const atanh = arrops.atanh;
pub const atanh_ = arrops.atanh_;

// Error and gamma functions
pub const erf = arrops.erf;
pub const erf_ = arrops.erf_;
pub const erfc = arrops.erfc;
pub const erfc_ = arrops.erfc_;
pub const gamma = arrops.gamma;
pub const gamma_ = arrops.gamma_;
pub const lgamma = arrops.lgamma;
pub const lgamma_ = arrops.lgamma_;

// Nearest integer operations
pub const ceil = arrops.ceil;
pub const ceil_ = arrops.ceil_;

pub const Broadcast = struct {
    ndim: usize,
    shape: [max_dimensions]usize,
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

        if (shape.len > max_dimensions) {
            return Error.TooManyDimensions;
        }

        if (shape.len > ndim) {
            ndim = shape.len;
        }
    }

    var result: [max_dimensions]usize = .{0} ** max_dimensions;
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

    if (ndim == 0 or ndim > max_dimensions) {
        return false; // empty or too many dimensions is not a valid permutation
    }

    var seen: [max_dimensions]bool = .{false} ** max_dimensions;
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

pub fn trivialPermutation(ndim: usize) [max_dimensions]usize {
    var result: [max_dimensions]usize = .{0} ** max_dimensions;

    if (ndim == 0 or ndim > max_dimensions)
        return result; // empty or too many dimensions, return trivial permutation

    for (result[0..ndim], 0..) |*axis, i| {
        axis.* = i;
    }

    return result;
}

pub fn trivialReversePermutation(ndim: usize) [max_dimensions]usize {
    var result: [max_dimensions]usize = .{0} ** max_dimensions;

    if (ndim == 0 or ndim > max_dimensions)
        return result; // empty or too many dimensions, return empty permutation

    for (result[0..ndim], 0..) |*axis, i| {
        axis.* = ndim - i - 1;
    }

    return result;
}

pub const Error = error{
    ArrayNotWriteable,
    TooManyDimensions,
    TooLittleDimensions,
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
    order: Order = .col_major,
    storage: Storage = .dense,
    owns_data: bool = true,
};

pub const Order = enum(u1) {
    row_major,
    col_major,

    pub fn toIterationOrder(self: Order) IterationOrder {
        return switch (self) {
            .row_major => .right_to_left,
            .col_major => .left_to_right,
        };
    }

    pub fn resolve2(self: Order, other: Order) Order {
        if (self == other) {
            return self;
        }

        return .col_major; // default order
    }

    pub fn resolve3(self: Order, other1: Order, other2: Order) Order {
        if (self == other1 and self == other2)
            return self;

        if (self == other1 or self == other2)
            return self;

        if (other1 == other2)
            return other1;

        return .col_major; // default order
    }
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
    //sparse: Sparse,

    pub const Dense = struct {
        strides: [max_dimensions]usize,
    };

    pub const Strided = struct {
        strides: [max_dimensions]isize,
        offset: usize,
    };

    pub const Sparse = struct {
        // data holds the data, with the first element (index 0) being the
        // default value, and the rest being the nonzero elements.
        //
        // number of nonzero elements
        nnz: usize,
        // todo: explained in TODO.md
    };
};

pub const Range = struct {
    start: usize,
    stop: usize,
    step: isize,

    pub const all: Range = .{ .start = 0, .stop = int.maxVal(usize), .step = 1 };

    pub const all_reverse: Range = .{ .start = int.maxVal(usize), .stop = int.maxVal(usize), .step = -1 };

    pub fn init(start: ?usize, stop: ?usize, step: ?isize) !Range {
        const range: Range = .{
            .start = start orelse int.maxVal(usize),
            .stop = stop orelse int.maxVal(usize),
            .step = step orelse 1,
        };

        if (step == 0) {
            return Error.ZeroStep;
        }

        if (((range.step > 0 and range.start >= range.stop) or
            (range.step < 0 and range.start <= range.stop)) and
            (range.start != int.maxVal(usize) and range.stop != int.maxVal(usize)))
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
