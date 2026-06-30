//! Namespace for array types and operations.

const std = @import("std");

const int = @import("int.zig");

const numeric = @import("numeric.zig");

pub const max_dimensions = @import("options").max_dimensions;

// Arrays shall no longer use Layout, use Order = enum { c, f }, but only on creation

const static = @import("array/static.zig");
pub const Static = static.Static;
const dense = @import("array/dense.zig");
pub const Dense = dense.Dense;
const sparse = @import("array/sparse.zig");
pub const Sparse = sparse.Sparse;

pub const builder = @import("array/builder.zig");

// sum, prod, min, max, argmin, argmax (over whole array, or specific axes if array is static)
// _Alloc (over specific axes)
// _Into

// const arrops = @import("array/ops.zig");
// pub const apply1 = arrops.apply1;
// pub const apply1Alloc = arrops.apply1Alloc;
// pub const apply1Into = arrops.apply1Into;
// pub const apply2 = arrops.apply2;
// pub const apply2Alloc = arrops.apply2Alloc;
// pub const apply2Into = arrops.apply2Into;

pub const Order = enum(u1) {
    c,
    f,

    pub const default: Order = .f;
};

pub const Range = struct {
    start: ?usize,
    stop: ?usize,
    step: isize,
    is_index: bool = false,

    /// Selects the entire dimension: `[:]`
    pub const all: Range = .{ .start = null, .stop = null, .step = 1 };

    /// Selects the entire dimension in reverse: `[::-1]`
    pub const all_reverse: Range = .{ .start = null, .stop = null, .step = -1 };

    /// Creates a slice range: `[start:stop:step]`
    pub fn slice(start: ?usize, stop: ?usize, step: ?isize) !Range {
        const s = step orelse 1;
        if (s == 0) return Error.ZeroStep;
        return .{ .start = start, .stop = stop, .step = s, .is_index = false };
    }

    /// Selects a single integer index, collapsing the dimension: `[index]`
    pub fn index(idx: usize) Range {
        return .{ .start = idx, .stop = idx + 1, .step = 1, .is_index = true };
    }

    /// Resolves concrete boundaries against a specific dimension size.
    pub fn resolve(self: Range, dim_size: usize) !struct { start: usize, stop: usize, len: usize } {
        if (dim_size == 0)
            return .{ .start = 0, .stop = 0, .len = 0 };

        var actual_start: usize = undefined;
        var actual_stop: usize = undefined;
        if (self.step > 0) {
            actual_start = int.min(self.start orelse 0, dim_size);
            actual_stop = int.min(self.stop orelse dim_size, dim_size);

            if (actual_start >= actual_stop)
                return .{ .start = actual_start, .stop = actual_stop, .len = 0 };

            const diff = actual_stop - actual_start;
            const ustep = numeric.cast(usize, self.step);
            return .{
                .start = actual_start,
                .stop = actual_stop,
                .len = try std.math.divCeil(usize, diff, ustep),
            };
        } else {
            actual_start = int.min(self.start orelse (dim_size - 1), dim_size - 1);
            const go_to_end = (self.stop == null);
            actual_stop = if (go_to_end) 0 else int.min(self.stop.?, dim_size - 1);

            if (!go_to_end and actual_start <= actual_stop)
                return .{ .start = actual_start, .stop = actual_stop, .len = 0 };

            const diff = if (go_to_end) (actual_start + 1) else (actual_start - actual_stop);
            const ustep = numeric.cast(usize, -self.step);
            return .{
                .start = actual_start,
                .stop = actual_stop,
                .len = try std.math.divCeil(usize, diff, ustep),
            };
        }
    }
};

pub const Error = error{
    ArrayNotWriteable,
    TooManyDimensions,
    TooLittleDimensions,
    InvalidShape,
    InvalidFlags,
    InvalidAxes,
    InvalidKind,
    IncompatibleLayout,
    ZeroDimension,
    NotImplemented,
    NotBroadcastable,
    NotConvertible,
    DimensionMismatch,
    PositionOutOfBounds,
    InvalidRange,
    RangeOutOfBounds,
    ZeroStep,
    NeedDense,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    noconj: bool = true,
};
