const std = @import("std");

const opts = @import("options");

const types = @import("types.zig");
const Order = types.Order;

const ops = @import("ops.zig");

const int = @import("int.zig");
const float = @import("float.zig");

const dense = @import("array/dense.zig");
pub const Dense = dense.Dense;
pub const strided = @import("array/strided.zig");
pub const Strided = strided.Strided;
const sparse = @import("array/sparse.zig");
pub const Sparse = sparse.Sparse;

pub fn AnyArray(comptime T: type) type {
    return union(enum) {
        dense: *const Dense(T),
        strided: *const Strided(T),
        sparse: *const Sparse(T),
    };
}

pub const Iterator = @import("array/iterators.zig").Iterator;
//pub const MultiIterator = @import("array/iterators.zig").MultiIterator;

pub const max_dimensions = opts.max_dimensions;

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
    ndim: u32,
    shape: [max_dimensions]u32,
};

pub fn broadcastShapes(
    shapes: []const []const u32,
) !Broadcast {
    if (shapes.len == 0) {
        return Error.ZeroDimension;
    }

    var ndim: u32 = 0;
    for (shapes) |shape| {
        if (shape.len == 0) {
            return Error.ZeroDimension;
        }

        if (shape.len > max_dimensions) {
            return Error.TooManyDimensions;
        }

        if (shape.len > ndim) {
            ndim = types.scast(u32, shape.len);
        }
    }

    var result: [max_dimensions]u32 = .{0} ** max_dimensions;
    var i: i32 = types.scast(i32, ndim - 1);
    while (i >= 0) : (i -= 1) {
        var max_dim: u32 = 1;
        for (shapes) |shape| {
            const diff: i32 = types.scast(i32, ndim - shape.len);
            if (i - diff >= 0) {
                if (shape[types.scast(u32, i - diff)] == 0) {
                    return Error.ZeroDimension;
                }

                if (shape[types.scast(u32, i - diff)] > max_dim) {
                    if (max_dim != 1) {
                        return Error.NotBroadcastable;
                    }

                    max_dim = shape[types.scast(u32, i - diff)];
                }

                if (shape[types.scast(u32, i - diff)] != 1 and
                    shape[types.scast(u32, i - diff)] != max_dim)
                {
                    return Error.NotBroadcastable;
                }
            }
        }

        result[types.scast(u32, i)] = max_dim;
    }

    return .{
        .ndim = ndim,
        .shape = result,
    };
}

/// Checks if the given axes form a valid permutation of `[0, ..., ndim - 1]`.
pub fn isPermutation(
    ndim: u32,
    axes: []const u32,
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

pub fn trivialPermutation(ndim: u32) [max_dimensions]u32 {
    var result: [max_dimensions]u32 = .{0} ** max_dimensions;

    if (ndim == 0 or ndim > max_dimensions)
        return result; // empty or too many dimensions, return trivial permutation

    for (result[0..ndim], 0..) |*axis, i| {
        axis.* = i;
    }

    return result;
}

pub fn trivialReversePermutation(ndim: u32) [max_dimensions]u32 {
    var result: [max_dimensions]u32 = .{0} ** max_dimensions;

    if (ndim == 0 or ndim > max_dimensions)
        return result; // empty or too many dimensions, return empty permutation

    var i: u32 = 0;
    for (result[0..ndim]) |*axis| {
        axis.* = ndim - i - 1;

        i += 1;
    }

    return result;
}

pub const Range = struct {
    start: u32,
    stop: u32,
    step: i32,

    pub const all: Range = .{ .start = 0, .stop = int.maxVal(u32), .step = 1 };

    pub const all_reverse: Range = .{ .start = int.maxVal(u32), .stop = int.maxVal(u32), .step = -1 };

    pub fn init(start: ?u32, stop: ?u32, step: ?i32) !Range {
        const range: Range = .{
            .start = start orelse int.maxVal(u32),
            .stop = stop orelse int.maxVal(u32),
            .step = step orelse 1,
        };

        if (step == 0) {
            return Error.ZeroStep;
        }

        if (((range.step > 0 and range.start >= range.stop) or
            (range.step < 0 and range.start <= range.stop)) and
            (range.start != int.maxVal(u32) and range.stop != int.maxVal(u32)))
        {
            return Error.RangeOutOfBounds;
        }

        return range;
    }

    pub fn single(index: u32) Range {
        return Range{ .start = index, .stop = index + 1, .step = 1 };
    }

    pub fn len(self: Range) u32 {
        if (self.start == self.stop) {
            return 0;
        }

        if (self.step > 0) {
            return (self.stop - self.start + types.scast(u32, self.step) - 1) / types.scast(u32, self.step);
        }

        return (self.start - self.stop + types.scast(u32, int.abs(self.step)) - 1) / types.scast(u32, int.abs(self.step));
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

// Flags common between Dense, Strided, and Sparse arrays.
pub const Flags = packed struct {
    order: Order = .col_major,
    owns_data: bool = true,
};
