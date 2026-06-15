//! Namespace for array types and operations.

const int = @import("int.zig");

const numeric = @import("numeric.zig");

pub const max_dimensions = @import("options").max_dimensions;

// Arrays shall no longer use Layout, use Order = enum { c, f }, but only on creation

const dense = @import("array/dense.zig");
pub const Dense = dense.Dense;
const strided = @import("array/strided.zig");
pub const Strided = strided.Strided;
const sparse = @import("array/sparse.zig");
pub const Sparse = sparse.Sparse;

// sum, prod, min, max, argmin, argmax (over whole array, or specific axes if array is static)
// _Alloc (over specific axes)
// _Into

// const arrops = @import("array/ops.zig");
// pub const apply1 = arrops.apply1;
// pub const apply1_ = arrops.apply1_;
// pub const apply2 = arrops.apply2;
// pub const apply2_ = arrops.apply2_;

// pub const add = arrops.add;
// pub const add_ = arrops.add_;
// pub const sub = arrops.sub;
// pub const sub_ = arrops.sub_;
// pub const mul = arrops.mul;
// pub const mul_ = arrops.mul_;
// pub const div = arrops.div;
// pub const div_ = arrops.div_;

// pub const eq = arrops.eq;
// pub const eq_ = arrops.eq_;
// pub const ne = arrops.ne;
// pub const ne_ = arrops.ne_;
// pub const lt = arrops.lt;
// pub const lt_ = arrops.lt_;
// pub const le = arrops.le;
// pub const le_ = arrops.le_;
// pub const gt = arrops.gt;
// pub const gt_ = arrops.gt_;
// pub const ge = arrops.ge;
// pub const ge_ = arrops.ge_;

// pub const max = arrops.max;
// pub const max_ = arrops.max_;
// pub const min = arrops.min;
// pub const min_ = arrops.min_;

// // Basic operations
// pub const abs = arrops.abs;
// pub const abs_ = arrops.abs_;
// pub const abs1 = arrops.abs1;
// pub const abs1_ = arrops.abs1_;
// pub const abs2 = arrops.abs2;
// pub const abs2_ = arrops.abs2_;

// // Exponential functions
// pub const exp = arrops.exp;
// pub const exp_ = arrops.exp_;
// pub const exp10 = arrops.exp10;
// pub const exp10_ = arrops.exp10_;
// pub const exp2 = arrops.exp2;
// pub const exp2_ = arrops.exp2_;
// pub const log = arrops.log;
// pub const log_ = arrops.log_;
// pub const log10 = arrops.log10;
// pub const log10_ = arrops.log10_;
// pub const log2 = arrops.log2;
// pub const log2_ = arrops.log2_;

// // Power functions
// pub const pow = arrops.pow;
// pub const pow_ = arrops.pow_;
// pub const sqrt = arrops.sqrt;
// pub const sqrt_ = arrops.sqrt_;
// pub const cbrt = arrops.cbrt;
// pub const cbrt_ = arrops.cbrt_;
// pub const hypot = arrops.hypot;
// pub const hypot_ = arrops.hypot_;

// // Trigonometric functions
// pub const sin = arrops.sin;
// pub const sin_ = arrops.sin_;
// pub const cos = arrops.cos;
// pub const cos_ = arrops.cos_;
// pub const tan = arrops.tan;
// pub const tan_ = arrops.tan_;
// pub const asin = arrops.asin;
// pub const asin_ = arrops.asin_;
// pub const acos = arrops.acos;
// pub const acos_ = arrops.acos_;
// pub const atan = arrops.atan;
// pub const atan_ = arrops.atan_;
// pub const atan2 = arrops.atan2;
// pub const atan2_ = arrops.atan2_;

// // Hyperbolic functions
// pub const sinh = arrops.sinh;
// pub const sinh_ = arrops.sinh_;
// pub const cosh = arrops.cosh;
// pub const cosh_ = arrops.cosh_;
// pub const tanh = arrops.tanh;
// pub const tanh_ = arrops.tanh_;
// pub const asinh = arrops.asinh;
// pub const asinh_ = arrops.asinh_;
// pub const acosh = arrops.acosh;
// pub const acosh_ = arrops.acosh_;
// pub const atanh = arrops.atanh;
// pub const atanh_ = arrops.atanh_;

// // Error and gamma functions
// pub const erf = arrops.erf;
// pub const erf_ = arrops.erf_;
// pub const erfc = arrops.erfc;
// pub const erfc_ = arrops.erfc_;
// pub const gamma = arrops.gamma;
// pub const gamma_ = arrops.gamma_;
// pub const lgamma = arrops.lgamma;
// pub const lgamma_ = arrops.lgamma_;

// // Nearest integer operations
// pub const ceil = arrops.ceil;
// pub const ceil_ = arrops.ceil_;

pub fn broadcastShapes(shapes: []const []const usize) !struct { ndim: usize, shape: [max_dimensions]usize } {
    if (shapes.len == 0)
        return Error.ZeroDimension;

    var ndim: usize = 0;
    for (shapes) |shape| {
        if (shape.len == 0)
            return Error.ZeroDimension;

        if (shape.len > max_dimensions)
            return Error.TooManyDimensions;

        if (shape.len > ndim)
            ndim = shape.len;
    }

    var result: [max_dimensions]usize = .{0} ** max_dimensions;
    var i: isize = numeric.cast(isize, ndim - 1);
    while (i >= 0) : (i -= 1) {
        var max_dim: usize = 1;
        for (shapes) |shape| {
            const diff: isize = numeric.cast(isize, ndim - shape.len);
            if (i - diff >= 0) {
                if (shape[numeric.cast(usize, i - diff)] == 0)
                    return Error.ZeroDimension;

                if (shape[numeric.cast(usize, i - diff)] > max_dim) {
                    if (max_dim != 1)
                        return Error.NotBroadcastable;

                    max_dim = shape[numeric.cast(usize, i - diff)];
                }

                if (shape[numeric.cast(usize, i - diff)] != 1 and
                    shape[numeric.cast(usize, i - diff)] != max_dim)
                    return Error.NotBroadcastable;
            }
        }

        result[numeric.cast(usize, i)] = max_dim;
    }

    return .{
        .ndim = ndim,
        .shape = result,
    };
}

/// Checks if the given axes form a valid permutation of `[0, ..., ndim - 1]`.
pub fn isPermutation(ndim: usize, axes: []const usize) bool {
    if (ndim != axes.len)
        return false; // axes must match the shape length

    if (ndim == 0 or ndim > max_dimensions)
        return false; // empty or too many dimensions is not a valid permutation

    var seen: [max_dimensions]bool = .{false} ** max_dimensions;
    for (axes) |axis| {
        if (axis >= ndim)
            return false; // axis out of bounds

        if (seen[axis])
            return false; // duplicate axis

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

    var i: usize = 0;
    for (result[0..ndim]) |*axis| {
        axis.* = ndim - i - 1;

        i += 1;
    }

    return result;
}

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
            return (self.stop - self.start + numeric.cast(usize, self.step) - 1) / numeric.cast(usize, self.step);
        }

        return (self.start - self.stop + numeric.cast(usize, int.abs(self.step)) - 1) / numeric.cast(usize, int.abs(self.step));
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

pub const Flags = packed struct {
    owns_data: bool = true,
};
