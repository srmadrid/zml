//! Namespace for vector types and operations.

const static = @import("vector/static.zig");
pub const Static = static.Static;
const dense = @import("vector/dense.zig");
pub const Dense = dense.Dense;
const sparse = @import("vector/sparse.zig");
pub const Sparse = sparse.Sparse;

const vecops = @import("vector/ops.zig");
pub const Add = vecops.Add;
pub const add = vecops.add;
pub const add_ = vecops.add_;
pub const Sub = vecops.Sub;
pub const sub = vecops.sub;
pub const sub_ = vecops.sub_;
pub const Mul = vecops.Mul;
pub const mul = vecops.mul;
pub const mul_ = vecops.mul_;
pub const Div = vecops.Div;
pub const div = vecops.div;
pub const div_ = vecops.div_;

pub const Error = error{
    ZeroLength,
    PositionOutOfBounds,
    DimensionMismatch,
    NonContiguousData,
    ZeroDimension,
    DataNotOwned,
};

pub const Flags = packed struct {
    owns_data: bool = true,
};
