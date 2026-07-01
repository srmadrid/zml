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
pub const addUnchecked = vecops.addUnchecked;
pub const addAlloc = vecops.addAlloc;
pub const addInto = vecops.addInto;
pub const addIntoUnchecked = vecops.addIntoUnchecked;
pub const Sub = vecops.Sub;
pub const sub = vecops.sub;
pub const subUnchecked = vecops.subUnchecked;
pub const subAlloc = vecops.subAlloc;
pub const subInto = vecops.subInto;
pub const subIntoUnchecked = vecops.subIntoUnchecked;
pub const Mul = vecops.Mul;
pub const mul = vecops.mul;
pub const mulAlloc = vecops.mulAlloc;
pub const mulInto = vecops.mulInto;
pub const mulIntoUnchecked = vecops.mulIntoUnchecked;
pub const Div = vecops.Div;
pub const div = vecops.div;
pub const divAlloc = vecops.divAlloc;
pub const divInto = vecops.divInto;
pub const divIntoUnchecked = vecops.divIntoUnchecked;

pub const Error = error{
    ZeroLength,
    PositionOutOfBounds,
    DimensionMismatch,
    NonContiguousData,
    ZeroDimension,
    DataNotOwned,
    InsufficientSpace,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    noconj: bool = true,
};
