//! Namespace for polynomial types and operations.

const static = @import("poly/static.zig");
pub const Static = static.Static;
const dense = @import("poly/dense.zig");
pub const Dense = dense.Dense;
const sparse = @import("poly/sparse.zig");
pub const Sparse = sparse.Sparse;

// const polops = @import("poly/ops.zig");
// pub const Add = polops.Add;
// pub const add = polops.add;
// pub const addUnchecked = polops.addUnchecked;
// pub const addAlloc = polops.addAlloc;
// pub const addInto = polops.addInto;
// pub const addIntoUnchecked = polops.addIntoUnchecked;
// pub const Sub = polops.Sub;
// pub const sub = polops.sub;
// pub const subUnchecked = polops.subUnchecked;
// pub const subAlloc = polops.subAlloc;
// pub const subInto = polops.subInto;
// pub const subIntoUnchecked = polops.subIntoUnchecked;
// pub const Mul = polops.Mul;
// pub const mul = polops.mul;
// pub const mulAlloc = polops.mulAlloc;
// pub const mulInto = polops.mulInto;
// pub const mulIntoUnchecked = polops.mulIntoUnchecked;
// pub const Div = polops.Div;
// pub const div = polops.div;
// pub const divAlloc = polops.divAlloc;
// pub const divInto = polops.divInto;
// pub const divIntoUnchecked = polops.divIntoUnchecked;

pub const Error = error{
    PositionOutOfBounds,
    DimensionMismatch,
    ZeroLength,
    DataNotOwned,
};

pub const Flags = packed struct {
    owns_data: bool = true,
};
