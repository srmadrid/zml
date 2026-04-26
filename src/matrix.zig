//! Namespace for matrix types and operations.

pub const general = @import("matrix/general.zig");
pub const symmetric = @import("matrix/symmetric.zig");
pub const hermitian = @import("matrix/hermitian.zig");
pub const triangular = @import("matrix/triangular.zig");
const diagonal = @import("matrix/diagonal.zig");
pub const Diagonal = diagonal.Diagonal;
const permutation = @import("matrix/permutation.zig");
pub const Permutation = permutation.Permutation;

pub const builder = @import("matrix/builder.zig");

const matops = @import("matrix/ops.zig");
pub const apply2_ = matops.apply2_;
// pub const Add = matops.Add;
// pub const add = matops.add;
pub const add_ = matops.add_;
// pub const Sub = matops.Sub;
// pub const sub = matops.sub;
pub const sub_ = matops.sub_;
// pub const Mul = matops.Mul;
// pub const mul = matops.mul;
pub const mul_ = matops.mul_;
// pub const Div = matops.Div;
// pub const div = matops.div;
pub const div_ = matops.div_;

pub const Error = error{
    ZeroDimension,
    PositionOutOfBounds,
    BreaksStructure,
    InvalidRange,
    DimensionMismatch,
    InvalidBandwidth,
    NotSquare,
    DataNotOwned,
    InsuficientSpace,
};

pub const Flags = packed struct {
    owns_data: bool = true,
};
