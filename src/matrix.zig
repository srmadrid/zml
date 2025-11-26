//! Namespace for matrix types and operations.

pub const general = @import("matrix/general.zig");
pub const symmetric = @import("matrix/symmetric.zig");
pub const hermitian = @import("matrix/hermitian.zig");
pub const triangular = @import("matrix/triangular.zig");
const diagonal = @import("matrix/diagonal.zig");
pub const Diagonal = diagonal.Diagonal;
const banded = @import("matrix/banded.zig");
pub const Banded = banded.Banded;
const tridiagonal = @import("matrix/tridiagonal.zig");
pub const Tridiagonal = tridiagonal.Tridiagonal;
const permutation = @import("matrix/permutation.zig");
pub const Permutation = permutation.Permutation;

pub const builder = @import("matrix/builder.zig");

const matops = @import("matrix/ops.zig");
pub const apply2 = matops.apply2;

pub const add = matops.add;
pub const sub = matops.sub;
pub const mul = matops.mul;
pub const div = matops.div;

pub const Error = error{
    ZeroDimension,
    PositionOutOfBounds,
    BreaksStructure,
    InvalidRange,
    DimensionMismatch,
    InvalidBandwidth,
    NotSquare,
    DataNotOwned,
};

pub const Flags = packed struct {
    owns_data: bool = true,
};
