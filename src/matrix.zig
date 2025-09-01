const std = @import("std");

const opts = @import("options");

const types = @import("types.zig");

const ops = @import("ops.zig");

const int = @import("int.zig");

const general = @import("matrix/general.zig");
pub const General = general.General;
const symmetric = @import("matrix/symmetric.zig");
pub const Symmetric = symmetric.Symmetric;
const hermitian = @import("matrix/hermitian.zig");
pub const Hermitian = hermitian.Hermitian;
const triangular = @import("matrix/triangular.zig");
pub const Triangular = triangular.Triangular;
const diagonal = @import("matrix/diagonal.zig");
pub const Diagonal = diagonal.Diagonal;
const banded = @import("matrix/banded.zig");
pub const Banded = banded.Banded;
const tridiagonal = @import("matrix/tridiagonal.zig");
pub const Tridiagonal = tridiagonal.Tridiagonal;
const permutation = @import("matrix/permutation.zig");
pub const Permutation = permutation.Permutation;
const sparse = @import("matrix/sparse.zig");
pub const Sparse = sparse.Sparse;

const matops = @import("matrix/ops.zig");
pub const apply2 = matops.apply2;

pub const add = matops.add;
pub const sub = matops.sub;
pub const mul = matops.mul;
pub const div = matops.div;

pub const Error = error{
    ZeroDimension,
    PositionOutOfBounds,
    InvalidRange,
    BreaksStructure,
    DimensionMismatch,
    InvalidBandwidth,
};

pub const Flags = packed struct {
    owns_data: bool = true,
};
