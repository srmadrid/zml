pub const dense = @import("matrix/dense.zig");
pub const sparse = @import("matrix/sparse.zig");

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
};

pub const Flags = packed struct {
    owns_data: bool = true,
};
