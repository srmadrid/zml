//! Storage scheme:
//!
//! If column-major order is used, CSC (Compressed Sparse Column) format is
//! used, if row-major order is used, CSR (Compressed Sparse Row) format is
//! used.

const std = @import("std");

const types = @import("../types.zig");
const Order = types.Order;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const matrix = @import("../matrix.zig");
const Flags = matrix.Flags;

pub fn Sparse(comptime T: type) type {
    if (!types.isNumeric(T))
        @compileError("Sparse requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        nnz: u32,

        pub const empty: Sparse(T) = Sparse(.{
            .data = &.{},
            .nnz = 0,
        });
    };
}
