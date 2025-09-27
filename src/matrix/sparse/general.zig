//! Storage scheme:
//!
//! If order is column major, CSC (Compressed Sparse Column), otherwise, i.e.,
//! row major, CSR (Compressed Sparse Row).

const std = @import("std");

const types = @import("../../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Order = types.Order;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn General(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        idx: [*]u32,
        ptr: [*]u32,
        nnz: u32,
        rows: u32,
        cols: u32,
        flags: Flags = .{},

        pub const empty = General(T, order){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };
    };
}
