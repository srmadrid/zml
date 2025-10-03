//! Storage scheme:
//!
//! If order is column major, BSC (Blocked Sparse Column), otherwise, i.e.,
//! row major, BSR (Blocked Sparse Row).

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

pub fn Block(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        flags: Flags = .{},

        pub const empty = Block(T, order){
            .data = &.{},
            .flags = .{ .owns_data = false },
        };
    };
}
