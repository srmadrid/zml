const std = @import("std");

const types = @import("../../types.zig");
const Order = types.Order;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

///
pub fn Banded(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("matrix.Banded requires a numeric type, got " ++ @typeName(T));

    _ = order;

    return struct {
        data: [*]T,
        rows: u32,
        cols: u32,
        ld: u32, // leading dimension
        bandwidth: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_matrix = {};
        pub const is_banded = {};
        pub const is_symmetric = {};

        /// Numeric type
        pub const Numeric = T;
    };
}
