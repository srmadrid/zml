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
pub fn Tridiagonal(T: type) type {
    if (!types.isNumeric(T))
        @compileError("matrix.Tridiagonal requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        size: u32,
        osize: u32,
        offset: u32,
        sdoffset: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_matrix = {};
        pub const is_tridiagonal = {};
        pub const is_symmetric = {};

        /// Numeric type
        pub const Numeric = T;
    };
}
