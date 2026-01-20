const std = @import("std");

const types = @import("../types.zig");
const Scalar = types.Scalar;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const scast = types.scast;
const cast = types.cast;
const validateContext = types.validateContext;

const ops = @import("../ops.zig");
const int = @import("../int.zig");

const array = @import("../array.zig");
const max_dimensions = array.max_dimensions;
const Order = types.Layout;
const Flags = array.Flags;
const Range = array.Range;

const dense = @import("dense.zig");
const Dense = dense.Dense;

pub fn Sparse(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("Strided requires a numeric type, got " ++ @typeName(T));

    return struct {
        nnz: usize,

        /// Type signatures
        pub const is_array = {};
        pub const is_sparse = {};

        /// Numeric type
        pub const Numeric = T;

        pub const empty: Dense(T, order) = .{
            .nnz = 0,
        };
    };
}
