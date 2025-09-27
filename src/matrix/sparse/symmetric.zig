//! Storage scheme:
//!
//! If order is column major, CSC (Compressed Sparse Column), otherwise, i.e.,
//! row major, CSR (Compressed Sparse Row), storing only the upper or lower
//! triangular part of the matrix.

const std = @import("std");

const types = @import("../../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn Symmetric(T: type, uplo: Uplo, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        idx: [*]u32,
        ptr: [*]u32,
        nnz: u32,
        size: u32,
        _dlen: u32, // allocated length of data
        _ilen: u32, // allocated length of idx
        _plen: u32, // allocated length of ptr
        flags: Flags = .{},

        pub const empty = Symmetric(T, uplo, order){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .size = 0,
            ._dlen = 0,
            ._ilen = 0,
            ._plen = 0,
            .flags = .{ .owns_data = false },
        };
    };
}
