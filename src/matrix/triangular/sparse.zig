//! Storage scheme:
//!
//! If order is column major, CSC (Compressed Sparse Column), otherwise, i.e.,
//! row major, CSR (Compressed Sparse Row), storing only the upper or lower
//! triangular part of the matrix, with or without the diagonal.

const std = @import("std");

const types = @import("../../types.zig");
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn Sparse(T: type, uplo: Uplo, diag: Diag, order: Order) type {
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

        pub const empty = Sparse(T, uplo, diag, order){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn deinit(self: *Sparse(T, uplo, diag, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.nnz]);
                allocator.free(self.idx[0..self.nnz]);
                allocator.free(self.ptr[0..(if (comptime order == .col_major) self.cols + 1 else self.rows + 1)]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Sparse(T, uplo, diag, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return constants.zero(T, .{}) catch unreachable;
            } else {
                if (r < c)
                    return constants.zero(T, .{}) catch unreachable;
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return constants.one(T, .{}) catch unreachable;
            }

            if (comptime order == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r)
                        return self.data[i]
                    else if (self.idx[i] > r)
                        break;
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c)
                        return self.data[j]
                    else if (self.idx[j] > c)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn at(self: *Sparse(T, uplo, diag, order), r: u32, c: u32) T {
            // Unchecked version of get. Assumes r and c are valid and on the correct
            // triangular part, and outside the diagonal if diag is unit.
            if (comptime order == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r)
                        return self.data[i]
                    else if (self.idx[i] > r)
                        break;
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c)
                        return self.data[j]
                    else if (self.idx[j] > c)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn set(self: *Sparse(T, uplo, diag, order), r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return matrix.Error.BreaksStructure;
            } else {
                if (r < c)
                    return matrix.Error.BreaksStructure;
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return matrix.Error.BreaksStructure;
            }

            // Find the position to update. If the position does not exist,
            // we return an error.
            if (comptime order == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r) {
                        self.data[i] = value;
                        return;
                    } else if (self.idx[i] > r) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c) {
                        self.data[j] = value;
                        return;
                    } else if (self.idx[j] > c) {
                        break;
                    }
                }
            }

            return matrix.Error.BreaksStructure;
        }

        pub fn put(self: *Sparse(T, uplo, diag, order), r: u32, c: u32, value: T) void {
            // Unchecked version of set. Assumes r and c are valid (i.e., within
            // bounds) and on the correct triangular part, and outside the diagonal
            // if diag is unit. Returns void if the position does not exist.
            if (comptime order == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r) {
                        self.data[i] = value;
                        return;
                    } else if (self.idx[i] > r) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c) {
                        self.data[j] = value;
                        return;
                    } else if (self.idx[j] > c) {
                        break;
                    }
                }
            }

            return;
        }
    };
}
