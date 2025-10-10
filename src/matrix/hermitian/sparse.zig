//! Storage scheme:
//!
//! If order is column major, CSC (Compressed Sparse Column), otherwise, i.e.,
//! row major, CSR (Compressed Sparse Row), storing only the upper or lower
//! triangular part of the matrix.
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
const Uplo = types.Uplo;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn Sparse(T: type, uplo: Uplo, order: Order) type {
    if (!types.isNumeric(T) or !types.isComplex(T))
        @compileError("T must be a complex numeric type");

    return struct {
        data: [*]T,
        idx: [*]u32,
        ptr: [*]u32,
        nnz: u32,
        size: u32,
        flags: Flags = .{},

        pub const empty = Sparse(T, uplo, order){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .size = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn deinit(self: *Sparse(T, uplo, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.nnz]);
                allocator.free(self.idx[0..self.nnz]);
                allocator.free(self.ptr[0 .. self.size + 1]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Sparse(T, uplo, order), r: u32, c: u32) !T {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var rr: u32 = r;
            var cc: u32 = c;
            var noconj: bool = true;
            if (comptime uplo == .upper) {
                if (rr > cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    noconj = false;
                }
            } else {
                if (rr < cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    noconj = false;
                }
            }

            if (comptime order == .col_major) {
                const col_start = self.ptr[cc];
                const col_end = self.ptr[cc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == rr)
                        return if (noconj)
                            self.data[i]
                        else
                            ops.conj(self.data[i], .{}) catch unreachable
                    else if (self.idx[i] > rr)
                        break;
                }
            } else {
                const row_start = self.ptr[rr];
                const row_end = self.ptr[rr + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == cc)
                        return if (noconj)
                            self.data[j]
                        else
                            ops.conj(self.data[j], .{}) catch unreachable
                    else if (self.idx[j] > cc)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn at(self: *Sparse(T, uplo, order), r: u32, c: u32) T {
            // Unchecked version of get. Assumes r and c are valid and on
            // the correct triangular part.
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

        pub fn set(self: *Sparse(T, uplo, order), r: u32, c: u32, value: T) !void {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var rr: u32 = r;
            var cc: u32 = c;
            var noconj: bool = true;
            if (comptime uplo == .upper) {
                if (rr > cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    noconj = false;
                }
            } else {
                if (rr < cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    noconj = false;
                }
            }

            // Find the position to update. If the position does not exist,
            // we return an error.
            if (comptime order == .col_major) {
                const col_start = self.ptr[cc];
                const col_end = self.ptr[cc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == rr) {
                        self.data[i] = if (noconj)
                            value
                        else
                            ops.conj(value, .{}) catch unreachable;
                        return;
                    } else if (self.idx[i] > rr) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[rr];
                const row_end = self.ptr[rr + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == cc) {
                        self.data[j] = if (noconj)
                            value
                        else
                            ops.conj(value, .{}) catch unreachable;
                        return;
                    } else if (self.idx[j] > cc) {
                        break;
                    }
                }
            }

            return matrix.Error.BreaksStructure;
        }

        pub fn put(self: *Sparse(T, uplo, order), r: u32, c: u32, value: T) void {
            // Unchecked version of set. Assumes r and c are valid (i.e., within
            // bounds) and returns void if the position does not exist.
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
