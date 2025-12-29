//! Storage scheme:
//!
//! If order is column major, BSC (Blocked Block Column), otherwise, i.e.,
//! row major, BSR (Blocked Block Row).

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

pub fn Block(T: type, border: Order, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        idx: [*]u32,
        ptr: [*]u32,
        nnzb: u32,
        bsize: u32,
        rows: u32,
        cols: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_matrix = {};
        pub const is_block = {};
        pub const is_general = {};

        /// Numeric type
        pub const Numeric = T;

        pub const empty = Block(T, border, order){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnzb = 0,
            .bsize = 0,
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn deinit(self: *Block(T, border, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.nnzb * self.bsize * self.bsize]);
                allocator.free(self.idx[0..self.nnzb]);
                allocator.free(self.ptr[0..if (comptime order == .col_major) self.cols / self.bsize + 1 else self.rows / self.bsize + 1]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Block(T, border, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            if (comptime order == .col_major) {
                const col_start = self.ptr[bc];
                const col_end = self.ptr[bc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == br)
                        return self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj]
                    else if (self.idx[i] > br)
                        break;
                }
            } else {
                const row_start = self.ptr[br];
                const row_end = self.ptr[br + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == bc)
                        return self.data[j * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj]
                    else if (self.idx[j] > bc)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn at(self: *Block(T, border, order), r: u32, c: u32) T {
            // Unchecked version of get. Assumes r and c are valid.
            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            if (comptime order == .col_major) {
                const col_start = self.ptr[bc];
                const col_end = self.ptr[bc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == br)
                        return self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj]
                    else if (self.idx[i] > br)
                        break;
                }
            } else {
                const row_start = self.ptr[br];
                const row_end = self.ptr[br + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == bc)
                        return self.data[j * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj]
                    else if (self.idx[j] > bc)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn set(self: *Block(T, border, order), r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            // Find the position to update. If the position does not exist,
            // we return an error.
            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            if (comptime order == .col_major) {
                const col_start = self.ptr[bc];
                const col_end = self.ptr[bc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == br) {
                        self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
                        return;
                    } else if (self.idx[i] > br) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[br];
                const row_end = self.ptr[br + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == bc) {
                        self.data[j * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
                        return;
                    } else if (self.idx[j] > bc) {
                        break;
                    }
                }
            }

            return matrix.Error.BreaksStructure;
        }

        pub fn put(self: *Block(T, border, order), r: u32, c: u32, value: T) void {
            // Unchecked version of set. Assumes r and c are valid (i.e., within
            // bounds) and returns void if the position does not exist.
            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            if (comptime order == .col_major) {
                const col_start = self.ptr[bc];
                const col_end = self.ptr[bc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == br) {
                        self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
                        return;
                    } else if (self.idx[i] > br) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[br];
                const row_end = self.ptr[br + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == bc) {
                        self.data[j * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
                        return;
                    } else if (self.idx[j] > bc) {
                        break;
                    }
                }
            }

            return;
        }
    };
}
