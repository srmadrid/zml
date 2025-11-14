//! Storage scheme:
//!
//! Block-COO (Coordinate List) format. Order chooses ordering of (row, col, data)
//! arrays:
//! - row_major: (row, col, data) sorted by row, then by col -> compiles to CSR
//! - col_major: (col, row, data) sorted by col, then by row -> compiles to CSC

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

pub fn Block(T: type, border: Order, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        row: [*]u32,
        col: [*]u32,
        nnzb: u32,
        bsize: u32,
        rows: u32,
        cols: u32,
        _dlen: u32, // allocated length of data
        _rlen: u32, // allocated length of rows
        _clen: u32, // allocated length of cols
        flags: Flags = .{},

        pub const empty = Block(T, border, order){
            .data = &.{},
            .row = &.{},
            .col = &.{},
            .nnzb = 0,
            .bsize = 0,
            .rows = 0,
            .cols = 0,
            ._dlen = 0,
            ._rlen = 0,
            ._clen = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(allocator: std.mem.Allocator, bsize: u32, rows: u32, cols: u32, nnzb: u32) !Block(T, border, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (rows % bsize != 0 or cols % bsize != 0 or nnzb == 0 or nnzb > (rows / bsize) * (cols / bsize))
                return matrix.Error.DimensionMismatch;

            const data: []T = try allocator.alloc(T, nnzb * bsize * bsize);
            errdefer allocator.free(data);

            const row: []u32 = try allocator.alloc(u32, nnzb);
            errdefer allocator.free(row);

            return .{
                .data = data.ptr,
                .row = row.ptr,
                .col = (try allocator.alloc(u32, nnzb)).ptr,
                .nnzb = 0,
                .bsize = bsize,
                .rows = rows,
                .cols = cols,
                ._dlen = nnzb * bsize * bsize,
                ._rlen = nnzb,
                ._clen = nnzb,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn deinit(self: *Block(T, border, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.row[0..self._rlen]);
                allocator.free(self.col[0..self._clen]);
            }

            self.* = undefined;
        }

        pub fn reserve(self: *Block(T, border, order), allocator: std.mem.Allocator, new_nnzb: u32) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnzb <= self._dlen / (self.bsize * self.bsize) and new_nnzb <= self._rlen and new_nnzb <= self._clen)
                return;

            if (new_nnzb > (self.rows / self.bsize) * (self.cols / self.bsize))
                return matrix.Error.DimensionMismatch;

            if (new_nnzb > self._dlen / (self.bsize * self.bsize)) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnzb * self.bsize * self.bsize)).ptr;
                self._dlen = new_nnzb * self.bsize * self.bsize;
            }

            if (new_nnzb > self._rlen) {
                self.row = (try allocator.realloc(self.row[0..self._rlen], new_nnzb)).ptr;
                self._rlen = new_nnzb;
            }

            if (new_nnzb > self._clen) {
                self.col = (try allocator.realloc(self.col[0..self._clen], new_nnzb)).ptr;
                self._clen = new_nnzb;
            }
        }

        pub fn get(self: *const Block(T, border, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc)
                    return self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj];

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
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

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc)
                    return self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj];

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn getBlock(self: *const Block(T, border, order), br: u32, bc: u32) !matrix.general.Dense(T, border) {
            if (br >= self.rows / self.bsize or bc >= self.cols / self.bsize)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc)
                    return .{
                        .data = self.data + i * self.bsize * self.bsize,
                        .rows = self.bsize,
                        .cols = self.bsize,
                        .ld = self.bsize,
                        .flags = .{ .owns_data = false },
                    };

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            return matrix.Error.PositionOutOfBounds;
        }

        pub fn set(self: *Block(T, border, order), allocator: std.mem.Allocator, r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc) {
                    self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            if (self.nnzb * self.bsize * self.bsize == self._dlen or self.nnzb == self._rlen or self.nnzb == self._clen) {
                // Need more space
                var new_nnzb = if (self.nnzb * self.bsize * self.bsize * 2 > self.rows * self.cols)
                    self.rows * self.cols / (self.bsize * self.bsize)
                else
                    self.nnzb * 2;

                if (new_nnzb == 0)
                    new_nnzb = 2;

                try self.reserve(allocator, new_nnzb);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnzb;
            while (j > i) : (j -= 1) {
                const bj: u32 = j * self.bsize * self.bsize;
                var k: u32 = 0;
                while (k < self.bsize * self.bsize) : (k += 1) {
                    self.data[bj + k] = self.data[bj - self.bsize * self.bsize + k];
                }

                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
            self.row[i] = br;
            self.col[i] = bc;
            self.nnzb += 1;
        }

        pub fn put(self: *Block(T, border, order), r: u32, c: u32, value: T) void {
            // Unchecked version of set. Assumes r and c are valid and there is space.
            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc) {
                    self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnzb;
            while (j > i) : (j -= 1) {
                const bj: u32 = j * self.bsize * self.bsize;
                var k: u32 = 0;
                while (k < self.bsize * self.bsize) : (k += 1) {
                    self.data[bj + k] = self.data[bj - self.bsize * self.bsize + k];
                }

                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
            self.row[i] = br;
            self.col[i] = bc;
            self.nnzb += 1;
        }

        pub fn setBlock(self: *Block(T, border, order), allocator: std.mem.Allocator, br: u32, bc: u32, block: matrix.general.Dense(T, border)) !void {
            if (br >= self.rows / self.bsize or bc >= self.cols / self.bsize)
                return matrix.Error.PositionOutOfBounds;

            if (block.rows != self.bsize or block.cols != self.bsize)
                return matrix.Error.DimensionMismatch;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc) {
                    if (comptime border == .col_major) {
                        var bj: u32 = 0;
                        while (bj < self.bsize) : (bj += 1) {
                            var bi: u32 = 0;
                            while (bi < self.bsize) : (bi += 1) {
                                self.data[i * self.bsize * self.bsize + bi + bj * self.bsize] = block.data[bi + bj * block.ld];
                            }
                        }
                    } else {
                        var bi: u32 = 0;
                        while (bi < self.bsize) : (bi += 1) {
                            var bj: u32 = 0;
                            while (bj < self.bsize) : (bj += 1) {
                                self.data[i * self.bsize * self.bsize + bi * self.bsize + bj] = block.data[bi * block.ld + bj];
                            }
                        }
                    }

                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            if (self.nnzb * self.bsize * self.bsize == self._dlen or self.nnzb == self._rlen or self.nnzb == self._clen) {
                // Need more space
                var new_nnzb = if (self.nnzb * self.bsize * self.bsize * 2 > self.rows * self.cols)
                    self.rows * self.cols / (self.bsize * self.bsize)
                else
                    self.nnzb * 2;

                if (new_nnzb == 0)
                    new_nnzb = 2;

                try self.reserve(allocator, new_nnzb);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnzb;
            while (j > i) : (j -= 1) {
                const bj: u32 = j * self.bsize * self.bsize;
                var k: u32 = 0;
                while (k < self.bsize * self.bsize) : (k += 1) {
                    self.data[bj + k] = self.data[bj - self.bsize * self.bsize + k];
                }

                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            if (comptime border == .col_major) {
                var bj: u32 = 0;
                while (bj < self.bsize) : (bj += 1) {
                    var bi: u32 = 0;
                    while (bi < self.bsize) : (bi += 1) {
                        self.data[i * self.bsize * self.bsize + bi + bj * self.bsize] = block.data[bi + bj * block.ld];
                    }
                }
            } else {
                var bi: u32 = 0;
                while (bi < self.bsize) : (bi += 1) {
                    var bj: u32 = 0;
                    while (bj < self.bsize) : (bj += 1) {
                        self.data[i * self.bsize * self.bsize + bi * self.bsize + bj] = block.data[bi * block.ld + bj];
                    }
                }
            }
            self.row[i] = br;
            self.col[i] = bc;
            self.nnzb += 1;
        }

        pub fn accumulate(self: *Block(T, border, order), allocator: std.mem.Allocator, r: u32, c: u32, value: anytype, ctx: anytype) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            const br: u32 = r / self.bsize;
            const bc: u32 = c / self.bsize;
            const ri: u32 = r % self.bsize;
            const cj: u32 = c % self.bsize;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc) {
                    try ops.add_(
                        &self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj],
                        self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj],
                        value,
                        ctx,
                    );

                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            if (self.nnzb * self.bsize * self.bsize == self._dlen or self.nnzb == self._rlen or self.nnzb == self._clen) {
                // Need more space
                var new_nnzb = if (self.nnzb * self.bsize * self.bsize * 2 > self.rows * self.cols)
                    self.rows * self.cols / (self.bsize * self.bsize)
                else
                    self.nnzb * 2;

                if (new_nnzb == 0)
                    new_nnzb = 2;

                try self.reserve(allocator, new_nnzb);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnzb;
            while (j > i) : (j -= 1) {
                const bj: u32 = j * self.bsize * self.bsize;
                var k: u32 = 0;
                while (k < self.bsize * self.bsize) : (k += 1) {
                    self.data[bj + k] = self.data[bj - self.bsize * self.bsize + k];
                }

                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i * self.bsize * self.bsize + if (comptime border == .col_major) ri + cj * self.bsize else ri * self.bsize + cj] = value;
            self.row[i] = br;
            self.col[i] = bc;
            self.nnzb += 1;
        }

        pub fn accumulateBlock(self: *Block(T, border, order), allocator: std.mem.Allocator, br: u32, bc: u32, block: matrix.general.Dense(T, border), ctx: anytype) !void {
            if (br >= self.rows / self.bsize or bc >= self.cols / self.bsize)
                return matrix.Error.PositionOutOfBounds;

            if (block.rows != self.bsize or block.cols != self.bsize)
                return matrix.Error.DimensionMismatch;

            var i: u32 = 0;
            while (i < self.nnzb) : (i += 1) {
                if (self.row[i] == br and self.col[i] == bc) {
                    if (comptime border == .col_major) {
                        var bj: u32 = 0;
                        while (bj < self.bsize) : (bj += 1) {
                            var bi: u32 = 0;
                            while (bi < self.bsize) : (bi += 1) {
                                try ops.add_(
                                    &self.data[i * self.bsize * self.bsize + bi + bj * self.bsize],
                                    self.data[i * self.bsize * self.bsize + bi + bj * self.bsize],
                                    block.data[bi + bj * block.ld],
                                    ctx,
                                );
                            }
                        }
                    } else {
                        var bi: u32 = 0;
                        while (bi < self.bsize) : (bi += 1) {
                            var bj: u32 = 0;
                            while (bj < self.bsize) : (bj += 1) {
                                try ops.add_(
                                    &self.data[i * self.bsize * self.bsize + bi * self.bsize + bj],
                                    self.data[i * self.bsize * self.bsize + bi * self.bsize + bj],
                                    block.data[bi * block.ld + bj],
                                    ctx,
                                );
                            }
                        }
                    }

                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > bc or (self.col[i] == bc and self.row[i] > br))
                        break;
                } else {
                    if (self.row[i] > br or (self.row[i] == br and self.col[i] > bc))
                        break;
                }
            }

            if (self.nnzb * self.bsize * self.bsize == self._dlen or self.nnzb == self._rlen or self.nnzb == self._clen) {
                // Need more space
                var new_nnzb = if (self.nnzb * self.bsize * self.bsize * 2 > self.rows * self.cols)
                    self.rows * self.cols / (self.bsize * self.bsize)
                else
                    self.nnzb * 2;
                if (new_nnzb == 0)
                    new_nnzb = 2;

                try self.reserve(allocator, new_nnzb);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnzb;
            while (j > i) : (j -= 1) {
                const bj: u32 = j * self.bsize * self.bsize;
                var k: u32 = 0;
                while (k < self.bsize * self.bsize) : (k += 1) {
                    self.data[bj + k] = self.data[bj - self.bsize * self.bsize + k];
                }

                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            if (comptime border == .col_major) {
                var bj: u32 = 0;
                while (bj < self.bsize) : (bj += 1) {
                    var bi: u32 = 0;
                    while (bi < self.bsize) : (bi += 1) {
                        self.data[i * self.bsize * self.bsize + bi + bj * self.bsize] = block.data[bi + bj * block.ld];
                    }
                }
            } else {
                var bi: u32 = 0;
                while (bi < self.bsize) : (bi += 1) {
                    var bj: u32 = 0;
                    while (bj < self.bsize) : (bj += 1) {
                        self.data[i * self.bsize * self.bsize + bi * self.bsize + bj] = block.data[bi * block.ld + bj];
                    }
                }
            }
            self.row[i] = br;
            self.col[i] = bc;
            self.nnzb += 1;
        }

        pub fn compile(self: *Block(T, border, order), allocator: std.mem.Allocator) !matrix.general.Block(T, border, order) {
            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols / self.bsize + 1 else self.rows / self.bsize + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnzb and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnzb and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnzb * self.bsize * self.bsize)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnzb * self.bsize * self.bsize)).ptr;

                if (self._rlen > self.nnzb)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnzb)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnzb * self.bsize * self.bsize)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnzb * self.bsize * self.bsize)).ptr;

                if (self._clen > self.nnzb)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnzb)).ptr;
            }

            const result = matrix.general.Block(T, border, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnzb = self.nnzb,
                .bsize = self.bsize,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        pub fn compileCopy(self: *Block(T, border, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.general.Block(T, border, order) {
            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols / self.bsize + 1 else self.rows / self.bsize + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnzb and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnzb and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            var data: []T = try allocator.alloc(T, self.nnzb * self.bsize * self.bsize);
            errdefer allocator.free(data);
            var idx: []u32 = try allocator.alloc(u32, self.nnzb);
            errdefer allocator.free(idx);

            i = 0;
            while (i < self.nnzb) : (i += 1) {
                data[i] = try ops.copy(self.data[i], ctx);
                idx[i] = if (comptime order == .col_major) self.row[i] else self.col[i];
            }

            while (i < self.nnzb * self.bsize * self.bsize) : (i += 1) {
                data[i] = try ops.copy(self.data[i], ctx);
            }

            const result = matrix.general.Block(T, border, order){
                .data = data.ptr,
                .idx = idx.ptr,
                .ptr = ptr.ptr,
                .nnzb = self.nnzb,
                .bsize = self.bsize,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }
    };
}
