//! Storage scheme:
//!
//! COO (Coordinate List) format. Order chooses ordering of (row, col, data)
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
const Diag = types.Diag;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

const array = @import("../../array.zig");

pub fn Sparse(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        row: [*]u32,
        col: [*]u32,
        nnz: u32,
        rows: u32,
        cols: u32,
        _dlen: u32, // allocated length of data
        _rlen: u32, // allocated length of rows
        _clen: u32, // allocated length of cols
        flags: Flags = .{},

        pub const empty = Sparse(T, order){
            .data = &.{},
            .row = &.{},
            .col = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            ._dlen = 0,
            ._rlen = 0,
            ._clen = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(allocator: std.mem.Allocator, rows: u32, cols: u32, nnz: u32) !Sparse(T, order) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (nnz == 0 or nnz > rows * cols)
                return matrix.Error.DimensionMismatch;

            const data: []T = try allocator.alloc(T, nnz);
            errdefer allocator.free(data);

            const row: []u32 = try allocator.alloc(u32, nnz);
            errdefer allocator.free(row);

            return .{
                .data = data.ptr,
                .row = row.ptr,
                .col = (try allocator.alloc(u32, nnz)).ptr,
                .nnz = 0,
                .rows = rows,
                .cols = cols,
                ._dlen = nnz,
                ._rlen = nnz,
                ._clen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn deinit(self: *Sparse(T, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.row[0..self._rlen]);
                allocator.free(self.col[0..self._clen]);
            }

            self.* = undefined;
        }

        pub fn reserve(self: *Sparse(T, order), allocator: std.mem.Allocator, new_nnz: u32) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._rlen and new_nnz <= self._clen)
                return;

            if (new_nnz > self.rows * self.cols)
                return matrix.Error.DimensionMismatch;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._rlen) {
                self.row = (try allocator.realloc(self.row[0..self._rlen], new_nnz)).ptr;
                self._rlen = new_nnz;
            }

            if (new_nnz > self._clen) {
                self.col = (try allocator.realloc(self.col[0..self._clen], new_nnz)).ptr;
                self._clen = new_nnz;
            }
        }

        pub fn get(self: *const Sparse(T, order), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c)
                    return self.data[i];

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn at(self: *Sparse(T, order), r: u32, c: u32) T {
            // Unchecked version of get. Assumes r and c are valid.
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c)
                    return self.data[i];

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn set(self: *Sparse(T, order), allocator: std.mem.Allocator, r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c) {
                    self.data[i] = value;
                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._rlen or self.nnz == self._clen) {
                // Need more space
                var new_nnz = if (self.nnz * 2 > self.rows * self.cols) self.rows * self.cols else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i] = value;
            self.row[i] = r;
            self.col[i] = c;
            self.nnz += 1;
        }

        pub fn put(self: *Sparse(T, order), r: u32, c: u32, value: T) void {
            // Unchecked version of set. Assumes r and c are valid and there is space.
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c) {
                    self.data[i] = value;
                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i] = value;
            self.row[i] = r;
            self.col[i] = c;
            self.nnz += 1;
        }

        pub fn accumulate(self: *Sparse(T, order), allocator: std.mem.Allocator, r: u32, c: u32, value: anytype, ctx: anytype) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.row[i] == r and self.col[i] == c) {
                    try ops.add_(
                        &self.data[i],
                        self.data[i],
                        value,
                        ctx,
                    );

                    return;
                }

                if (comptime order == .col_major) {
                    if (self.col[i] > c or (self.col[i] == c and self.row[i] > r))
                        break;
                } else {
                    if (self.row[i] > r or (self.row[i] == r and self.col[i] > c))
                        break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._rlen or self.nnz == self._clen) {
                // Need more space
                var new_nnz = if (self.nnz * 2 > self.rows * self.cols) self.rows * self.cols else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to make space for new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.row[j] = self.row[j - 1];
                self.col[j] = self.col[j - 1];
            }

            self.data[i] = value;
            self.row[i] = r;
            self.col[i] = c;
            self.nnz += 1;
        }

        pub fn compile(self: *Sparse(T, order), allocator: std.mem.Allocator) !matrix.general.Sparse(T, order) {
            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.general.Sparse(T, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        pub fn compileCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, ctx: anytype) !matrix.general.Sparse(T, order) {
            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            var data: []T = try allocator.alloc(T, self.nnz);
            errdefer allocator.free(data);
            var idx: []u32 = try allocator.alloc(u32, self.nnz);
            errdefer allocator.free(idx);

            i = 0;
            while (i < self.nnz) : (i += 1) {
                data[i] = try ops.copy(self.data[i], ctx);
                idx[i] = if (comptime order == .col_major) self.row[i] else self.col[i];
            }

            const result = matrix.general.Sparse(T, order){
                .data = data.ptr,
                .idx = idx.ptr,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        fn removeTriangle(self: *Sparse(T, order), comptime uplo: Uplo, comptime diagonal: bool) void {
            var i: u32 = 0;
            var j: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                const r = self.row[i];
                const c = self.col[i];

                const keep = switch (comptime uplo) {
                    .upper => if (diagonal) r > c else r >= c,
                    .lower => if (diagonal) r < c else r <= c,
                };

                if (keep) {
                    if (i != j) {
                        self.data[j] = self.data[i];
                        self.row[j] = self.row[i];
                        self.col[j] = self.col[i];
                    }

                    j += 1;
                }
            }

            self.nnz = j;
        }

        pub fn compileSymmetric(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo) !matrix.symmetric.Sparse(T, uplo, order) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            self.removeTriangle(comptime uplo.invert(), false);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.symmetric.Sparse(T, uplo, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .size = self.rows,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        pub fn compileSymmetricCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, ctx: anytype) !matrix.symmetric.Sparse(T, uplo, order) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            var cpy = try self.copy(allocator, ctx);
            errdefer cpy.deinit(allocator);

            cpy.removeTriangle(comptime uplo.invert(), false);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) cpy.cols + 1 else cpy.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < cpy.nnz and cpy.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < cpy.nnz and cpy.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(cpy.col[0..cpy._clen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._rlen > cpy.nnz)
                    cpy.row = (try allocator.realloc(cpy.row[0..cpy._rlen], cpy.nnz)).ptr;
            } else {
                allocator.free(cpy.row[0..cpy._rlen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._clen > cpy.nnz)
                    cpy.col = (try allocator.realloc(cpy.col[0..cpy._clen], cpy.nnz)).ptr;
            }

            const result = matrix.symmetric.Sparse(T, uplo, order){
                .data = cpy.data,
                .idx = if (comptime order == .col_major) cpy.row else cpy.col,
                .ptr = ptr.ptr,
                .nnz = cpy.nnz,
                .size = cpy.rows,
                .flags = .{ .owns_data = true },
            };

            cpy.* = undefined;

            return result;
        }

        pub fn compileHermitian(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo) !matrix.hermitian.Sparse(T, uplo, order) {
            if (comptime !types.isComplex(T))
                @compileError("T must be a complex type");

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            self.removeTriangle(comptime uplo.invert(), false);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.hermitian.Sparse(T, uplo, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .size = self.rows,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        pub fn compileHermitianCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, ctx: anytype) !matrix.hermitian.Sparse(T, uplo, order) {
            if (comptime !types.isComplex(T))
                @compileError("T must be a complex type");

            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            var cpy = try self.copy(allocator, ctx);
            errdefer cpy.deinit(allocator);

            cpy.removeTriangle(comptime uplo.invert(), false);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) cpy.cols + 1 else cpy.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < cpy.nnz and cpy.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < cpy.nnz and cpy.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(cpy.col[0..cpy._clen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._rlen > cpy.nnz)
                    cpy.row = (try allocator.realloc(cpy.row[0..cpy._rlen], cpy.nnz)).ptr;
            } else {
                allocator.free(cpy.row[0..cpy._rlen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._clen > cpy.nnz)
                    cpy.col = (try allocator.realloc(cpy.col[0..cpy._clen], cpy.nnz)).ptr;
            }

            const result = matrix.hermitian.Sparse(T, uplo, order){
                .data = cpy.data,
                .idx = if (comptime order == .col_major) cpy.row else cpy.col,
                .ptr = ptr.ptr,
                .nnz = cpy.nnz,
                .size = cpy.rows,
                .flags = .{ .owns_data = true },
            };

            cpy.* = undefined;

            return result;
        }

        pub fn compileTriangular(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime diag: Diag) !matrix.triangular.Sparse(T, uplo, diag, order) {
            self.removeTriangle(comptime uplo.invert(), diag == .unit);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < self.nnz and self.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < self.nnz and self.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(self.col[0..self._clen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._rlen > self.nnz)
                    self.row = (try allocator.realloc(self.row[0..self._rlen], self.nnz)).ptr;
            } else {
                allocator.free(self.row[0..self._rlen]);

                if (self._dlen > self.nnz)
                    self.data = (try allocator.realloc(self.data[0..self._dlen], self.nnz)).ptr;

                if (self._clen > self.nnz)
                    self.col = (try allocator.realloc(self.col[0..self._clen], self.nnz)).ptr;
            }

            const result = matrix.triangular.Sparse(T, uplo, diag, order){
                .data = self.data,
                .idx = if (comptime order == .col_major) self.row else self.col,
                .ptr = ptr.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                .flags = .{ .owns_data = true },
            };

            self.* = undefined;

            return result;
        }

        pub fn compileTriangularCopy(self: *Sparse(T, order), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime diag: Diag, ctx: anytype) !matrix.triangular.Sparse(T, uplo, diag, order) {
            var cpy = try self.copy(allocator, ctx);
            errdefer cpy.deinit(allocator);

            cpy.removeTriangle(comptime uplo.invert(), diag == .unit);

            var ptr: []u32 = try allocator.alloc(u32, if (comptime order == .col_major) cpy.cols + 1 else cpy.rows + 1);
            errdefer allocator.free(ptr);
            ptr[0] = 0;

            var p: u32 = 0;
            var i: u32 = 0;
            while (p < ptr.len - 1) : (p += 1) {
                if (comptime order == .col_major) {
                    while (i < cpy.nnz and cpy.col[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                } else {
                    while (i < cpy.nnz and cpy.row[i] == p) : (i += 1) {}
                    ptr[p + 1] = i;
                }
            }

            if (comptime order == .col_major) {
                allocator.free(cpy.col[0..cpy._clen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._rlen > cpy.nnz)
                    cpy.row = (try allocator.realloc(cpy.row[0..cpy._rlen], cpy.nnz)).ptr;
            } else {
                allocator.free(cpy.row[0..cpy._rlen]);

                if (cpy._dlen > cpy.nnz)
                    cpy.data = (try allocator.realloc(cpy.data[0..cpy._dlen], cpy.nnz)).ptr;

                if (cpy._clen > cpy.nnz)
                    cpy.col = (try allocator.realloc(cpy.col[0..cpy._clen], cpy.nnz)).ptr;
            }

            const result = matrix.triangular.Sparse(T, uplo, diag, order){
                .data = cpy.data,
                .idx = if (comptime order == .col_major) cpy.row else cpy.col,
                .ptr = ptr.ptr,
                .nnz = cpy.nnz,
                .rows = cpy.rows,
                .cols = cpy.cols,
                .flags = .{ .owns_data = true },
            };

            cpy.* = undefined;

            return result;
        }

        fn copy(self: *const Sparse(T, order), allocator: std.mem.Allocator, ctx: anytype) !Sparse(T, order) {
            var data: []T = try allocator.alloc(T, self.nnz);
            errdefer allocator.free(data);
            var row: []u32 = try allocator.alloc(u32, self.nnz);
            errdefer allocator.free(row);
            var col: []u32 = try allocator.alloc(u32, self.nnz);
            errdefer allocator.free(col);

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                data[i] = try ops.copy(self.data[i], ctx);
                row[i] = self.row[i];
                col[i] = self.col[i];
            }

            return .{
                .data = data.ptr,
                .row = row.ptr,
                .col = col.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                ._dlen = self.nnz,
                ._rlen = self.nnz,
                ._clen = self.nnz,
                .flags = .{ .owns_data = true },
            };
        }
    };
}
