const std = @import("std");

const types = @import("../types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const int = @import("../int.zig");
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const linalg = @import("../linalg.zig");

const vector = @import("../vector.zig");
const Flags = vector.Flags;

pub fn Sparse(T: type) type {
    if (!types.isNumeric(T))
        @compileError("vector.Sparse requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        idx: [*]u32,
        nnz: u32,
        len: u32,
        _dlen: u32, // allocated length of data
        _ilen: u32, // allocated length of idx
        flags: Flags = .{},

        pub const empty = Sparse(T){
            .data = &.{},
            .idx = &.{},
            .nnz = 0,
            .len = 0,
            ._dlen = 0,
            ._ilen = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(allocator: std.mem.Allocator, len: u32, nnz: u32) !Sparse(T) {
            if (len == 0)
                return vector.Error.ZeroLength;

            if (nnz == 0 or nnz > len)
                return vector.Error.DimensionMismatch;

            const data: []T = try allocator.alloc(T, nnz);
            errdefer allocator.free(data);

            return .{
                .data = data.ptr,
                .idx = (try allocator.alloc(u32, nnz)).ptr,
                .nnz = 0,
                .len = len,
                ._dlen = nnz,
                ._ilen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn deinit(self: *Sparse(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.idx[0..self._ilen]);
            }

            self.* = undefined;
        }

        pub fn reserve(self: *Sparse(T), allocator: std.mem.Allocator, new_nnz: u32) !void {
            if (self.flags.owns_data == false)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._ilen)
                return;

            if (new_nnz > self.len)
                return vector.Error.DimensionMismatch;

            if (new_nnz > self._dlen) {
                self.data = try allocator.realloc(self.data[0..self._dlen], new_nnz).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._ilen) {
                self.idx = try allocator.realloc(self.idx[0..self._ilen], new_nnz).ptr;
                self._ilen = new_nnz;
            }
        }

        pub fn get(self: *const Sparse(T), index: u32) !T {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index)
                    return self.data[i]
                else if (self.idx[i] > index)
                    break;
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn at(self: *Sparse(T), index: u32) T {
            // Unchecked version of get. Assumes index is valid.
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index)
                    return self.data[i]
                else if (self.idx[i] > index)
                    break;
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        pub fn set(self: *Sparse(T), allocator: std.mem.Allocator, index: u32, value: T) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            if (self.flags.owns_data == false)
                return;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index) {
                    self.data[i] = value;
                    return;
                } else if (self.idx[i] > index) {
                    break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._ilen) {
                const new_nnz = if (self._dlen == 0) 1 else self._dlen * 2;
                try self.reserve(allocator, new_nnz);
            }

            // Shift elements to the right to make space for the new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.idx[j] = self.idx[j - 1];
            }

            self.data[i] = value;
            self.idx[i] = index;
            self.nnz += 1;
        }

        pub fn put(self: *Sparse(T), index: u32, value: T) void {
            // Unchecked version of set. Assumes index is valid and there is space.
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index) {
                    self.data[i] = value;
                    return;
                } else if (self.idx[i] > index) {
                    break;
                }
            }

            // Shift elements to the right to make space for the new element
            var j: u32 = self.nnz;
            while (j > i) : (j -= 1) {
                self.data[j] = self.data[j - 1];
                self.idx[j] = self.idx[j - 1];
            }

            self.data[i] = value;
            self.idx[i] = index;
            self.nnz += 1;
        }
    };
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !Sparse(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = Sparse(ReturnType2(op, Numeric(X), Numeric(Y)));

    if (comptime !types.isSparseVector(@TypeOf(x))) {
        var result: R = try .init(allocator, y.len, y.nnz);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < result.nnz) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x, y.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x, y.data[i], ctx);
            }

            result.idx[i] = y.idx[i];
            result.nnz += 1;
        }

        return result;
    } else if (comptime !types.isSparseVector(@TypeOf(y))) {
        var result: R = try .init(allocator, x.len);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < result.nnz) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[i], y);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[i], y, ctx);
            }

            result.idx[i] = x.idx[i];
            result.nnz += 1;
        }

        return result;
    }

    if (x.len != y.len)
        return vector.Error.DimensionMismatch;

    var result: R = try .init(allocator, x.len, int.min(x.nnz + y.nnz, x.len));
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var i: u32 = 0;
    var j: u32 = 0;
    while (i < x.nnz and j < y.nnz) {
        if (x.idx[i] == y.idx[j]) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[result.nnz] = op(x.data[i], y.data[j]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[result.nnz] = try op(x.data[i], y.data[j], ctx);
            }

            result.idx[result.nnz] = x.idx[i];
            result.nnz += 1;
            i += 1;
            j += 1;
        } else if (x.idx[i] < y.idx[j]) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[result.nnz] = op(x.data[i], constants.zero(Numeric(Y), .{}) catch unreachable);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[result.nnz] = try op(x.data[i], constants.zero(Numeric(Y), .{}) catch unreachable, ctx);
            }

            result.idx[result.nnz] = x.idx[i];
            result.nnz += 1;
            i += 1;
        } else {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[result.nnz] = op(constants.zero(Numeric(X), .{}) catch unreachable, y.data[j]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[result.nnz] = try op(constants.zero(Numeric(X), .{}) catch unreachable, y.data[j], ctx);
            }

            result.idx[result.nnz] = y.idx[j];
            result.nnz += 1;
            j += 1;
        }
    }

    while (i < x.nnz) : (i += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[result.nnz] = op(x.data[i], constants.zero(Numeric(Y), .{}) catch unreachable);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[result.nnz] = try op(x.data[i], constants.zero(Numeric(Y), .{}) catch unreachable, ctx);
        }

        result.idx[result.nnz] = x.idx[i];
        result.nnz += 1;
    }

    while (j < y.nnz) : (j += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[result.nnz] = op(constants.zero(Numeric(X), .{}) catch unreachable, y.data[j]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[result.nnz] = try op(constants.zero(Numeric(X), .{}) catch unreachable, y.data[j], ctx);
        }

        result.idx[result.nnz] = y.idx[j];
        result.nnz += 1;
    }

    return result;
}
