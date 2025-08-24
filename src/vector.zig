const std = @import("std");

const types = @import("types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;

pub const Flags = packed struct {
    owns_data: bool = true,
};

pub fn Vector(T: type) type {
    if (!types.isNumeric(T))
        @compileError("Vector requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        len: usize,
        flags: Flags = .{},

        pub const empty = Vector(T){
            .data = &.{},
            .len = 0,
            .flags = .{ .owns_data = false },
        };

        pub fn init(allocator: std.mem.Allocator, len: usize) !Vector(T) {
            if (len == 0)
                return Error.ZeroLength;

            return .{
                .data = (try allocator.alloc(T, len)).ptr,
                .len = len,
                .flags = .{ .owns_data = true },
            };
        }

        pub fn deinit(self: *Vector(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.len]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Vector(T), index: usize) !T {
            if (index >= self.len)
                return Error.PositionOutOfBounds;

            return self.data[index];
        }

        pub inline fn at(self: *const Vector(T), index: usize) T {
            // Unchecked version of get. Assumes index is valid.
            return self.data[index];
        }

        pub fn set(self: *Vector(T), index: usize, value: T) !void {
            if (index >= self.len)
                return Error.PositionOutOfBounds;

            self.data[index] = value;
        }

        pub inline fn put(self: *Vector(T), index: usize, value: T) void {
            // Unchecked version of set. Assumes index is valid.
            self.data[index] = value;
        }
    };
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !Vector(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = Vector(ReturnType2(op, X, Y));

    if (comptime !types.isGeneralMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.len);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < result.len) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x, y.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x, y.data[i], ctx);
            }
        }

        return result;
    } else if (comptime !types.isGeneralMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.len);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < result.len) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[i], y);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[i], y, ctx);
            }
        }

        return result;
    }

    if (x.len != y.len)
        return Error.DimensionMismatch;

    var result: R = try .init(allocator, x.len);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    var i: u32 = 0;
    while (i < result.len) : (i += 1) {
        if (comptime opinfo.@"fn".params.len == 2) {
            result.data[i] = op(x.data[i], y.data[i]);
        } else if (comptime opinfo.@"fn".params.len == 3) {
            result.data[i] = try op(x.data[i], y.data[i], ctx);
        }
    }

    return result;
}

pub const Error = error{
    ZeroLength,
    PositionOutOfBounds,
    DimensionMismatch,
};
