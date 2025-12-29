const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureArray = types.EnsureArray;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const orderOf = types.orderOf;

const ops = @import("../ops.zig");
const int = @import("../int.zig");

const array = @import("../array.zig");
const Order = types.Order;
const Flags = array.Flags;
const Range = array.Range;

const dense = @import("dense.zig");
const Dense = dense.Dense;

pub fn Strided(T: type, base_order: Order) type {
    if (!types.isNumeric(T))
        @compileError("Strided requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        ndim: u32,
        shape: [array.max_dimensions]u32,
        strides: [array.max_dimensions]i32,
        size: u64,
        offset: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_array = {};
        pub const is_strided = {};

        /// Numeric type
        pub const Numeric = T;

        pub const empty: Strided(T, base_order) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** array.max_dimensions,
            .strides = .{0} ** array.max_dimensions,
            .size = 0,
            .offset = 0,
            .flags = .{},
        };

        pub fn get(self: *const Strided(T, base_order), position: []const u32) !T {
            try self._checkPosition(position);

            return self.data[self._index(position)];
        }

        pub inline fn at(self: *const Strided(T, base_order), position: []const u32) T {
            // Unchecked version of get. Assumes position is valid.
            return self.data[self._index(position)];
        }

        pub fn set(self: *const Strided(T, base_order), position: []const u32, value: T) !void {
            try self._checkPosition(position);

            self.data[self._index(position)] = value;
        }

        pub inline fn put(self: *const Strided(T, base_order), position: []const u32, value: T) void {
            // Unchecked version of set. Assumes position is valid.
            self.data[self._index(position)] = value;
        }

        pub fn transpose(self: *const Strided(T, base_order), axes: ?[]const u32) !Strided(T, base_order) {
            const axes_: []const u32 =
                axes orelse
                array.trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axes_.len == 0)
                return array.Error.ZeroDimension;

            if (axes_.len != self.ndim)
                return array.Error.DimensionMismatch;

            if (!array.isPermutation(self.ndim, axes_))
                return array.Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]

            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var size: u32 = 1;

            var i: u32 = 0;
            while (i < self.ndim) : (i += 1) {
                const idx: u32 = axes_[i];

                new_shape[i] = self.shape[idx];
                new_strides[i] = self.strides[idx];
                size *= new_shape[i];
            }

            return Strided(T, base_order){
                .data = self.data,
                .ndim = self.ndim,
                .shape = new_shape,
                .strides = new_strides,
                .size = size,
                .offset = self.offset,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn broadcast(self: *const Strided(T, base_order), shape: []const u32) !Strided(T, base_order) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len < self.ndim)
                return array.Error.TooLittleDimensions;

            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var size: u32 = 1;

            var i: i32 = types.scast(i32, shape.len - 1);
            const diff: i32 = types.scast(i32, shape.len - self.ndim);
            while (i >= 0) : (i -= 1) {
                if (i - diff >= 0) {
                    new_shape[types.scast(u32, i)] = int.max(self.shape[types.scast(u32, i - diff)], shape[types.scast(u32, i)]);
                    strides[types.scast(u32, i)] = self.strides[types.scast(u32, i - diff)];
                } else {
                    new_shape[types.scast(u32, i)] = shape[types.scast(u32, i)];
                    strides[types.scast(u32, i)] = 0; // No stride for the new dimensions.
                }

                size *= new_shape[types.scast(u32, i)];
            }

            return Strided(T, base_order){
                .data = self.data,
                .ndim = types.scast(u32, shape.len),
                .shape = new_shape,
                .strides = strides,
                .size = size,
                .offset = self.offset,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn slice(self: *const Strided(T, base_order), ranges: []const Range) !Strided(T, base_order) {
            if (ranges.len == 0 or ranges.len > self.ndim)
                return error.DimensionMismatch;

            var ndim: u32 = self.ndim;
            var size: u32 = 1;
            var shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var offset: u32 = self.offset;

            var i: u32 = 0;
            var j: u32 = 0;
            while (i < self.ndim) {
                const stride: i32 = self.strides[i];

                if (i >= ranges.len) {
                    shape[j] = self.shape[i];
                    strides[j] = stride;
                    size *= self.shape[i];
                    j += 1;
                    i += 1;
                    continue;
                }

                var range: Range = ranges[i];
                if (range.start != int.maxVal(u32) and range.start == range.stop) {
                    return array.Error.InvalidRange;
                } else if (range.step > 0) {
                    if (range.start != int.maxVal(u32) and range.start >= self.shape[i] or
                        (range.stop != int.maxVal(u32) and range.stop > self.shape[i]))
                        return array.Error.RangeOutOfBounds;
                } else if (range.step < 0) {
                    if ((range.stop != int.maxVal(u32) and range.stop >= self.shape[i]) or
                        (range.start != int.maxVal(u32) and range.start > self.shape[i]))
                        return array.Error.RangeOutOfBounds;
                }

                var len_adjustment: u32 = 0;
                if (range.step > 0) {
                    if (range.start == int.maxVal(u32)) {
                        range.start = 0;
                    }

                    if (range.stop == int.maxVal(u32)) {
                        range.stop = self.shape[i];
                    }
                } else if (range.step < 0) {
                    if (range.start == int.maxVal(u32)) {
                        range.start = self.shape[i] - 1;
                    }

                    if (range.stop == int.maxVal(u32)) {
                        range.stop = 0;
                        len_adjustment = 1;
                    }
                }

                const len: u32 = range.len() + len_adjustment;
                if (len == 1) {
                    ndim -= 1;
                } else {
                    shape[j] = len;
                    strides[j] = stride * range.step;
                    size *= len;
                    j += 1;
                }

                if (stride < 0) {
                    offset -= range.start * types.scast(u32, int.abs(stride));
                } else {
                    offset += range.start * types.scast(u32, stride);
                }

                i += 1;
            }

            return Strided(T, base_order){
                .data = self.data,
                .ndim = ndim,
                .shape = shape,
                .strides = strides,
                .size = size,
                .offset = offset,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        inline fn _index(
            self: *const Strided(T, base_order),
            position: []const u32,
        ) u32 {
            var idx: u32 = self.offset;
            var i: u32 = 0;
            while (i < position.len) : (i += 1) {
                const stride: i32 = self.strides[i];
                if (stride < 0) {
                    idx -= position[i] * types.scast(u32, int.abs(stride));
                } else {
                    idx += position[i] * types.scast(u32, stride);
                }
            }

            return idx;
        }

        inline fn _checkPosition(self: *const Strided(T, base_order), position: []const u32) !void {
            if (position.len > self.ndim)
                return array.Error.DimensionMismatch;

            var i: u32 = 0;
            while (i < position.len) : (i += 1) {
                if (position[i] >= self.shape[i]) {
                    return array.Error.PositionOutOfBounds;
                }
            }
        }
    };
}

// pub inline fn cleanup(
//     comptime T: type,
//     arr: *const Strided(T, base_order),
//     num_elements: u32,
//     iteration_order: array.IterationOrder,
//     ctx: anytype,
// ) void {
//     comptime if (types.isArbitraryPrecision(T)) {
//         types.validateContext(
//             @TypeOf(ctx),
//             .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
//         );
//     } else {
//         types.validateContext(@TypeOf(ctx), .{});
//     };

//     comptime if (!types.isArbitraryPrecision(T))
//         return;

//     var iter: array.Iterator(T) = .init(arr);
//     var axis: u32 = 0;
//     if (iteration_order == .right_to_left) {
//         axis = arr.ndim - 1;
//     }
//     for (0..num_elements) |_| {
//         ops.deinit(arr.data[iter.index], ctx);

//         _ = iter.nextAO(axis, iteration_order);
//     }
// }

fn loop1(
    result: anytype,
    x: anytype,
    comptime op: anytype,
    comptime @"inline": bool,
    depth: u32,
    comptime order: types.IterationOrder,
    ir: if (types.isStridedArray(@TypeOf(result))) i32 else u32,
    ix: if (types.isStridedArray(@TypeOf(x))) i32 else u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jx: if (types.isStridedArray(@TypeOf(x))) i32 else u32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 1) {
                    result.data[types.scast(u32, jr)] = op(x.data[types.scast(u32, jx)]);
                } else if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[types.scast(u32, jr)] = try op(x.data[types.scast(u32, jx)], ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 2) {
                    op(&result.data[types.scast(u32, jr)], x.data[types.scast(u32, jx)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try op(&result.data[types.scast(u32, jr)], x.data[types.scast(u32, jx)], ctx);
                }
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jx: if (types.isStridedArray(@TypeOf(x))) i32 else u32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            try loop1(
                result,
                x,
                op,
                @"inline",
                depth - 1,
                order,
                jr,
                jx,
                ctx,
            );

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    }
}

pub fn apply1(
    allocator: std.mem.Allocator,
    x: anytype,
    comptime op: anytype,
    ctx: anytype,
) !Dense(ReturnType1(op, Numeric(@TypeOf(x))), orderOf(@TypeOf(x))) {
    const X: type = Numeric(@TypeOf(x));

    var result: Dense(ReturnType1(op, X), orderOf(@TypeOf(x))) = try .init(allocator, x.shape[0..x.ndim]);
    errdefer result.deinit(allocator);

    try loop1(
        &result,
        &x,
        op,
        false,
        result.ndim - 1,
        comptime orderOf(@TypeOf(x)).toIterationOrder(),
        0,
        types.scast(i32, x.offset),
        ctx,
    );

    return result;
}

pub fn apply1_(
    o: anytype,
    x: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    const X: type = Numeric(@TypeOf(x));

    var xx: Strided(X) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
        // Same shape
        xx = x;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });

        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim]))
            return array.Error.NotBroadcastable;

        xx = try x.broadcast(bct.shape[0..bct.ndim]);
    }

    try loop1(
        o,
        &xx,
        op_,
        true,
        o.ndim - 1,
        comptime orderOf(@TypeOf(o)).toIterationOrder(),
        if (comptime types.isStridedArray(@TypeOf(o))) types.scast(i32, o.offset) else 0,
        if (comptime types.isStridedArray(@TypeOf(xx))) types.scast(i32, xx.offset) else 0,
        ctx,
    );

    return;
}

fn loop2_left(
    result: anytype,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    comptime @"inline": bool,
    depth: u32,
    comptime order: types.IterationOrder,
    ir: if (types.isStridedArray(@TypeOf(result))) i32 else u32,
    iy: if (types.isStridedArray(@TypeOf(y))) i32 else u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jy: if (types.isStridedArray(@TypeOf(y))) i32 else u32 = iy;
        var j: u32 = 0;
        while (j < y.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[types.scast(u32, jr)] = op(x, y.data[types.scast(u32, jy)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[types.scast(u32, jr)] = try op(x, y.data[types.scast(u32, jy)], ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op(&result.data[types.scast(u32, jr)], x, y.data[types.scast(u32, jy)]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op(&result.data[types.scast(u32, jr)], x, y.data[types.scast(u32, jy)], ctx);
                }
            }

            jr += result.strides[idx];
            jy += y.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jy: if (types.isStridedArray(@TypeOf(y))) i32 else u32 = iy;
        var j: u32 = 0;
        while (j < y.shape[idx]) : (j += 1) {
            try loop2_left(
                result,
                x,
                y,
                op,
                @"inline",
                depth - 1,
                order,
                jr,
                jy,
                ctx,
            );

            jr += result.strides[idx];
            jy += y.strides[idx];
        }
    }
}

fn loop2_right(
    result: anytype,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    comptime @"inline": bool,
    depth: u32,
    comptime order: types.IterationOrder,
    ir: if (types.isStridedArray(@TypeOf(result))) i32 else u32,
    ix: if (types.isStridedArray(@TypeOf(x))) i32 else u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jx: if (types.isStridedArray(@TypeOf(x))) i32 else u32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[types.scast(u32, jr)] = op(x.data[types.scast(u32, jx)], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[types.scast(u32, jr)] = try op(x.data[types.scast(u32, jx)], y, ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op(&result.data[types.scast(u32, jr)], x.data[types.scast(u32, jx)], y);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op(&result.data[types.scast(u32, jr)], x.data[types.scast(u32, jx)], y, ctx);
                }
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jx: if (types.isStridedArray(@TypeOf(x))) i32 else u32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            try loop2_right(
                result,
                x,
                y,
                op,
                @"inline",
                depth - 1,
                order,
                jr,
                jx,
                ctx,
            );

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    }
}

fn loop2(
    result: anytype,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    comptime @"inline": bool,
    depth: u32,
    comptime order: types.IterationOrder,
    ir: if (types.isStridedArray(@TypeOf(result))) i32 else u32,
    ix: if (types.isStridedArray(@TypeOf(x))) i32 else u32,
    iy: if (types.isStridedArray(@TypeOf(y))) i32 else u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jx: if (types.isStridedArray(@TypeOf(x))) i32 else u32 = ix;
        var jy: if (types.isStridedArray(@TypeOf(y))) i32 else u32 = iy;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[types.scast(u32, jr)] = op(x.data[types.scast(u32, jx)], y.data[types.scast(u32, jy)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[types.scast(u32, jr)] = try op(x.data[types.scast(u32, jx)], y.data[types.scast(u32, jy)], ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op(&result.data[types.scast(u32, jr)], x.data[types.scast(u32, jx)], y.data[types.scast(u32, jy)]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op(&result.data[types.scast(u32, jr)], x.data[types.scast(u32, jx)], y.data[types.scast(u32, jy)], ctx);
                }
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
            jy += y.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: if (types.isStridedArray(@TypeOf(result))) i32 else u32 = ir;
        var jx: if (types.isStridedArray(@TypeOf(x))) i32 else u32 = ix;
        var jy: if (types.isStridedArray(@TypeOf(y))) i32 else u32 = iy;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            try loop2(
                result,
                x,
                y,
                op,
                @"inline",
                depth - 1,
                order,
                jr,
                jx,
                jy,
                ctx,
            );

            jr += result.strides[idx];
            jx += x.strides[idx];
            jy += y.strides[idx];
        }
    }
}

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));

    if (comptime !types.isStridedArray(@TypeOf(x))) {
        var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, y.shape[0..y.ndim]);
        errdefer result.deinit(allocator);

        try loop2_left(
            &result,
            x,
            &y,
            op,
            false,
            result.ndim - 1,
            comptime orderOf(@TypeOf(result)).toIterationOrder(),
            0,
            if (comptime types.isStridedArray(@TypeOf(y))) types.scast(i32, y.offset) else 0,
            ctx,
        );

        return result;
    } else if (comptime !types.isStridedArray(@TypeOf(y))) {
        var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, x.shape[0..x.ndim]);
        errdefer result.deinit(allocator);

        try loop2_right(
            &result,
            &x,
            y,
            op,
            false,
            result.ndim - 1,
            comptime orderOf(@TypeOf(result)).toIterationOrder(),
            0,
            if (comptime types.isStridedArray(@TypeOf(x))) types.scast(i32, x.offset) else 0,
            ctx,
        );

        return result;
    }

    var xx: Strided(X) = undefined;
    var yy: Strided(Y) = undefined;
    if (std.mem.eql(u32, x.shape[0..x.ndim], y.shape[0..y.ndim])) {
        // Same shape
        xx = x;
        yy = y;
    } else {
        const bct = try array.broadcastShapes(&.{ x.shape[0..x.ndim], y.shape[0..y.ndim] });
        xx = try x.broadcast(bct.shape[0..bct.ndim]);
        yy = try y.broadcast(bct.shape[0..bct.ndim]);
    }

    var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, xx.shape[0..xx.ndim]);
    errdefer result.deinit(allocator);

    try loop2(
        &result,
        &xx,
        &yy,
        op,
        false,
        result.ndim - 1,
        comptime Order.resolve3(orderOf(@TypeOf(x)), orderOf(@TypeOf(y)), orderOf(@TypeOf(result))).toIterationOrder(),
        0,
        if (comptime types.isStridedArray(@TypeOf(xx))) types.scast(i32, xx.offset) else 0,
        if (comptime types.isStridedArray(@TypeOf(yy))) types.scast(i32, yy.offset) else 0,
        ctx,
    );

    return result;
}

pub fn apply2_(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));

    if (comptime !types.isStridedArray(@TypeOf(x))) {
        var yy: Strided(Y) = undefined;
        if (std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
            yy = y;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
            if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            yy = try y.broadcast(bct.shape[0..bct.ndim]);
        }

        try loop2_left(
            o,
            x,
            &yy,
            op_,
            true,
            o.ndim - 1,
            comptime orderOf(@TypeOf(o)).toIterationOrder(),
            if (comptime types.isStridedArray(@TypeOf(o))) types.scast(i32, o.offset) else 0,
            if (comptime types.isStridedArray(@TypeOf(yy))) types.scast(i32, yy.offset) else 0,
            ctx,
        );

        return;
    } else if (comptime !types.isStridedArray(@TypeOf(y))) {
        var xx: Strided(X) = undefined;
        if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
            xx = x;
        } else {
            const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
            if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
                return array.Error.NotBroadcastable;
            }

            xx = try x.broadcast(bct.shape[0..bct.ndim]);
        }

        try loop2_right(
            o,
            &xx,
            y,
            op_,
            true,
            o.ndim - 1,
            comptime orderOf(@TypeOf(o)).toIterationOrder(),
            if (comptime types.isStridedArray(@TypeOf(o))) types.scast(i32, o.offset) else 0,
            if (comptime types.isStridedArray(@TypeOf(xx))) types.scast(i32, xx.offset) else 0,
            ctx,
        );

        return;
    }

    var xx: Strided(X) = undefined;
    var yy: Strided(Y) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]))
    {
        xx = x;
        yy = y;
    } else {
        const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim], y.shape[0..y.ndim] });

        if (!std.mem.eql(u32, bct.shape[0..bct.ndim], o.shape[0..o.ndim]))
            return array.Error.NotBroadcastable;

        xx = try x.broadcast(bct.shape[0..bct.ndim]);
        yy = try y.broadcast(bct.shape[0..bct.ndim]);
    }

    try loop2(
        o,
        &xx,
        &yy,
        op_,
        true,
        o.ndim - 1,
        comptime Order.resolve3(orderOf(@TypeOf(x)), orderOf(@TypeOf(y)), orderOf(@TypeOf(o))).toIterationOrder(),
        if (comptime types.isStridedArray(@TypeOf(o))) types.scast(i32, o.offset) else 0,
        if (comptime types.isStridedArray(@TypeOf(xx))) types.scast(i32, xx.offset) else 0,
        if (comptime types.isStridedArray(@TypeOf(yy))) types.scast(i32, yy.offset) else 0,
        ctx,
    );

    return;
}
