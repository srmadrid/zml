const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const int = @import("../../int.zig");

const array = @import("../../array.zig");
const Dense = @import("../../array.zig").Dense;
const Range = @import("../../array.zig").Range;

const dense = @import("../dense.zig");
const strided = @import("../strided.zig");
const Strided = strided.Strided;

pub inline fn init(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    order: array.Order,
) !Dense(T) {
    var size: u32 = 1;
    var array_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    var array_strides: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    if (shape.len > 0) {
        for (0..shape.len) |i| {
            const idx: usize = if (order == .row_major) shape.len - i - 1 else i;

            if (shape[idx] == 1) {
                array_strides[idx] = 0; // No stride for the new dimension.
            } else {
                array_strides[idx] = size;
            }

            size *= shape[idx];

            array_shape[i] = shape[i];
        }
    }

    return Dense(T){
        .data = try allocator.alloc(T, size),
        .ndim = types.scast(u32, shape.len),
        .shape = array_shape,
        .strides = array_strides,
        .size = size,
        .base = null,
        .flags = .{
            .order = order,
            .owns_data = true,
        },
        .kind = .{ .general = .{} },
    };
}

pub inline fn full(
    comptime T: type,
    allocator: std.mem.Allocator,
    shape: []const u32,
    value: anytype,
    order: array.Order,
    ctx: anytype,
) !Dense(T) {
    var arr: Dense(T) = try init(T, allocator, shape, order);
    errdefer arr.deinit(allocator);

    if (comptime !types.isArbitraryPrecision(T)) {
        const casted_value: T = types.scast(T, value);

        for (0..arr.size) |i| {
            arr.data[i] = casted_value;
        }
    } else {
        // Orientative for arbitrary precision types.
        arr.data[0] = try types.cast(T, value, ctx);

        var j: u32 = 1;
        errdefer cleanup(T, arr.data[0..j], ctx);

        for (1..arr.size) |i| {
            arr.data[i] = ops.copy(arr.data[0], ctx);

            j += 1;
        }
    }

    return arr;
}

pub inline fn arange(
    comptime T: type,
    allocator: std.mem.Allocator,
    start: anytype,
    stop: anytype,
    step: anytype,
    ctx: anytype,
) !Dense(T) {
    var arr: Dense(T) = undefined;
    if (comptime !types.isArbitraryPrecision(T)) {
        const positive_step: bool = try ops.gt(step, 0, .{});
        if (ops.eq(step, 0, .{}) catch unreachable or
            (ops.lt(stop, start, .{}) and positive_step) catch unreachable or
            (ops.gt(stop, start, .{}) and !positive_step) catch unreachable)
            return array.Error.InvalidRange;

        const start_casted: T = types.scast(T, start);
        const stop_casted: T = types.scast(T, ctx);
        const step_casted: T = types.scast(T, ctx);
        const diff: T = if (positive_step)
            ops.sub( // diff = stop_casted - start_casted
                stop_casted,
                start_casted,
                ctx,
            ) catch unreachable
        else
            ops.sub( // diff = start_casted - stop_casted
                start_casted,
                stop_casted,
                ctx,
            ) catch unreachable;

        var len_T: T = ops.div( // len_T = diff / step_casted
            diff,
            step_casted,
            ctx,
        ) catch unreachable;
        ops.abs_( // len_T = abs(len_T)
            &len_T,
            len_T,
            ctx,
        ) catch unreachable;
        ops.ceil_( // len_T = ceil(len_T)
            &len_T,
            len_T,
            ctx,
        ) catch unreachable;

        const len: u32 = types.scast(u32, len_T);

        if (len == 0)
            return array.Error.InvalidRange;

        arr = try init(T, allocator, &.{len}, .col_major);

        arr.data[0] = start_casted;
        for (0..len - 1) |i| {
            ops.add_( // arr.data[i] = arr.data[i - 1] + step_casted
                &arr.data[i],
                arr.data[i - 1],
                step_casted,
                ctx,
            ) catch unreachable;
        }
        arr.data[len - 1] = stop_casted;
    } else {
        // Orientative for arbitrary precision types.
        const positive_step: bool = try ops.gt(step, 0, .{});
        if (try ops.eq(step, 0, .{}) or
            (try ops.lt(stop, start, .{}) and positive_step) or
            (try ops.gt(stop, start, .{}) and !positive_step))
        {
            return array.Error.InvalidRange;
        }

        var start_casted: T = try ops.init(T, ctx);
        errdefer ops.deinit(&start_casted, ctx);
        try ops.set(&start_casted, start, ctx);

        var stop_casted: T = try ops.init(T, ctx);
        errdefer ops.deinit(&stop_casted, ctx);
        try ops.set(&stop_casted, stop, ctx);

        var step_casted: T = try ops.init(T, ctx);
        errdefer ops.deinit(&step_casted, ctx);
        try ops.set(&step_casted, step, ctx);

        var diff: T = if (positive_step)
            try ops.sub(stop_casted, start_casted, ctx)
        else
            try ops.sub(start_casted, stop_casted, ctx);
        errdefer ops.deinit(&diff, ctx);

        var len_T: T = try ops.div(diff, step_casted, ctx);
        errdefer ops.deinit(&len_T, ctx);
        try ops.abs_(&len_T, len_T, ctx);
        try ops.ceil_(&len_T, len_T, ctx);
        const len: u32 = types.scast(u32, len_T);

        if (len == 0) {
            return array.Error.InvalidRange;
        }

        arr = try init(T, allocator, &.{len}, .col_major);
        errdefer arr.deinit(allocator);
        switch (len) {
            1 => {
                arr.data[0] = start_casted;

                ops.deinit(&stop_casted, ctx);
                ops.deinit(&step_casted, ctx);
                ops.deinit(&diff, ctx);
                ops.deinit(&len_T, ctx);

                return arr;
            },
            2 => {
                arr.data[0] = start_casted;
                try ops.add_(&step_casted, step_casted, arr.data[0], ctx);
                arr.data[1] = step_casted;

                ops.deinit(&stop_casted, ctx);
                ops.deinit(&diff, ctx);
                ops.deinit(&len_T, ctx);

                return arr;
            },
            3 => {
                arr.data[0] = start_casted;
                try ops.add_(&stop_casted, arr.data[0], step_casted, ctx);
                arr.data[1] = stop_casted;
                try ops.add_(&step_casted, arr.data[1], step_casted, ctx);
                arr.data[2] = step_casted;

                ops.deinit(&diff, ctx);
                ops.deinit(&len_T, ctx);

                return arr;
            },
            4 => {
                arr.data[0] = start_casted;
                try ops.add_(&diff, arr.data[0], step_casted, ctx);
                arr.data[1] = diff;
                try ops.add_(&stop_casted, arr.data[1], step_casted, ctx);
                arr.data[2] = stop_casted;
                try ops.add_(&step_casted, arr.data[2], step_casted, ctx);
                arr.data[3] = step_casted;

                ops.deinit(&len_T, ctx);

                return arr;
            },
            else => {
                arr.data[0] = start_casted;
                try ops.add_(&len_T, arr.data[0], step_casted, ctx);
                arr.data[1] = len_T;
                try ops.add_(&diff, arr.data[1], step_casted, ctx);
                arr.data[2] = diff;
                try ops.add_(&stop_casted, arr.data[2], step_casted, ctx);
                arr.data[3] = stop_casted;
            },
        }

        var j: u32 = 4;
        errdefer cleanup(T, arr.data[4..j], ctx);

        for (4..len - 1) |i| {
            arr.data[i] = try ops.add(arr.data[i - 1], step_casted, ctx);
            j += 1;
        }

        try ops.add_(&step_casted, step_casted, arr.data[len - 2], ctx);
        arr.data[len - 1] = step_casted;
    }

    return arr;
}

pub inline fn linspace(
    comptime T: type,
    allocator: std.mem.Allocator,
    start: anytype,
    stop: anytype,
    num: u32,
    endpoint: bool,
    retstep: ?*T,
    ctx: anytype,
) !Dense(T) {
    var arr: Dense(T) = try init(T, allocator, &.{num}, .row_major);
    errdefer arr.deinit(allocator);

    if (comptime !types.isArbitraryPrecision(T)) {
        if (num == 1) {
            ops.set( // arr.data[0] = start
                &arr.data[0],
                start,
                ctx,
            ) catch unreachable;

            return arr;
        } else if (num == 2) {
            if (endpoint) {
                ops.set( // arr.data[0] = start
                    &arr.data[0],
                    start,
                    ctx,
                ) catch unreachable;

                ops.set( // arr.data[1] = stop
                    &arr.data[1],
                    stop,
                    ctx,
                ) catch unreachable;

                return arr;
            } else {
                ops.set( // arr.data[0] = start
                    &arr.data[0],
                    start,
                    ctx,
                ) catch unreachable;

                ops.div_( // arr.data[1] += (start + stop) / 2
                    &arr.data[1],
                    ops.add(
                        arr.data[0],
                        stop,
                        ctx,
                    ) catch unreachable,
                    2,
                    ctx,
                ) catch unreachable;

                return arr;
            }
        }

        const start_casted: T = types.scast(T, ctx);
        const stop_casted: T = types.scast(T, ctx);
        var step: T = ops.sub( // step = stop_casted - start_casted
            stop_casted,
            start_casted,
            ctx,
        ) catch unreachable;

        if (endpoint) {
            ops.div_(
                &step,
                step,
                num - 1,
                ctx,
            ) catch unreachable;
        } else {
            ops.div_(
                &step,
                step,
                num,
                ctx,
            ) catch unreachable;
        }

        if (retstep) |*s| {
            ops.set(
                s,
                step,
                ctx,
            ) catch unreachable;
        }

        if (num == 3 and endpoint) {
            arr.data[0] = start_casted;
            ops.add_( // arr.data[1] = start_casted + step
                &arr.data[1],
                start_casted,
                step,
                ctx,
            ) catch unreachable;
            arr.data[2] = stop_casted;

            return arr;
        } else if (num == 3 and !endpoint) {
            arr.data[0] = start_casted;
            ops.add_( // arr.data[1] = start_casted + step
                &arr.data[1],
                start_casted,
                step,
                ctx,
            ) catch unreachable;
            ops.sub_( // arr.data[2] = stop_casted - step
                &arr.data[2],
                stop_casted,
                step,
                ctx,
            ) catch unreachable;

            return arr;
        }

        arr.data[0] = start_casted;
        for (1..num - 2) |i| {
            ops.add_( // arr.data[i] = arr.data[i - 1] + step
                &arr.data[i],
                arr.data[i - 1],
                step,
                ctx,
            ) catch unreachable;
        }

        if (endpoint) {
            ops.add_( // arr.data[num - 2] = arr.data[num - 3] + step
                &arr.data[num - 2],
                arr.data[num - 3],
                step,
                ctx,
            );
            arr.data[num - 1] = stop_casted;
        } else {
            ops.add_( // arr.data[num - 2] = arr.data[num - 3] + step
                &arr.data[num - 2],
                arr.data[num - 3],
                step,
                ctx,
            ) catch unreachable;
            ops.sub_( // arr.data[num - 1] = stop_casted - step
                &arr.data[num - 1],
                stop_casted,
                step,
                ctx,
            ) catch unreachable;
        }
    } else {
        if (num == 1) {
            arr.data[0] = try ops.init(T, ctx);
            errdefer ops.deinit(&arr.data[0], ctx);
            try ops.set(&arr.data[0], start, ctx);

            return arr;
        } else if (num == 2 and endpoint) {
            arr.data[0] = try ops.init(T, ctx);
            errdefer ops.deinit(&arr.data[0], ctx);
            try ops.set(&arr.data[0], start, ctx);
            arr.data[1] = try ops.init(T, ctx);
            errdefer ops.deinit(&arr.data[1], ctx);
            try ops.set(&arr.data[1], stop, ctx);

            return arr;
        } else if (num == 2 and !endpoint) {
            arr.data[0] = try ops.init(T, ctx);
            errdefer ops.deinit(&arr.data[0], ctx);
            try ops.set(&arr.data[0], start, ctx);
            arr.data[1] = try ops.init(T, ctx);
            errdefer ops.deinit(&arr.data[1], ctx);
            try ops.set(&arr.data[1], stop, ctx);
            try ops.add_(&arr.data[1], arr.data[1], arr.data[0], ctx);
            try ops.div_(&arr.data[1], arr.data[1], 2, ctx);

            return arr;
        }

        var start_casted: T = try ops.init(T, ctx);
        errdefer ops.deinit(&start_casted, ctx);
        try ops.set(&start_casted, start, ctx);

        var stop_casted: T = try ops.init(T, ctx);
        errdefer ops.deinit(&stop_casted, ctx);
        try ops.set(&stop_casted, stop, ctx);

        var step: T = try ops.sub(stop_casted, start_casted, ctx);
        errdefer ops.deinit(&step, ctx);
        if (endpoint) {
            try ops.div_(&step, step, num - 1, ctx);
        } else {
            try ops.div_(&step, step, num, ctx);
        }

        if (retstep) |*s| {
            try ops.set(s, step, ctx);
        }

        if (num == 3 and endpoint) {
            arr.data[0] = start_casted;
            try ops.add_(&step, step, arr.data[0], ctx);
            arr.data[1] = step;
            arr.data[2] = stop_casted;

            return arr;
        } else if (num == 3 and !endpoint) {
            arr.data[0] = start_casted;
            try ops.sub_(&stop_casted, stop_casted, step, ctx);
            arr.data[2] = stop_casted;
            try ops.add_(&step, step, arr.data[0], ctx);
            arr.data[1] = step;

            return arr;
        }

        arr.data[0] = start_casted;

        var j: u32 = 1;
        errdefer cleanup(T, allocator, arr.data[1..j]);
        for (1..num - 2) |i| {
            arr.data[i] = try ops.add(arr.data[i - 1], step, ctx);

            j += 1;
        }

        if (endpoint) {
            try ops.add_(&step, step, arr.data[num - 3], ctx);
            arr.data[num - 2] = step;
            arr.data[num - 1] = stop_casted;
        } else {
            try ops.sub_(&stop_casted, stop_casted, step, ctx);
            try ops.add_(&step, step, arr.data[num - 3], ctx);
            arr.data[num - 2] = step;
            arr.data[num - 1] = stop_casted;
        }
    }

    return arr;
}

pub inline fn logspace(
    comptime T: type,
    allocator: std.mem.Allocator,
    start: anytype,
    stop: anytype,
    base: anytype,
    num: u32,
    endpoint: bool,
    ctx: anytype,
) !Dense(T) {
    var arr: Dense(T) = try linspace(T, allocator, start, stop, num, endpoint, null, ctx);
    errdefer arr.deinit(allocator);

    try ops.pow_(&arr, base, arr, ctx);

    return arr;
}

pub inline fn index(comptime T: type, arr: *const Dense(T), position: []const u32) u32 {
    var idx: u32 = 0;
    for (0..position.len) |i| {
        idx += position[i] * arr.strides[i];
    }

    return idx;
}

pub inline fn checkPosition(comptime T: type, arr: *const Dense(T), position: []const u32) !void {
    if (position.len > arr.ndim)
        return array.Error.DimensionMismatch;

    for (0..position.len) |i| {
        if (position[i] >= arr.shape[i]) {
            return array.Error.PositionOutOfBounds;
        }
    }
}

pub inline fn set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) !void {
    try checkPosition(T, arr, position);

    arr.data[index(T, arr, position)] = value;
}

pub inline fn _set(comptime T: type, arr: *const Dense(T), position: []const u32, value: T) void {
    // Assumes position is valid.
    arr.data[index(T, arr, position)] = value;
}

pub inline fn get(comptime T: type, arr: *const Dense(T), position: []const u32) !T {
    try checkPosition(T, arr, position);

    return arr.data[index(T, arr, position)];
}

pub inline fn _get(comptime T: type, arr: *const Dense(T), position: []const u32) T {
    // Assumes position is valid.
    return arr.data[index(T, arr, position)];
}

pub inline fn reshape(comptime T: type, arr: *const Dense(T), shape: []const u32) !Strided(T) {
    var new_size: u32 = 1;
    var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    var new_strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;

    for (0..shape.len) |i| {
        const idx: u32 = if (arr.flags.order == .row_major) shape.len - i - 1 else i;

        new_strides[idx] = types.scast(i32, new_size);
        new_size *= shape[idx];

        new_shape[i] = shape[i];
    }

    if (new_size != arr.size) {
        return error.DimensionMismatch;
    }

    return Strided(T){
        .data = arr.data,
        .ndim = shape.len,
        .shape = new_shape,
        .strides = new_strides,
        .size = new_size,
        .offset = 0, // No offset in reshaping
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order,
            .owns_data = false,
        },
    };
}

pub inline fn ravel(comptime T: type, arr: *const Dense(T)) Strided(T) {
    return Strided(T){
        .data = arr.data,
        .ndim = 1,
        .shape = .{arr.size} ++ .{0} ** (array.max_dimensions - 1),
        .strides = .{1} ++ .{0} ** (array.max_dimensions - 1),
        .size = arr.size,
        .offset = 0, // No offset in raveling.
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order,
            .owns_data = false,
        },
    };
}

pub inline fn transpose(
    comptime T: type,
    arr: *Dense(T),
    axes: []const u32,
) !Dense(T) {
    var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    var new_strides: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    var size: u32 = 1;

    for (0..arr.ndim) |i| {
        const idx: u32 = axes[i];

        new_shape[i] = arr.shape[idx];
        new_strides[i] = arr.strides[idx];
        size *= new_shape[i];
    }

    return Dense(T){
        .data = arr.data,
        .ndim = arr.ndim,
        .shape = new_shape,
        .strides = new_strides,
        .size = size,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Underlying order remains the same.
            .owns_data = false,
        },
        .kind = .{ .general = .{} },
    };
}

pub inline fn slice(comptime T: type, arr: *const Dense(T), ranges: []const Range) !Strided(T) {
    var ndim: u32 = arr.ndim;
    var size: u32 = 1;
    var shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
    var offset: u32 = 0;

    var i: u32 = 0;
    var j: u32 = 0;
    while (i < arr.ndim) {
        const stride: i32 = types.scast(i32, arr.metadata.dense.strides[i]);

        if (i >= ranges.len) {
            shape[j] = arr.shape[i];
            strides[j] = stride;
            size *= arr.shape[i];
            j += 1;
            i += 1;
            continue;
        }

        var range: Range = ranges[i];
        if (range.start != int.maxVal(u32) and range.start == range.stop) {
            return array.Error.InvalidRange;
        } else if (range.step > 0) {
            if (range.start != int.maxVal(u32) and range.start >= arr.shape[i] or
                (range.stop != int.maxVal(u32) and range.stop > arr.shape[i]))
                return array.Error.RangeOutOfBounds;
        } else if (range.step < 0) {
            if ((range.stop != int.maxVal(u32) and range.stop >= arr.shape[i]) or
                (range.start != int.maxVal(u32) and range.start > arr.shape[i]))
                return array.Error.RangeOutOfBounds;
        }

        var len_adjustment: u32 = 0;
        if (range.step > 0) {
            if (range.start == int.maxVal(u32)) {
                range.start = 0;
            }

            if (range.stop == int.maxVal(u32)) {
                range.stop = arr.shape[i];
            }
        } else if (range.step < 0) {
            if (range.start == int.maxVal(u32)) {
                range.start = arr.shape[i] - 1;
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

    return Strided(T){
        .data = arr.data,
        .ndim = ndim,
        .shape = shape,
        .strides = strides,
        .size = size,
        .offset = offset,
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .owns_data = false,
        },
    };
}

pub inline fn broadcast(
    comptime T: type,
    arr: *const Dense(T),
    shape: []const u32,
) !Strided(T) {
    var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
    var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
    var size: u32 = 1;

    var i: i32 = types.scast(i32, shape.len - 1);
    const diff: i32 = types.scast(i32, shape.len - arr.ndim);
    while (i >= 0) : (i -= 1) {
        if (shape[types.scast(u32, i)] == 0)
            return array.Error.ZeroDimension;

        if (i - diff >= 0) {
            if (arr.shape[types.scast(u32, i - diff)] != 1 and
                arr.shape[types.scast(u32, i - diff)] != shape[types.scast(u32, i)])
                return array.Error.NotBroadcastable; // Broadcasting is not possible if the shapes do not match or are not compatible.

            new_shape[types.scast(u32, i)] = try ops.max(arr.shape[types.scast(u32, i - diff)], shape[types.scast(u32, i)], .{});
            strides[types.scast(u32, i)] = types.scast(i32, arr.metadata.dense.strides[types.scast(u32, i - diff)]);
        } else {
            new_shape[types.scast(u32, i)] = shape[types.scast(u32, i)];
            strides[types.scast(u32, i)] = 0; // No stride for the new dimensions.
        }

        size *= new_shape[types.scast(u32, i)];
    }

    return Strided(T){
        .data = arr.data,
        .ndim = shape.len,
        .shape = new_shape,
        .strides = strides,
        .size = size,
        .offset = 0, // No offset in broadcasting.
        .base = if (arr.flags.owns_data) arr else arr.base,
        .flags = .{
            .order = arr.flags.order, // Although it is strided, knowing the underlying order is useful for efficient iteration.
            .owns_data = false,
        },
    };
}

pub inline fn cleanup(
    comptime T: type,
    data: []T,
    ctx: anytype,
) void {
    comptime if (types.isArbitraryPrecision(T)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    comptime if (!types.isArbitraryPrecision(T))
        return;

    for (data) |*value| {
        ops.deinit(value, ctx);
    }
}
