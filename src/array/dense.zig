//! Dense arrays are multi-dimensional arrays that store their data in a single
//! 1-dimensional array. Any dense array guarantees that:
//! - The entire data array is accessed.
//! - The strides are ordered in ascending order if the array is in column-major
//! order, or in descending order if it is in row-major order.

const std = @import("std");

const types = @import("../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureArray = types.EnsureArray;
const Numeric = types.Numeric;
const orderOf = types.orderOf;
const Coerce = types.Coerce;

const ops = @import("../ops.zig");
const int = @import("../int.zig");

const array = @import("../array.zig");
const Order = types.Order;
const Flags = array.Flags;
const Range = array.Range;

const strided = @import("strided.zig");
const Strided = strided.Strided;

const matrix = @import("../matrix.zig");

pub fn Dense(T: type, order: Order) type {
    if (!types.isNumeric(T))
        @compileError("Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        ndim: u32,
        shape: [array.max_dimensions]u32,
        strides: [array.max_dimensions]u32,
        size: u64,
        flags: Flags = .{},

        pub const empty: Dense(T, order) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** array.max_dimensions,
            .strides = .{0} ** array.max_dimensions,
            .size = 0,
            .flags = .{},
        };

        pub fn init(
            allocator: std.mem.Allocator,
            shape: []const u32,
        ) !Dense(T, order) {
            if (shape.len > array.max_dimensions) {
                return array.Error.TooManyDimensions;
            }

            for (shape) |dim| {
                if (dim == 0) {
                    return array.Error.ZeroDimension;
                }
            }

            var size: u64 = 1;
            var array_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var array_strides: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            if (shape.len > 0) {
                var i: u32 = 0;
                while (i < shape.len) : (i += 1) {
                    const idx: u32 = types.scast(u32, if (comptime order == .row_major) shape.len - i - 1 else i);

                    if (shape[idx] == 1) {
                        array_strides[idx] = 0; // No stride for the new dimension.
                    } else {
                        array_strides[idx] = types.scast(u32, size);
                    }

                    size *= shape[idx];

                    array_shape[i] = shape[i];
                }
            }

            return .{
                .data = (try allocator.alloc(T, size)).ptr,
                .ndim = types.scast(u32, shape.len),
                .shape = array_shape,
                .strides = array_strides,
                .size = size,
                .flags = .{
                    .owns_data = true,
                },
            };
        }

        pub fn full(
            allocator: std.mem.Allocator,
            shape: []const u32,
            value: anytype,
            ctx: anytype,
        ) !Dense(T, order) {
            comptime if (types.isArbitraryPrecision(T)) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            };

            var arr: Dense(T, order) = try .init(T, allocator, shape);
            errdefer arr.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                const casted_value: T = types.scast(T, value);

                var i: u32 = 0;
                while (i < arr.size) : (i += 1) {
                    arr.data[i] = casted_value;
                }
            } else {
                // Orientative for arbitrary precision types.
                arr.data[0] = try types.cast(T, value, ctx);

                var i: u32 = 1;
                errdefer _cleanup(T, arr.data[0..i], ctx);

                while (i < arr.size) : (i += 1) {
                    arr.data[i] = ops.copy(arr.data[0], ctx);
                }
            }

            return arr;
        }

        pub fn arange(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            step: anytype,
            ctx: anytype,
        ) !Dense(T, order) {
            comptime if (types.isComplex(T))
                @compileError("array.arange does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.arange: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.arange: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.arange: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            comptime if (!types.isNumeric(@TypeOf(step)))
                @compileError("array.arange: step must be numeric, got " ++ @typeName(@TypeOf(step)));

            comptime if (types.isComplex(@TypeOf(step)))
                @compileError("array.arange: step cannot be complex, got " ++ @typeName(@TypeOf(step)));

            comptime if (types.isArbitraryPrecision(T)) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            };

            var arr: Dense(T, order) = undefined;
            if (comptime !types.isArbitraryPrecision(T)) {
                const positive_step: bool = try ops.gt(step, 0, .{});
                if (ops.eq(step, 0, .{}) catch unreachable or
                    (ops.lt(stop, start, .{}) catch unreachable and positive_step) or
                    (ops.gt(stop, start, .{}) catch unreachable and !positive_step))
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

                arr = try init(T, allocator, &.{len});

                arr.data[0] = start_casted;
                var i: u32 = 0;
                while (i < len - 1) : (i += 1) {
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

                arr = try .init(T, allocator, &.{len});
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

                var i: u32 = 4;
                errdefer cleanup(T, arr.data[4..i], ctx);

                while (i < len - 1) : (i += 1) {
                    arr.data[i] = try ops.add(arr.data[i - 1], step_casted, ctx);
                }

                try ops.add_(&step_casted, step_casted, arr.data[len - 2], ctx);
                arr.data[len - 1] = step_casted;
            }

            return arr;
        }

        pub fn linspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            opts: struct {
                num: u32 = 50,
                endpoint: bool = true,
                retstep: ?*T = null,
            },
            ctx: anytype,
        ) !Dense(T, order) {
            comptime if (types.isComplex(T))
                @compileError("array.linspace does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.linspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.linspace: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.linspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (opts.num == 0)
                return array.Error.ZeroDimension;

            comptime if (types.isArbitraryPrecision(T)) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            };

            var arr: Dense(T, order) = try init(allocator, &.{opts.num});
            errdefer arr.deinit(allocator);

            if (comptime !types.isArbitraryPrecision(T)) {
                if (opts.num == 1) {
                    ops.set( // arr.data[0] = start
                        &arr.data[0],
                        start,
                        ctx,
                    ) catch unreachable;

                    if (opts.retstep) |*r| {
                        ops.set(
                            r,
                            0,
                            ctx,
                        ) catch unreachable;
                    }

                    return arr;
                } else if (opts.num == 2) {
                    if (opts.endpoint) {
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
                    }

                    if (opts.retstep) |*r| {
                        ops.set( // r = (arr.data[1] - arr.data[0]) / 2
                            r,
                            ops.div(
                                ops.sub(
                                    arr.data[1],
                                    arr.data[0],
                                    ctx,
                                ) catch unreachable,
                                2,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;
                    }

                    return arr;
                }

                const start_casted: T = types.scast(T, ctx);
                const stop_casted: T = types.scast(T, ctx);
                var step: T = ops.sub( // step = stop_casted - start_casted
                    stop_casted,
                    start_casted,
                    ctx,
                ) catch unreachable;

                if (opts.endpoint) {
                    ops.div_(
                        &step,
                        step,
                        opts.num - 1,
                        ctx,
                    ) catch unreachable;
                } else {
                    ops.div_(
                        &step,
                        step,
                        opts.num,
                        ctx,
                    ) catch unreachable;
                }

                if (opts.retstep) |*r| {
                    ops.set(
                        r,
                        step,
                        ctx,
                    ) catch unreachable;
                }

                if (opts.num == 3 and opts.endpoint) {
                    arr.data[0] = start_casted;
                    ops.add_( // arr.data[1] = start_casted + step
                        &arr.data[1],
                        start_casted,
                        step,
                        ctx,
                    ) catch unreachable;
                    arr.data[2] = stop_casted;

                    return arr;
                } else if (opts.num == 3 and !opts.endpoint) {
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
                var i: u32 = 1;
                while (i < opts.num - 2) : (i += 1) {
                    ops.add_( // arr.data[i] = arr.data[i - 1] + step
                        &arr.data[i],
                        arr.data[i - 1],
                        step,
                        ctx,
                    ) catch unreachable;
                }

                if (opts.endpoint) {
                    ops.add_( // arr.data[num - 2] = arr.data[num - 3] + step
                        &arr.data[opts.num - 2],
                        arr.data[opts.num - 3],
                        step,
                        ctx,
                    );
                    arr.data[opts.num - 1] = stop_casted;
                } else {
                    ops.add_( // arr.data[num - 2] = arr.data[num - 3] + step
                        &arr.data[opts.num - 2],
                        arr.data[opts.num - 3],
                        step,
                        ctx,
                    ) catch unreachable;
                    ops.sub_( // arr.data[num - 1] = stop_casted - step
                        &arr.data[opts.num - 1],
                        stop_casted,
                        step,
                        ctx,
                    ) catch unreachable;
                }
            } else {
                if (opts.num == 1) {
                    arr.data[0] = try ops.init(T, ctx);
                    errdefer ops.deinit(&arr.data[0], ctx);
                    try ops.set(&arr.data[0], start, ctx);

                    return arr;
                } else if (opts.num == 2 and opts.endpoint) {
                    arr.data[0] = try ops.init(T, ctx);
                    errdefer ops.deinit(&arr.data[0], ctx);
                    try ops.set(&arr.data[0], start, ctx);
                    arr.data[1] = try ops.init(T, ctx);
                    errdefer ops.deinit(&arr.data[1], ctx);
                    try ops.set(&arr.data[1], stop, ctx);

                    return arr;
                } else if (opts.num == 2 and !opts.endpoint) {
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
                if (opts.endpoint) {
                    try ops.div_(&step, step, opts.num - 1, ctx);
                } else {
                    try ops.div_(&step, step, opts.num, ctx);
                }

                if (opts.retstep) |*r| {
                    try ops.set(r, step, ctx);
                }

                if (opts.num == 3 and opts.endpoint) {
                    arr.data[0] = start_casted;
                    try ops.add_(&step, step, arr.data[0], ctx);
                    arr.data[1] = step;
                    arr.data[2] = stop_casted;

                    return arr;
                } else if (opts.num == 3 and !opts.endpoint) {
                    arr.data[0] = start_casted;
                    try ops.sub_(&stop_casted, stop_casted, step, ctx);
                    arr.data[2] = stop_casted;
                    try ops.add_(&step, step, arr.data[0], ctx);
                    arr.data[1] = step;

                    return arr;
                }

                arr.data[0] = start_casted;

                var i: u32 = 1;
                errdefer cleanup(T, allocator, arr.data[1..i]);
                while (i < opts.num - 2) : (i += 1) {
                    arr.data[i] = try ops.add(arr.data[i - 1], step, ctx);
                }

                if (opts.endpoint) {
                    try ops.add_(&step, step, arr.data[opts.num - 3], ctx);
                    arr.data[opts.num - 2] = step;
                    arr.data[opts.num - 1] = stop_casted;
                } else {
                    try ops.sub_(&stop_casted, stop_casted, step, ctx);
                    try ops.add_(&step, step, arr.data[opts.num - 3], ctx);
                    arr.data[opts.num - 2] = step;
                    arr.data[opts.num - 1] = stop_casted;
                }
            }

            return arr;
        }

        pub fn logspace(
            allocator: std.mem.Allocator,
            start: anytype,
            stop: anytype,
            base: anytype,
            opts: struct {
                num: u32 = 50,
                endpoint: bool = true,
            },
            ctx: anytype,
        ) !Dense(T, order) {
            comptime if (types.isComplex(T))
                @compileError("array.logspace does not support " ++ @typeName(T));

            comptime if (!types.isNumeric(@TypeOf(start)))
                @compileError("array.arange: start must be numeric, got " ++ @typeName(@TypeOf(start)));

            comptime if (types.isComplex(@TypeOf(start)))
                @compileError("array.logspace: start cannot be complex, got " ++ @typeName(@TypeOf(start)));

            comptime if (!types.isNumeric(@TypeOf(stop)))
                @compileError("array.logspace: stop must be numeric, got " ++ @typeName(@TypeOf(stop)));

            comptime if (types.isComplex(@TypeOf(stop)))
                @compileError("array.logspace: stop cannot be complex, got " ++ @typeName(@TypeOf(stop)));

            if (opts.num == 0)
                return array.Error.ZeroDimension;

            comptime if (types.isArbitraryPrecision(T)) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            };

            var arr: Dense(T, order) = try linspace(
                T,
                allocator,
                start,
                stop,
                .{
                    .num = opts.num,
                    .endpoint = opts.endpoint,
                    .retstep = null,
                },
                ctx,
            );
            errdefer arr.deinit(allocator);

            try ops.pow_(&arr, base, arr, ctx);

            return arr;
        }

        /// Cleans up the array by deinitializing its elements. If the array holds
        /// a fixed precision type this is does not do anything.
        pub fn cleanup(
            self: *Dense(T, order),
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

            if (comptime types.isArbitraryPrecision(T)) {
                _cleanup(self.data, ctx);
            } // For fixed precision types, this is a no-op.
        }

        /// Deinitializes the array, freeing its data if it owns it.
        ///
        /// If the array holds an arbitrary precision type, it will not free the
        /// elements; call `cleanup` first to do that.
        pub fn deinit(self: *Dense(T, order), allocator: ?std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.?.free(self.data[0..self.size]);
            }

            self.* = undefined;
        }

        pub fn get(self: *const Dense(T, order), position: []const u32) !T {
            try self._checkPosition(position);

            return self.data[self._index(position)];
        }

        pub inline fn at(self: *const Dense(T, order), position: []const u32) T {
            // Unchecked version of get. Assumes position is valid.
            return self.data[self._index(position)];
        }

        pub fn set(self: *const Dense(T, order), position: []const u32, value: T) !void {
            try self._checkPosition(position);

            self.data[self._index(position)] = value;
        }

        pub inline fn put(self: *const Dense(T, order), position: []const u32, value: T) void {
            // Unchecked version of set. Assumes position is valid.
            self.data[self._index(position)] = value;
        }

        pub fn asGeneralDenseMatrix(self: *const Dense(T, order), axes: ?[2]u32, position: []const u32) !matrix.dense.General(T, order) {
            if (self.ndim < 2)
                return array.Error.NotConvertible;

            const axes_: [2]u32 = axes orelse if (comptime order == .col_major) .{ 0, 1 } else .{ self.ndim - 2, self.ndim - 1 };

            if (axes_[0] >= self.ndim or axes_[1] >= self.ndim or axes_[0] == axes_[1])
                return array.Error.InvalidAxes;

            if (self.strides[axes_[0]] != 1 and self.strides[axes_[1]] != 1)
                return array.Error.NotConvertible; // At least one of the axes must have a stride of 1.

            if (position.len > self.ndim - 2)
                return array.Error.TooManyDimensions;

            const position_: [8]u32 = blk: {
                var pos: [8]u32 = .{0} ** 8;
                var passed: u32 = 0;
                var i: u32 = 0;
                while (i < self.ndim) : (i += 1) {
                    if (i == axes_[0] or i == axes_[1]) {
                        pos[i] = 0;
                        passed += 1;
                    } else {
                        pos[i] = position[i - passed];
                    }
                }

                break :blk pos;
            };

            try self._checkPosition(position_[0..self.ndim]);

            return .{
                .data = self.data + self._index(position_[0..self.ndim]),
                .rows = self.shape[axes_[0]],
                .cols = self.shape[axes_[1]],
                .strides = .{ self.strides[axes_[0]], self.strides[axes_[1]] },
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn reshape(self: *const Dense(T, order), shape: []const u32) !Dense(T, order) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len == 0)
                return array.Error.ZeroDimension;

            var new_size: u64 = 1;
            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var i: u32 = 0;
            while (i < shape.len) : (i += 1) {
                const idx: u32 = if (comptime order == .row_major) shape.len - i - 1 else i;

                new_strides[idx] = types.scast(u32, new_size);
                new_size *= shape[idx];

                new_shape[i] = shape[i];
            }

            if (new_size != self.size) {
                return array.Error.DimensionMismatch;
            }

            return .{
                .data = self.data,
                .ndim = shape.len,
                .shape = new_shape,
                .strides = new_strides,
                .size = new_size,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn ravel(self: *const Dense(T, order)) !Dense(T, order) {
            if (self.ndim == 0) {
                return array.Error.ZeroDimension;
            }

            return .{
                .data = self.data,
                .ndim = 1,
                .shape = .{self.size} ++ .{0} ** (array.max_dimensions - 1),
                .strides = .{1} ++ .{0} ** (array.max_dimensions - 1),
                .size = self.size,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn transpose(self: *const Dense(T, order), axes: ?[]const u32) !Strided(T, order) {
            const axes_: []const u32 =
                axes orelse
                array.trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axes_.len == 0) {
                return array.Error.ZeroDimension;
            }

            if (axes_.len != self.ndim) {
                return array.Error.DimensionMismatch;
            }

            if (!array.isPermutation(self.ndim, axes_)) {
                return array.Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]
            }

            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var size: u32 = 1;

            var i: u32 = 0;
            while (i < self.ndim) : (i += 1) {
                const idx: u32 = axes_[i];

                new_shape[i] = self.shape[idx];
                new_strides[i] = types.scast(i32, self.strides[idx]);
                size *= new_shape[i];
            }

            return .{
                .data = self.data,
                .ndim = self.ndim,
                .shape = new_shape,
                .strides = new_strides,
                .size = size,
                .offset = 0,
                .base_type = Dense(T, order),
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn broadcast(self: *const Dense(T, order), shape: []const u32) !Dense(T, order) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len < self.ndim)
                return array.Error.TooLittleDimensions;

            var new_shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var size: u32 = 1;

            var i: i32 = types.scast(i32, shape.len - 1);
            const diff: i32 = types.scast(i32, shape.len - self.ndim);
            while (i >= 0) : (i -= 1) {
                if (shape[types.scast(u32, i)] == 0)
                    return array.Error.ZeroDimension;

                if (i - diff >= 0) {
                    if (self.shape[types.scast(u32, i - diff)] != 1 and
                        self.shape[types.scast(u32, i - diff)] != shape[types.scast(u32, i)])
                        return array.Error.NotBroadcastable; // Broadcasting is not possible if the shapes do not match or are not compatible.

                    new_shape[types.scast(u32, i)] = int.max(self.shape[types.scast(u32, i - diff)], shape[types.scast(u32, i)]);
                    strides[types.scast(u32, i)] = self.strides[types.scast(u32, i - diff)];
                } else {
                    new_shape[types.scast(u32, i)] = shape[types.scast(u32, i)];
                    strides[types.scast(u32, i)] = 0; // No stride for the new dimensions.
                }

                size *= new_shape[types.scast(u32, i)];
            }

            return .{
                .data = self.data,
                .ndim = types.scast(u32, shape.len),
                .shape = new_shape,
                .strides = strides,
                .size = size,
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn slice(self: *const Dense(T, order), ranges: []const Range) !Strided(T, order) {
            if (ranges.len == 0 or ranges.len > self.ndim) {
                return error.DimensionMismatch;
            }

            var ndim: u32 = self.ndim;
            var size: u32 = 1;
            var shape: [array.max_dimensions]u32 = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]i32 = .{0} ** array.max_dimensions;
            var offset: u32 = 0;

            var i: u32 = 0;
            var j: u32 = 0;
            while (i < self.ndim) {
                const stride: i32 = types.scast(i32, self.strides[i]);

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

            return .{
                .data = self.data,
                .ndim = ndim,
                .shape = shape,
                .strides = strides,
                .size = size,
                .offset = offset,
                .base_type = Dense(T, order),
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        inline fn _cleanup(
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

        inline fn _index(
            self: *const Dense(T, order),
            position: []const u32,
        ) u32 {
            var idx: u32 = 0;
            var i: u32 = 0;
            while (i < position.len) : (i += 1) {
                idx += position[i] * self.strides[i];
            }

            return idx;
        }

        inline fn _checkPosition(self: *const Dense(T, order), position: []const u32) !void {
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

fn loop1(
    result: anytype,
    x: anytype,
    comptime op: anytype,
    comptime @"inline": bool,
    depth: u32,
    comptime order: types.IterationOrder,
    ir: u32,
    ix: u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 1) {
                    result.data[jr] = op(x.data[jx]);
                } else if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[jr] = try op(x.data[jx], ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 2) {
                    op(&result.data[jr], x.data[jx]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try op(&result.data[jr], x.data[jx], ctx);
                }
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
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

    if (std.mem.eql(u32, result.strides[0..result.ndim], x.strides[0..x.ndim])) {
        // Trivial loop
        //errdefer cleanup(ReturnType1(op, X), allocator, result.data[0..j]);

        const opinfo = @typeInfo(@TypeOf(op));
        var i: u32 = 0;
        while (i < result.size) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 1) {
                result.data[i] = op(x.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = try op(x.data[i], ctx);
            }
        }
    } else {
        try loop1(
            &result,
            &x,
            op,
            false,
            result.ndim - 1,
            comptime orderOf(@TypeOf(x)).toIterationOrder(),
            0,
            0,
            ctx,
        );
    }

    return result;
}

pub fn apply1_(
    o: anytype,
    x: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    const X: type = Numeric(@TypeOf(x));

    var xx: Dense(X) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
        if (std.mem.eql(u32, o.strides[0..o.ndim], x.strides[0..x.ndim])) {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 2) {
                    op_(&o.data[i], x.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    try op_(&o.data[i], x.data[i], ctx);
                }
            }

            return;
        } else {
            // Different order, but same shape
            xx = x;
        }
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
        0,
        0,
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
    ir: u32,
    iy: u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jy: u32 = iy;
        var j: u32 = 0;
        while (j < y.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[jr] = op(x, y.data[jy]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[jr] = try op(x, y.data[jy], ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op(&result.data[jr], x, y.data[jy]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op(&result.data[jr], x, y.data[jy], ctx);
                }
            }

            jr += result.strides[idx];
            jy += y.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jy: u32 = iy;
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
    ir: u32,
    ix: u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[jr] = op(x.data[jx], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[jr] = try op(x.data[jx], y, ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op(&result.data[jr], x.data[jx], y);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op(&result.data[jr], x.data[jx], y, ctx);
                }
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
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
    ir: u32,
    ix: u32,
    iy: u32,
    ctx: anytype,
) !void {
    if (depth == 0) {
        const opinfo = @typeInfo(@TypeOf(op));
        const idx: u32 = if (comptime order == .left_to_right) 0 else result.ndim - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
        var jy: u32 = iy;
        var j: u32 = 0;
        while (j < x.shape[idx]) : (j += 1) {
            if (comptime !@"inline") {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[jr] = op(x.data[jx], y.data[jy]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[jr] = try op(x.data[jx], y.data[jy], ctx);
                }
            } else {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op(&result.data[jr], x.data[jx], y.data[jy]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op(&result.data[jr], x.data[jx], y.data[jy], ctx);
                }
            }

            jr += result.strides[idx];
            jx += x.strides[idx];
            jy += y.strides[idx];
        }
    } else {
        const idx: u32 = if (comptime order == .left_to_right) depth else result.ndim - depth - 1;

        var jr: u32 = ir;
        var jx: u32 = ix;
        var jy: u32 = iy;
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

    if (comptime !types.isDenseArray(@TypeOf(x))) {
        var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, y.shape[0..y.ndim]);
        errdefer result.deinit(allocator);

        if (std.mem.eql(u32, result.strides[0..result.ndim], y.strides[0..y.ndim])) {
            // Trivial loop
            //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

            const opinfo = @typeInfo(@TypeOf(op));

            var i: u32 = 0;
            while (i < result.size) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[i], ctx);
                }
            }
        } else {
            try loop2_left(
                &result,
                x,
                &y,
                op,
                false,
                result.ndim - 1,
                comptime orderOf(@TypeOf(result)).toIterationOrder(),
                0,
                0,
                ctx,
            );
        }

        return result;
    } else if (comptime !types.isDenseArray(@TypeOf(y))) {
        var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, x.shape[0..x.ndim]);
        errdefer result.deinit(allocator);

        if (std.mem.eql(u32, result.strides[0..result.ndim], x.strides[0..x.ndim])) {
            //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

            const opinfo = @typeInfo(@TypeOf(op));

            var i: u32 = 0;
            while (i < result.size) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y, ctx);
                }
            }
        } else {
            try loop2_right(
                &result,
                &x,
                y,
                op,
                false,
                result.ndim - 1,
                comptime orderOf(@TypeOf(result)).toIterationOrder(),
                0,
                0,
                ctx,
            );
        }

        return result;
    }

    var xx: Dense(X) = undefined;
    var yy: Dense(Y) = undefined;
    if (std.mem.eql(u32, x.shape[0..x.ndim], y.shape[0..y.ndim])) {
        if (std.mem.eql(u32, x.strides[0..x.ndim], y.strides[0..y.ndim])) {
            // Trivial loop
            var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, x.shape[0..x.ndim]);
            errdefer result.deinit(allocator);

            //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

            const opinfo = @typeInfo(@TypeOf(op));

            var i: u32 = 0;
            while (i < result.size) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y.data[i], ctx);
                }
            }

            return result;
        } else {
            // Different order, but same shape
            xx = x;
            yy = y;
        }
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
        0,
        0,
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

    if (comptime !types.isDenseArray(@TypeOf(x))) {
        if (std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]) and
            std.mem.eql(u32, o.strides[0..o.ndim], y.strides[0..y.ndim]))
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op_(&o.data[i], x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op_(&o.data[i], x, y.data[i], ctx);
                }
            }

            return;
        }

        var yy: Dense(Y, orderOf(@TypeOf(y))) = undefined;
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
            0,
            0,
            ctx,
        );

        return;
    } else if (comptime !types.isDenseArray(@TypeOf(y))) {
        if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
            std.mem.eql(u32, o.strides[0..o.ndim], x.strides[0..x.ndim]))
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op_(&o.data[i], x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op_(&o.data[i], x.data[i], y, ctx);
                }
            }

            return;
        }

        // Different order, but same shape
        var xx: Dense(X, orderOf(@TypeOf(x))) = undefined;
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
            0,
            0,
            ctx,
        );

        return;
    }

    var xx: Dense(X) = undefined;
    var yy: Dense(Y) = undefined;
    if (std.mem.eql(u32, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
        std.mem.eql(u32, o.shape[0..o.ndim], y.shape[0..y.ndim]))
    {
        if (std.mem.eql(u32, o.strides[0..o.ndim], x.strides[0..x.ndim]) and
            std.mem.eql(u32, o.strides[0..o.ndim], y.strides[0..y.ndim]))
        {
            // Trivial loop
            const opinfo = @typeInfo(@TypeOf(op_));
            for (0..o.size) |i| {
                if (comptime opinfo.@"fn".params.len == 3) {
                    op_(&o.data[i], x.data[i], y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 4) {
                    try op_(&o.data[i], x.data[i], y.data[i], ctx);
                }
            }

            return;
        } else {
            // Different order, but same shape
            xx = x;
            yy = y;
        }
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
        0,
        0,
        0,
        ctx,
    );

    return;
}
