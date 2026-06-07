const std = @import("std");

const meta = @import("../meta.zig");
const Layout = meta.Layout;

const int = @import("../int.zig");

const numeric = @import("../numeric.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

/// Dense `n`-dimensional array type, represented as a contiguous array of
/// elements of type `N` and a set of strides.
pub fn Dense(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.array.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        ndim: usize,
        shape: [array.max_dimensions]usize,
        strides: [array.max_dimensions]isize,
        flags: array.Flags,

        // Type signatures
        pub const is_array = true;
        pub const is_dense = true;

        // Numeric type
        pub const Numeric = N;

        pub const empty: array.Dense(N) = .{
            .data = &.{},
            .ndim = 0,
            .shape = .{0} ** array.max_dimensions,
            .strides = .{0} ** array.max_dimensions,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `array.Dense(N)` with the specified shape. Strides
        /// are ordered in ascending order if `layout` is column major, or in
        /// descending order otherwise.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `layout` (`Layout`): The memory layout of the array.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The newly initialized array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.TooManyDimensions`: If `shape.len` is larger than
        ///   `array.max_dimensions`.
        /// * `array.Error.ZeroDimension`: If `shape.len` or any dimension is
        ///   zero.
        pub fn init(allocator: std.mem.Allocator, shape: []const usize, layout: Layout) !array.Dense(N) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len == 0)
                return array.Error.ZeroDimension;

            for (shape) |dim| {
                if (dim == 0)
                    return array.Error.ZeroDimension;
            }

            var size: usize = 1;
            var array_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var array_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
            var i: usize = 0;
            while (i < shape.len) : (i += 1) {
                const idx: usize = if (layout == .col_major) i else shape.len - i - 1;

                if (shape[idx] == 1)
                    array_strides[idx] = 0 // No stride for the new dimension.
                else
                    array_strides[idx] = numeric.cast(isize, size);

                size *= shape[idx];
                array_shape[i] = shape[i];
            }

            return .{
                .data = (try allocator.alloc(N, size)).ptr,
                .ndim = shape.len,
                .shape = array_shape,
                .strides = array_strides,
                .flags = .{ .owns_data = true },
            };
        }

        // initBuffer

        /// Initializes a new `array.Dense(N)` with the specified shape, with
        /// all elements set to the specified value. Strides are ordered in
        /// ascending order if `layout` is column major, or in descending order
        /// otherwise.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `value` (`N`): The value to fill the array with.
        /// * `layout` (`Layout`): The memory layout of the array.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The newly initialized array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.TooManyDimensions`: If `shape.len` is larger than
        ///   `array.max_dimensions`.
        /// * `array.Error.ZeroDimension`: If `shape.len` or any dimension is
        ///   zero.
        pub fn initValue(allocator: std.mem.Allocator, shape: []const usize, value: N, layout: Layout) !array.Dense(N) {
            const arr: array.Dense(N) = try .init(allocator, shape, layout);

            const size = arr._size();
            var i: usize = 0;
            while (i < size) : (i += 1) {
                arr.data[i] = value;
            }

            return arr;
        }

        /// Initializes a new `array.Dense(N)` with the specified shape, with
        /// all elements set by calling the specified function with the given
        /// arguments. Strides are ordered in ascending order if `layout` is
        /// column major, or in descending order otherwise.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `@"fn"` (`anytype`): The function to call to fill the array.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        /// * `layout` (`Layout`): The memory layout of the array.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The newly initialized array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.TooManyDimensions`: If `shape.len` is larger than
        ///   `array.max_dimensions`.
        /// * `array.Error.ZeroDimension`: If `shape.len` or any dimension is
        ///   zero.
        pub fn initFn(allocator: std.mem.Allocator, shape: []const usize, comptime @"fn": anytype, args: anytype, layout: Layout) !array.Dense(N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.array.Dense(N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var arr: array.Dense(N) = try .init(allocator, shape, layout);
            errdefer arr.deinit(allocator);

            const size = arr._size();
            var i: usize = 0;
            while (i < size) : (i += 1) {
                arr.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return arr;
        }

        /// Initializes a new `array.Dense(N)` with all elements set as evenly
        /// spaced values within a given interval.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `start` (`N`): Start of the interval (inclusive).
        /// * `stop` (`N`): End of the interval (exclusive).
        /// * `step` (`N`): Spacing between values.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The newly initialized array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.InvalidRange`: If `start`, `stop` and `step` form an
        ///   impossible interval.
        pub fn initArange(allocator: std.mem.Allocator, start: N, stop: N, step: N) !array.Dense(N) {
            comptime if (meta.isComplex(N))
                @compileError("zsl.array.Dense(N).initArange: not defined for complex N, got \n\tN = " ++ @typeName(N) ++ "\n");

            const positive_step: bool = numeric.gt(step, 0);
            if (numeric.eq(step, numeric.zero(N)) or
                (numeric.lt(stop, start) and positive_step) or
                (numeric.gt(stop, start) and !positive_step))
                return array.Error.InvalidRange;

            const diff: N = if (positive_step)
                numeric.sub(stop, start)
            else
                numeric.sub(start, stop);

            const len: usize = numeric.cast(usize, numeric.ceil(numeric.abs(numeric.div(diff, step))));
            if (len == 0)
                return array.Error.InvalidRange;

            const arr: array.Dense(N) = try .init(allocator, &.{len}, .col_major);

            arr.data[0] = start;
            var i: usize = 1;
            while (i < len) : (i += 1) {
                numeric.addInto(&arr.data[i], arr.data[i - 1], step);
            }

            return arr;
        }

        /// Initializes a new `array.Dense(N)` with all elements set as
        /// `opts.num` evenly spaced samples, calculated over a given interval.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `start` (`N`): Start of the interval (inclusive).
        /// * `stop` (`N`): End of the interval (inclusive if `opts.endpoint`,
        ///   exclusive otherwise).
        /// * `opts`: Optional parameters:
        ///   * `num` (`usize = 50`): Number of samples to generate.
        ///   * `endpoint` (`bool = true`): Wether the range is inclusive or
        ///     exclusive.
        ///   * `retstep` (`?*N = null`): Output pointer to optionally store the
        ///     step.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The newly initialized array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.ZeroDimension`: If `opts.num` is zero.
        pub fn initLinspace(
            allocator: std.mem.Allocator,
            start: N,
            stop: N,
            opts: struct {
                num: usize = 50,
                endpoint: bool = true,
                retstep: ?*N = null,
            },
        ) !array.Dense(N) {
            comptime if (meta.isComplex(N))
                @compileError("zsl.array.Dense(N).initLinspace: not defined for complex N, got \n\tN = " ++ @typeName(N) ++ "\n");

            if (opts.num == 0)
                return array.Error.ZeroDimension;

            const arr: array.Dense(N) = try .init(allocator, &.{opts.num}, .col_major);

            if (opts.num == 1) {
                arr.data[0] = start;

                if (opts.retstep) |r|
                    r.* = numeric.zero(N);

                return arr;
            } else if (opts.num == 2) {
                if (opts.endpoint) {
                    arr.data[0] = start;
                    arr.data[1] = stop;
                } else {
                    arr.data[0] = start;
                    numeric.divInto(&arr.data[1], numeric.add(arr.data[0], stop), numeric.two(N));
                }

                if (opts.retstep) |r|
                    numeric.divInto(r, numeric.sub(arr.data[1], arr.data[0]), numeric.two(N));

                return arr;
            }

            var step: N = numeric.sub(stop, start);

            if (opts.endpoint)
                numeric.divInto(&step, step, opts.num - 1)
            else
                numeric.divInto(&step, step, opts.num);

            if (opts.retstep) |r|
                r.* = step;

            if (opts.num == 3 and opts.endpoint) {
                arr.data[0] = start;
                numeric.addInto(&arr.data[1], start, step);
                arr.data[2] = stop;

                return arr;
            } else if (opts.num == 3 and !opts.endpoint) {
                arr.data[0] = start;
                numeric.addInto(&arr.data[1], start, step);
                numeric.subInto(&arr.data[2], stop, step);

                return arr;
            }

            arr.data[0] = start;
            var i: usize = 1;
            while (i < opts.num - 2) : (i += 1) {
                numeric.addInto(&arr.data[i], arr.data[i - 1], step);
            }

            if (opts.endpoint) {
                numeric.addInto(&arr.data[opts.num - 2], arr.data[opts.num - 3], step);
                arr.data[opts.num - 1] = stop;
            } else {
                numeric.addInto(&arr.data[opts.num - 2], arr.data[opts.num - 3], step);
                numeric.subInto(&arr.data[opts.num - 1], stop, step);
            }

            return arr;
        }

        /// Initializes a new `array.Dense(N)` with all elements set as
        /// `opts.num` logarithmically spaced samples, calculated over a given
        /// interval. It is calculated by calling `initLogspace` and applying
        /// `numeric.pow` elementwise with `base`.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `start` (`N`): Start of the interval (inclusive). The actual start
        ///   is `base^start`.
        /// * `stop` (`N`): End of the interval (inclusive if `opts.endpoint`,
        ///   exclusive otherwise). The actual end is `base^end`.
        /// * `opts`: Optional parameters:
        ///   * `num` (`usize = 50`): Number of samples to generate.
        ///   * `endpoint` (`bool = true`): Wether the range is inclusive or
        ///     exclusive.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The newly initialized array.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `array.Error.ZeroDimension`: If `opts.num` is zero.
        pub fn initLogspace(
            allocator: std.mem.Allocator,
            start: N,
            stop: N,
            base: N,
            opts: struct {
                num: usize = 50,
                endpoint: bool = true,
            },
        ) !array.Dense(N) {
            comptime if (meta.isComplex(N))
                @compileError("zsl.array.Dense(N).initLogspace: not defined for complex N, got \n\tN = " ++ @typeName(N) ++ "\n");

            var arr: array.Dense(N) = try .initLinspace(
                allocator,
                start,
                stop,
                .{
                    .num = opts.num,
                    .endpoint = opts.endpoint,
                },
            );

            numeric.powInto(&arr, base, arr) catch unreachable;

            return arr;
        }

        /// Deinitializes the array, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*array.Dense(N)`): A pointer to the array to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *array.Dense(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._size()]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the element from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn get(self: array.Dense(N), index: []const usize) !N {
            try self._checkIndex(index);

            return self.data[self._index(index)];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the element from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn getAssumeInBounds(self: array.Dense(N), index: []const usize) N {
            return self.data[self._index(index)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*array.Dense(N)`): A pointer to the array to set the
        ///   element in.
        /// * `index` (`[]const usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn set(self: *array.Dense(N), index: []const usize, value: N) !void {
            try self._checkIndex(index);

            self.data[self._index(index)] = value;
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*array.Dense(N)`): A pointer to the array to set the
        ///   element in.
        /// * `index` (`[]const usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn put(self: *array.Dense(N), index: []const usize, value: N) void {
            self.data[self._index(index)] = value;
        }

        /// Returns a view of the 2d-array as a general dense matrix. If
        /// `layout` is column major the first stride must be one, otherwise the
        /// last stride must be one.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        /// * `layout` (`comptime Layout`): Specifies the layout of the matrix.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The general dense matrix view.
        ///
        /// ## Errors
        /// * `array.Error.NotConvertible`: If the array does not have 2
        ///   dimensions, or if the strides are incorrect.
        pub fn asGeneralDenseMatrix(self: array.Dense(N), comptime layout: Layout) !matrix.general.Dense(N, layout) {
            if (self.ndim != 2 or
                ((comptime layout == .col_major) and self.strides[0] != 1) or
                ((comptime layout == .row_major) and self.strides[1] != 1))
                return array.Error.NotConvertible;

            return .{
                .data = self.data,
                .rows = self.shape[0],
                .cols = self.shape[1],
                .ld = if (comptime layout == .col_major) self.strides[1] else self.strides[0],
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a view of the array with a different size. Strides are
        /// ordered in ascending order if `layout` is column major, or in
        /// descending order otherwise.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        /// * `layout` (`Layout`): Specifies the layout of the array.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The array view.
        ///
        /// ## Errors
        /// * `array.Error.TooManyDimensions`: If `shape.len` is larger than
        ///   `array.max_dimensions`.
        /// * `array.Error.ZeroDimension`: If any dimension is zero.
        /// * `array.Error.DimensionMismatch`: If the size of both shapes does
        ///   not match.
        pub fn reshape(self: array.Dense(N), shape: []const usize, layout: Layout) !array.Dense(N) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len == 0)
                return array.Error.ZeroDimension;

            var new_size: usize = 1;
            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
            var i: usize = 0;
            while (i < shape.len) : (i += 1) {
                const idx: usize = if (layout == .col_major) i else shape.len - i - 1;

                new_strides[idx] = new_size;
                new_size *= shape[idx];

                new_shape[i] = shape[i];
            }

            if (new_size != self._size())
                return array.Error.DimensionMismatch;

            return .{
                .data = self.data,
                .ndim = shape.len,
                .shape = new_shape,
                .strides = new_strides,
                .flags = .{ .owns_data = false },
            };
        }

        /// Returns a 1d view of the array.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The array view.
        pub fn ravel(self: array.Dense(N)) array.Dense(N) {
            return .{
                .data = self.data,
                .ndim = 1,
                .shape = .{self._size()} ++ .{0} ** (array.max_dimensions - 1),
                .strides = .{1} ++ .{0} ** (array.max_dimensions - 1),
                .flags = .{
                    .owns_data = false,
                },
            };
        }

        pub fn transpose(self: array.Dense(N), axes: ?[]const usize) !array.Dense(N) {
            const axes_: []const usize =
                axes orelse
                array.trivialReversePermutation(self.ndim)[0..self.ndim];

            if (axes_.len == 0)
                return array.Error.ZeroDimension;

            if (axes_.len != self.ndim)
                return array.Error.DimensionMismatch;

            if (!array.isPermutation(self.ndim, axes_))
                return array.Error.InvalidAxes; // axes must be a valid permutation of [0, ..., ndim - 1]

            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
            var size: usize = 1;

            var i: usize = 0;
            while (i < self.ndim) : (i += 1) {
                const idx: usize = axes_[i];

                new_shape[i] = self.shape[idx];
                new_strides[i] = numeric.cast(isize, self.strides[idx]);
                size *= new_shape[i];
            }

            return .{
                .data = self.data,
                .ndim = self.ndim,
                .shape = new_shape,
                .strides = new_strides,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn broadcast(self: array.Dense(N), shape: []const usize) !array.Dense(N) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len < self.ndim)
                return array.Error.TooLittleDimensions;

            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
            var size: usize = 1;

            var i: isize = numeric.cast(isize, shape.len - 1);
            const diff: isize = numeric.cast(isize, shape.len - self.ndim);
            while (i >= 0) : (i -= 1) {
                if (shape[numeric.cast(usize, i)] == 0)
                    return array.Error.ZeroDimension;

                if (i - diff >= 0) {
                    if (self.shape[numeric.cast(usize, i - diff)] != 1 and
                        self.shape[numeric.cast(usize, i - diff)] != shape[numeric.cast(usize, i)])
                        return array.Error.NotBroadcastable; // Broadcasting is not possible if the shapes do not match or are not compatible.

                    new_shape[numeric.cast(usize, i)] = int.max(self.shape[numeric.cast(usize, i - diff)], shape[numeric.cast(usize, i)]);
                    strides[numeric.cast(usize, i)] = self.strides[numeric.cast(usize, i - diff)];
                } else {
                    new_shape[numeric.cast(usize, i)] = shape[numeric.cast(usize, i)];
                    strides[numeric.cast(usize, i)] = 0; // No stride for the new dimensions.
                }

                size *= new_shape[numeric.cast(usize, i)];
            }

            return .{
                .data = self.data,
                .ndim = numeric.cast(usize, shape.len),
                .shape = new_shape,
                .strides = strides,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn slice(self: array.Dense(N), ranges: []const array.Range) !array.Dense(N) {
            if (ranges.len == 0 or ranges.len > self.ndim) {
                return error.DimensionMismatch;
            }

            var ndim: usize = self.ndim;
            var size: usize = 1;
            var shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
            var offset: usize = 0;

            var i: usize = 0;
            var j: usize = 0;
            while (i < self.ndim) {
                const stride: isize = numeric.cast(isize, self.strides[i]);

                if (i >= ranges.len) {
                    shape[j] = self.shape[i];
                    strides[j] = stride;
                    size *= self.shape[i];
                    j += 1;
                    i += 1;
                    continue;
                }

                var range: array.Range = ranges[i];
                if (range.start != int.maxVal(usize) and range.start == range.stop) {
                    return array.Error.InvalidRange;
                } else if (range.step > 0) {
                    if (range.start != int.maxVal(usize) and range.start >= self.shape[i] or
                        (range.stop != int.maxVal(usize) and range.stop > self.shape[i]))
                        return array.Error.RangeOutOfBounds;
                } else if (range.step < 0) {
                    if ((range.stop != int.maxVal(usize) and range.stop >= self.shape[i]) or
                        (range.start != int.maxVal(usize) and range.start > self.shape[i]))
                        return array.Error.RangeOutOfBounds;
                }

                var len_adjustment: usize = 0;
                if (range.step > 0) {
                    if (range.start == int.maxVal(usize)) {
                        range.start = 0;
                    }

                    if (range.stop == int.maxVal(usize)) {
                        range.stop = self.shape[i];
                    }
                } else if (range.step < 0) {
                    if (range.start == int.maxVal(usize)) {
                        range.start = self.shape[i] - 1;
                    }

                    if (range.stop == int.maxVal(usize)) {
                        range.stop = 0;
                        len_adjustment = 1;
                    }
                }

                const len: usize = range.len() + len_adjustment;
                if (len == 1) {
                    ndim -= 1;
                } else {
                    shape[j] = len;
                    strides[j] = stride * range.step;
                    size *= len;
                    j += 1;
                }

                if (stride < 0) {
                    offset -= range.start * numeric.cast(usize, int.abs(stride));
                } else {
                    offset += range.start * numeric.cast(usize, stride);
                }

                i += 1;
            }

            return .{
                .data = self.data + offset,
                .ndim = ndim,
                .shape = shape,
                .strides = strides,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn _size(self: array.Dense(N)) usize {
            var size: usize = 1;
            for (self.shape[0..self.ndim]) |dim| {
                size *= dim;
            }
            return size;
        }

        fn _index(self: array.Dense(N), index: []const usize) usize {
            var idx: isize = 0;
            var i: usize = 0;
            while (i < index.len) : (i += 1) {
                idx += numeric.cast(isize, index[i]) * self.strides[i];
            }

            return numeric.cast(usize, idx);
        }

        fn _checkIndex(self: array.Dense(N), index: []const usize) !void {
            if (index.len > self.ndim)
                return array.Error.DimensionMismatch;

            var i: usize = 0;
            while (i < index.len) : (i += 1) {
                if (index[i] >= self.shape[i]) {
                    return array.Error.PositionOutOfBounds;
                }
            }
        }
    };
}

// fn loop1(
//     result: anytype,
//     x: anytype,
//     comptime op: anytype,
//     comptime @"inline": bool,
//     depth: usize,
//     comptime layout: types.IterationLayout,
//     ir: usize,
//     ix: usize,
//     ctx: anytype,
// ) !void {
//     if (depth == 0) {
//         const opinfo = @typeInfo(@TypeOf(op));
//         const idx: usize = if (comptime layout == .left_to_right) 0 else result.ndim - 1;

//         var jr: usize = ir;
//         var jx: usize = ix;
//         var j: usize = 0;
//         while (j < x.shape[idx]) : (j += 1) {
//             if (comptime !@"inline") {
//                 if (comptime opinfo.@"fn".params.len == 1) {
//                     result.data[jr] = op(x.data[jx]);
//                 } else if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[jr] = try op(x.data[jx], ctx);
//                 }
//             } else {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     op(&result.data[jr], x.data[jx]);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     try op(&result.data[jr], x.data[jx], ctx);
//                 }
//             }

//             jr += result.strides[idx];
//             jx += x.strides[idx];
//         }
//     } else {
//         const idx: usize = if (comptime layout == .left_to_right) depth else result.ndim - depth - 1;

//         var jr: usize = ir;
//         var jx: usize = ix;
//         var j: usize = 0;
//         while (j < x.shape[idx]) : (j += 1) {
//             try loop1(
//                 result,
//                 x,
//                 op,
//                 @"inline",
//                 depth - 1,
//                 layout,
//                 jr,
//                 jx,
//                 ctx,
//             );

//             jr += result.strides[idx];
//             jx += x.strides[idx];
//         }
//     }
// }

// pub fn apply1(
//     allocator: std.mem.Allocator,
//     x: anytype,
//     comptime op: anytype,
//     ctx: anytype,
// ) !Dense(ReturnType1(op, Numeric(@TypeOf(x))), layoutOf(@TypeOf(x))) {
//     const X: type = Numeric(@TypeOf(x));

//     var result: Dense(ReturnType1(op, X), layoutOf(@TypeOf(x))) = try .init(allocator, x.shape[0..x.ndim]);
//     errdefer result.deinit(allocator);

//     if (std.mem.eql(usize, result.strides[0..result.ndim], x.strides[0..x.ndim])) {
//         // Trivial loop
//         //errdefer cleanup(ReturnType1(op, X), allocator, result.data[0..j]);

//         const opinfo = @typeInfo(@TypeOf(op));
//         var i: usize = 0;
//         while (i < result.size) : (i += 1) {
//             if (comptime opinfo.@"fn".params.len == 1) {
//                 result.data[i] = op(x.data[i]);
//             } else if (comptime opinfo.@"fn".params.len == 2) {
//                 result.data[i] = try op(x.data[i], ctx);
//             }
//         }
//     } else {
//         try loop1(
//             &result,
//             &x,
//             op,
//             false,
//             result.ndim - 1,
//             comptime layoutOf(@TypeOf(x)).toIterationLayout(),
//             0,
//             0,
//             ctx,
//         );
//     }

//     return result;
// }

// pub fn apply1_(
//     o: anytype,
//     x: anytype,
//     comptime op_: anytype,
//     ctx: anytype,
// ) !void {
//     const X: type = Numeric(@TypeOf(x));

//     var xx: Dense(X) = undefined;
//     if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
//         if (std.mem.eql(usize, o.strides[0..o.ndim], x.strides[0..x.ndim])) {
//             // Trivial loop
//             const opinfo = @typeInfo(@TypeOf(op_));
//             for (0..o.size) |i| {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     op_(&o.data[i], x.data[i]);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     try op_(&o.data[i], x.data[i], ctx);
//                 }
//             }

//             return;
//         } else {
//             // Different layout, but same shape
//             xx = x;
//         }
//     } else {
//         const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });

//         if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim]))
//             return array.Error.NotBroadcastable;

//         xx = try x.broadcast(bct.shape[0..bct.ndim]);
//     }

//     try loop1(
//         o,
//         &xx,
//         op_,
//         true,
//         o.ndim - 1,
//         comptime layoutOf(@TypeOf(o)).toIterationLayout(),
//         0,
//         0,
//         ctx,
//     );

//     return;
// }

// fn loop2_left(
//     result: anytype,
//     x: anytype,
//     y: anytype,
//     comptime op: anytype,
//     comptime @"inline": bool,
//     depth: usize,
//     comptime layout: types.IterationLayout,
//     ir: usize,
//     iy: usize,
//     ctx: anytype,
// ) !void {
//     if (depth == 0) {
//         const opinfo = @typeInfo(@TypeOf(op));
//         const idx: usize = if (comptime layout == .left_to_right) 0 else result.ndim - 1;

//         var jr: usize = ir;
//         var jy: usize = iy;
//         var j: usize = 0;
//         while (j < y.shape[idx]) : (j += 1) {
//             if (comptime !@"inline") {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[jr] = op(x, y.data[jy]);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     result.data[jr] = try op(x, y.data[jy], ctx);
//                 }
//             } else {
//                 if (comptime opinfo.@"fn".params.len == 3) {
//                     op(&result.data[jr], x, y.data[jy]);
//                 } else if (comptime opinfo.@"fn".params.len == 4) {
//                     try op(&result.data[jr], x, y.data[jy], ctx);
//                 }
//             }

//             jr += result.strides[idx];
//             jy += y.strides[idx];
//         }
//     } else {
//         const idx: usize = if (comptime layout == .left_to_right) depth else result.ndim - depth - 1;

//         var jr: usize = ir;
//         var jy: usize = iy;
//         var j: usize = 0;
//         while (j < y.shape[idx]) : (j += 1) {
//             try loop2_left(
//                 result,
//                 x,
//                 y,
//                 op,
//                 @"inline",
//                 depth - 1,
//                 layout,
//                 jr,
//                 jy,
//                 ctx,
//             );

//             jr += result.strides[idx];
//             jy += y.strides[idx];
//         }
//     }
// }

// fn loop2_right(
//     result: anytype,
//     x: anytype,
//     y: anytype,
//     comptime op: anytype,
//     comptime @"inline": bool,
//     depth: usize,
//     comptime layout: types.IterationLayout,
//     ir: usize,
//     ix: usize,
//     ctx: anytype,
// ) !void {
//     if (depth == 0) {
//         const opinfo = @typeInfo(@TypeOf(op));
//         const idx: usize = if (comptime layout == .left_to_right) 0 else result.ndim - 1;

//         var jr: usize = ir;
//         var jx: usize = ix;
//         var j: usize = 0;
//         while (j < x.shape[idx]) : (j += 1) {
//             if (comptime !@"inline") {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[jr] = op(x.data[jx], y);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     result.data[jr] = try op(x.data[jx], y, ctx);
//                 }
//             } else {
//                 if (comptime opinfo.@"fn".params.len == 3) {
//                     op(&result.data[jr], x.data[jx], y);
//                 } else if (comptime opinfo.@"fn".params.len == 4) {
//                     try op(&result.data[jr], x.data[jx], y, ctx);
//                 }
//             }

//             jr += result.strides[idx];
//             jx += x.strides[idx];
//         }
//     } else {
//         const idx: usize = if (comptime layout == .left_to_right) depth else result.ndim - depth - 1;

//         var jr: usize = ir;
//         var jx: usize = ix;
//         var j: usize = 0;
//         while (j < x.shape[idx]) : (j += 1) {
//             try loop2_right(
//                 result,
//                 x,
//                 y,
//                 op,
//                 @"inline",
//                 depth - 1,
//                 layout,
//                 jr,
//                 jx,
//                 ctx,
//             );

//             jr += result.strides[idx];
//             jx += x.strides[idx];
//         }
//     }
// }

// fn loop2(
//     result: anytype,
//     x: anytype,
//     y: anytype,
//     comptime op: anytype,
//     comptime @"inline": bool,
//     depth: usize,
//     comptime layout: types.IterationLayout,
//     ir: usize,
//     ix: usize,
//     iy: usize,
//     ctx: anytype,
// ) !void {
//     if (depth == 0) {
//         const opinfo = @typeInfo(@TypeOf(op));
//         const idx: usize = if (comptime layout == .left_to_right) 0 else result.ndim - 1;

//         var jr: usize = ir;
//         var jx: usize = ix;
//         var jy: usize = iy;
//         var j: usize = 0;
//         while (j < x.shape[idx]) : (j += 1) {
//             if (comptime !@"inline") {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[jr] = op(x.data[jx], y.data[jy]);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     result.data[jr] = try op(x.data[jx], y.data[jy], ctx);
//                 }
//             } else {
//                 if (comptime opinfo.@"fn".params.len == 3) {
//                     op(&result.data[jr], x.data[jx], y.data[jy]);
//                 } else if (comptime opinfo.@"fn".params.len == 4) {
//                     try op(&result.data[jr], x.data[jx], y.data[jy], ctx);
//                 }
//             }

//             jr += result.strides[idx];
//             jx += x.strides[idx];
//             jy += y.strides[idx];
//         }
//     } else {
//         const idx: usize = if (comptime layout == .left_to_right) depth else result.ndim - depth - 1;

//         var jr: usize = ir;
//         var jx: usize = ix;
//         var jy: usize = iy;
//         var j: usize = 0;
//         while (j < x.shape[idx]) : (j += 1) {
//             try loop2(
//                 result,
//                 x,
//                 y,
//                 op,
//                 @"inline",
//                 depth - 1,
//                 layout,
//                 jr,
//                 jx,
//                 jy,
//                 ctx,
//             );

//             jr += result.strides[idx];
//             jx += x.strides[idx];
//             jy += y.strides[idx];
//         }
//     }
// }

// pub fn apply2(
//     allocator: std.mem.Allocator,
//     x: anytype,
//     y: anytype,
//     comptime op: anytype,
//     ctx: anytype,
// ) !EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
//     const X: type = Numeric(@TypeOf(x));
//     const Y: type = Numeric(@TypeOf(y));

//     if (comptime !types.isDenseArray(@TypeOf(x))) {
//         var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, y.shape[0..y.ndim]);
//         errdefer result.deinit(allocator);

//         if (std.mem.eql(usize, result.strides[0..result.ndim], y.strides[0..y.ndim])) {
//             // Trivial loop
//             //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

//             const opinfo = @typeInfo(@TypeOf(op));

//             var i: usize = 0;
//             while (i < result.size) : (i += 1) {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[i] = op(x, y.data[i]);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     result.data[i] = try op(x, y.data[i], ctx);
//                 }
//             }
//         } else {
//             try loop2_left(
//                 &result,
//                 x,
//                 &y,
//                 op,
//                 false,
//                 result.ndim - 1,
//                 comptime layoutOf(@TypeOf(result)).toIterationLayout(),
//                 0,
//                 0,
//                 ctx,
//             );
//         }

//         return result;
//     } else if (comptime !types.isDenseArray(@TypeOf(y))) {
//         var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, x.shape[0..x.ndim]);
//         errdefer result.deinit(allocator);

//         if (std.mem.eql(usize, result.strides[0..result.ndim], x.strides[0..x.ndim])) {
//             //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

//             const opinfo = @typeInfo(@TypeOf(op));

//             var i: usize = 0;
//             while (i < result.size) : (i += 1) {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[i] = op(x.data[i], y);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     result.data[i] = try op(x.data[i], y, ctx);
//                 }
//             }
//         } else {
//             try loop2_right(
//                 &result,
//                 &x,
//                 y,
//                 op,
//                 false,
//                 result.ndim - 1,
//                 comptime layoutOf(@TypeOf(result)).toIterationLayout(),
//                 0,
//                 0,
//                 ctx,
//             );
//         }

//         return result;
//     }

//     var xx: Dense(X) = undefined;
//     var yy: Dense(Y) = undefined;
//     if (std.mem.eql(usize, x.shape[0..x.ndim], y.shape[0..y.ndim])) {
//         if (std.mem.eql(usize, x.strides[0..x.ndim], y.strides[0..y.ndim])) {
//             // Trivial loop
//             var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, x.shape[0..x.ndim]);
//             errdefer result.deinit(allocator);

//             //errdefer cleanup(ReturnType2(op, X, Y), result.data[0..j], ctx);

//             const opinfo = @typeInfo(@TypeOf(op));

//             var i: usize = 0;
//             while (i < result.size) : (i += 1) {
//                 if (comptime opinfo.@"fn".params.len == 2) {
//                     result.data[i] = op(x.data[i], y.data[i]);
//                 } else if (comptime opinfo.@"fn".params.len == 3) {
//                     result.data[i] = try op(x.data[i], y.data[i], ctx);
//                 }
//             }

//             return result;
//         } else {
//             // Different layout, but same shape
//             xx = x;
//             yy = y;
//         }
//     } else {
//         const bct = try array.broadcastShapes(&.{ x.shape[0..x.ndim], y.shape[0..y.ndim] });
//         xx = try x.broadcast(bct.shape[0..bct.ndim]);
//         yy = try y.broadcast(bct.shape[0..bct.ndim]);
//     }

//     var result: EnsureArray(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y)) = try .init(allocator, xx.shape[0..xx.ndim]);
//     errdefer result.deinit(allocator);

//     try loop2(
//         &result,
//         &xx,
//         &yy,
//         op,
//         false,
//         result.ndim - 1,
//         comptime Layout.resolve3(layoutOf(@TypeOf(x)), layoutOf(@TypeOf(y)), layoutOf(@TypeOf(result))).toIterationLayout(),
//         0,
//         0,
//         0,
//         ctx,
//     );

//     return result;
// }

// pub fn apply2_(
//     o: anytype,
//     x: anytype,
//     y: anytype,
//     comptime op_: anytype,
//     ctx: anytype,
// ) !void {
//     const X: type = Numeric(@TypeOf(x));
//     const Y: type = Numeric(@TypeOf(y));

//     if (comptime !types.isDenseArray(@TypeOf(x))) {
//         if (std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim]) and
//             std.mem.eql(usize, o.strides[0..o.ndim], y.strides[0..y.ndim]))
//         {
//             // Trivial loop
//             const opinfo = @typeInfo(@TypeOf(op_));
//             for (0..o.size) |i| {
//                 if (comptime opinfo.@"fn".params.len == 3) {
//                     op_(&o.data[i], x, y.data[i]);
//                 } else if (comptime opinfo.@"fn".params.len == 4) {
//                     try op_(&o.data[i], x, y.data[i], ctx);
//                 }
//             }

//             return;
//         }

//         var yy: Dense(Y, layoutOf(@TypeOf(y))) = undefined;
//         if (std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim])) {
//             yy = y;
//         } else {
//             const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], y.shape[0..y.ndim] });
//             if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
//                 return array.Error.NotBroadcastable;
//             }

//             yy = try y.broadcast(bct.shape[0..bct.ndim]);
//         }

//         try loop2_left(
//             o,
//             x,
//             &yy,
//             op_,
//             true,
//             o.ndim - 1,
//             comptime layoutOf(@TypeOf(o)).toIterationLayout(),
//             0,
//             0,
//             ctx,
//         );

//         return;
//     } else if (comptime !types.isDenseArray(@TypeOf(y))) {
//         if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
//             std.mem.eql(usize, o.strides[0..o.ndim], x.strides[0..x.ndim]))
//         {
//             // Trivial loop
//             const opinfo = @typeInfo(@TypeOf(op_));
//             for (0..o.size) |i| {
//                 if (comptime opinfo.@"fn".params.len == 3) {
//                     op_(&o.data[i], x.data[i], y);
//                 } else if (comptime opinfo.@"fn".params.len == 4) {
//                     try op_(&o.data[i], x.data[i], y, ctx);
//                 }
//             }

//             return;
//         }

//         // Different layout, but same shape
//         var xx: Dense(X, layoutOf(@TypeOf(x))) = undefined;
//         if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim])) {
//             xx = x;
//         } else {
//             const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim] });
//             if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim])) {
//                 return array.Error.NotBroadcastable;
//             }

//             xx = try x.broadcast(bct.shape[0..bct.ndim]);
//         }

//         try loop2_right(
//             o,
//             &xx,
//             y,
//             op_,
//             true,
//             o.ndim - 1,
//             comptime layoutOf(@TypeOf(o)).toIterationLayout(),
//             0,
//             0,
//             ctx,
//         );

//         return;
//     }

//     var xx: Dense(X) = undefined;
//     var yy: Dense(Y) = undefined;
//     if (std.mem.eql(usize, o.shape[0..o.ndim], x.shape[0..x.ndim]) and
//         std.mem.eql(usize, o.shape[0..o.ndim], y.shape[0..y.ndim]))
//     {
//         if (std.mem.eql(usize, o.strides[0..o.ndim], x.strides[0..x.ndim]) and
//             std.mem.eql(usize, o.strides[0..o.ndim], y.strides[0..y.ndim]))
//         {
//             // Trivial loop
//             const opinfo = @typeInfo(@TypeOf(op_));
//             for (0..o.size) |i| {
//                 if (comptime opinfo.@"fn".params.len == 3) {
//                     op_(&o.data[i], x.data[i], y.data[i]);
//                 } else if (comptime opinfo.@"fn".params.len == 4) {
//                     try op_(&o.data[i], x.data[i], y.data[i], ctx);
//                 }
//             }

//             return;
//         } else {
//             // Different layout, but same shape
//             xx = x;
//             yy = y;
//         }
//     } else {
//         const bct = try array.broadcastShapes(&.{ o.shape[0..o.ndim], x.shape[0..x.ndim], y.shape[0..y.ndim] });

//         if (!std.mem.eql(usize, bct.shape[0..bct.ndim], o.shape[0..o.ndim]))
//             return array.Error.NotBroadcastable;

//         xx = try x.broadcast(bct.shape[0..bct.ndim]);
//         yy = try y.broadcast(bct.shape[0..bct.ndim]);
//     }

//     try loop2(
//         o,
//         &xx,
//         &yy,
//         op_,
//         true,
//         o.ndim - 1,
//         comptime Layout.resolve3(layoutOf(@TypeOf(x)), layoutOf(@TypeOf(y)), layoutOf(@TypeOf(o))).toIterationLayout(),
//         0,
//         0,
//         0,
//         ctx,
//     );

//     return;
// }
