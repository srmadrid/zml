const std = @import("std");

const meta = @import("../meta.zig");
const int = @import("../int.zig");
const numeric = @import("../numeric.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

const arrutils = @import("utils.zig");

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
        /// are ordered in ascending order if `order` is `f`, or in descending
        /// order otherwise.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `order` (`array.Order`): The stride order of the array.
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
        pub fn init(allocator: std.mem.Allocator, shape: []const usize, order: array.Order) !array.Dense(N) {
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
                const idx: usize = if (order == .f) i else shape.len - i - 1;

                if (shape[idx] == 1)
                    array_strides[idx] = 0
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
        /// ascending order if `order` is `f`, or in descending order otherwise.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `value` (`N`): The value to fill the array with.
        /// * `order` (`array.Order`): The stride order of the array.
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
        pub fn initValue(allocator: std.mem.Allocator, shape: []const usize, value: N, order: array.Order) !array.Dense(N) {
            const arr: array.Dense(N) = try .init(allocator, shape, order);

            const size = arr._size();
            var i: usize = 0;
            while (i < size) : (i += 1) {
                arr.data[i] = value;
            }

            return arr;
        }

        /// Initializes a new `array.Dense(N)` with the specified shape, with
        /// all elements set by calling the specified function with the given
        /// arguments. Strides are ordered in ascending order if `order` is `f`,
        /// or in descending order otherwise.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `shape` (`[]const usize`): The shape of the array.
        /// * `@"fn"` (`anytype`): The function to call to fill the array.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        /// * `order` (`array.Order`): The stride order of the array.
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
        pub fn initFn(allocator: std.mem.Allocator, shape: []const usize, comptime @"fn": anytype, args: anytype, order: array.Order) !array.Dense(N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.array.Dense(N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var arr: array.Dense(N) = try .init(allocator, shape, order);
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
            if (numeric.eq(step, 0) or
                (numeric.lt(stop, start) and positive_step) or
                (numeric.gt(stop, start) and !positive_step))
                return array.Error.InvalidRange;

            const diff = if (positive_step)
                numeric.sub(stop, start)
            else
                numeric.sub(start, stop);

            const len: usize = numeric.cast(usize, numeric.ceil(numeric.abs(numeric.div(diff, step))));
            if (len == 0)
                return array.Error.InvalidRange;

            const arr: array.Dense(N) = try .init(allocator, &.{len}, .f);

            numeric.set(&arr.data[0], start);
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

            const arr: array.Dense(N) = try .init(allocator, &.{opts.num}, .f);

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
                    numeric.divInto(&arr.data[1], numeric.add(arr.data[0], stop), 2);
                }

                if (opts.retstep) |r|
                    numeric.divInto(r, numeric.sub(arr.data[1], arr.data[0]), 2);

                return arr;
            }

            var step = numeric.sub(stop, start);

            if (opts.endpoint)
                numeric.divInto(&step, step, opts.num - 1)
            else
                numeric.divInto(&step, step, opts.num);

            if (opts.retstep) |r|
                numeric.set(r, step);

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

            array.apply2IntoUnchecked(&arr, base, arr);

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
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            return self.getAssumeInBounds(index);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the element from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: array.Dense(N), index: []const usize) N {
            const offset = self._index(index);

            const val = if (offset >= 0)
                self.data[numeric.cast(usize, offset)]
            else
                (self.data - numeric.cast(usize, -offset))[0];

            return if (self.flags.noconj) val else numeric.conj(val);
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
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            self.setAssumeInBounds(index, value);
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
        pub fn setAssumeInBounds(self: *array.Dense(N), index: []const usize, value: N) void {
            const offset = self._index(index);
            const val_eff = if (self.flags.noconj) value else numeric.conj(value);

            if (offset >= 0)
                self.data[numeric.cast(usize, offset)] = val_eff
            else
                (self.data - numeric.cast(usize, -offset))[0] = val_eff;
        }

        /// Returns a view of the 1d-array as a dense vector.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The dense vector view.
        ///
        /// ## Errors
        /// * `array.Error.NotConvertible`: If the array does not have 1
        ///   dimension.
        pub fn vectorView(self: array.Dense(N)) !vector.Dense(N) {
            if (self.ndim != 1)
                return array.Error.NotConvertible;

            return .{
                .data = self.data,
                .len = self.shape[0],
                .inc = self.strides[0],
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns a view of the 2d-array as a general dense matrix. If
        /// `layout` is column major the first stride must be one, otherwise the
        /// last stride must be one. The other stride must be positive.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        /// * `layout` (`comptime matrix.Layout`): Specifies the layout of the
        ///   matrix.
        ///
        /// ## Returns
        /// `matrix.general.Dense(N, layout)`: The general dense matrix view.
        ///
        /// ## Errors
        /// * `array.Error.NotConvertible`: If the array does not have 2
        ///   dimensions, or if the strides are incorrect.
        pub fn matrixView(self: array.Dense(N), comptime layout: matrix.Layout) !matrix.general.Dense(N, layout) {
            if (self.ndim != 2)
                return array.Error.NotConvertible;

            if (comptime layout == .col_major) {
                if (self.strides[0] != 1 or self.strides[1] < self.shape[0])
                    return array.Error.NotConvertible;
            } else {
                if (self.strides[1] != 1 or self.strides[0] < self.shape[1])
                    return array.Error.NotConvertible;
            }

            return .{
                .data = self.data,
                .rows = self.shape[0],
                .cols = self.shape[1],
                .ld = numeric.cast(usize, if (comptime layout == .col_major) self.strides[1] else self.strides[0]),
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns a view of the array with a different size. Strides are
        /// ordered in ascending order if `order` is `f`, or in descending order
        /// otherwise.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        /// * `order` (`array.Order`): The stride order of the view.
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
        pub fn reshapeView(self: array.Dense(N), shape: []const usize, order: array.Order) !array.Dense(N) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (!self.isContiguous(order))
                return array.Error.IncompatibleLayout;

            var new_size: usize = 1;
            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;

            var i: usize = 0;
            while (i < shape.len) : (i += 1) {
                const dim_size = shape[i];
                if (dim_size == 0) return array.Error.ZeroDimension;

                const idx: usize = if (order == .f) i else shape.len - 1 - i;

                new_strides[idx] = numeric.cast(isize, new_size);
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
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns a 1d view of the array.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The array view.
        ///
        /// ## Errors
        /// * `array.Error.IncompatibleLayout`: If the N-D array is not contiguous
        ///   relative to the requested order.
        pub fn ravelView(self: array.Dense(N), order: array.Order) !array.Dense(N) {
            if (self.ndim == 1) {
                var view = self;
                view.flags.owns_data = false;
                return view;
            }

            if (!self.isContiguous(order))
                return array.Error.IncompatibleLayout;

            return .{
                .data = self.data,
                .ndim = 1,
                .shape = .{self._size()} ++ .{0} ** (array.max_dimensions - 1),
                .strides = .{1} ++ .{0} ** (array.max_dimensions - 1),
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns a view of the array with its axes permuted.
        ///
        /// By default, reverses the order of the dimensions. Pass an explicit
        /// slice of dimension indices to perform a custom axis permutation.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to get the view of.
        /// * `axes` (`?[]const usize`): The target permutation of axes. If
        ///   `null`, defaults to `.{ self.ndim - 1, ..., 0 }`.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The transposed array view.
        ///
        /// ## Errors
        /// * `array.Error.DimensionMismatch`: If `axes.len != self.ndim`.
        /// * `array.Error.InvalidAxes`: If `axes` is not a valid permutation of
        ///   the array's dimension indices.
        pub fn transposeView(self: array.Dense(N), axes: ?[]const usize) !array.Dense(N) {
            if (axes) |usr_axes| {
                if (usr_axes.len != self.ndim)
                    return array.Error.DimensionMismatch;

                if (!arrutils.isPermutation(usr_axes))
                    return array.Error.InvalidAxes;
            }

            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;

            var i: usize = 0;
            while (i < self.ndim) : (i += 1) {
                const src_dim: usize = if (axes) |usr_axes| usr_axes[i] else (self.ndim - 1 - i);

                new_shape[i] = self.shape[src_dim];
                new_strides[i] = self.strides[src_dim];
            }

            return .{
                .data = self.data,
                .ndim = self.ndim,
                .shape = new_shape,
                .strides = new_strides,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns a view of the array broadcasted to a new shape.
        ///
        /// Broadcasting allows arrays with different shapes to be combined in
        /// arithmetic operations. Dimensions are aligned from right to left.
        /// Two dimensions are compatible if they are equal, or if one of them
        /// is `1`. Stretched or prepended dimensions are mapped to a stride of
        /// `0`.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to broadcast.
        /// * `shape` (`[]const usize`): The target shape to broadcast into.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The broadcasted array view.
        ///
        /// ## Errors
        /// * `array.Error.TooManyDimensions`: If `shape.len` exceeds
        ///   `array.max_dimensions`.
        /// * `array.Error.TooLittleDimensions`: If `shape.len < self.ndim`.
        /// * `array.Error.ZeroDimension`: If any dimension in the target shape
        ///   is `0`.
        /// * `array.Error.NotBroadcastable`: If the existing dimensions cannot
        ///   be aligned to the target shape according to broadcasting rules.
        pub fn broadcastView(self: array.Dense(N), shape: []const usize) !array.Dense(N) {
            if (shape.len > array.max_dimensions)
                return array.Error.TooManyDimensions;

            if (shape.len < self.ndim)
                return array.Error.TooLittleDimensions;

            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;

            const dim_offset: usize = shape.len - self.ndim;

            var i: usize = 0;
            while (i < shape.len) : (i += 1) {
                const target_dim = shape[i];
                if (target_dim == 0) return array.Error.ZeroDimension;

                new_shape[i] = target_dim;

                if (i < dim_offset) {
                    new_strides[i] = 0;
                    continue;
                }

                const src_idx = i - dim_offset;
                const src_dim = self.shape[src_idx];

                if (src_dim == target_dim) {
                    new_strides[i] = self.strides[src_idx];
                } else if (src_dim == 1) {
                    new_strides[i] = 0;
                } else {
                    return array.Error.NotBroadcastable;
                }
            }

            return .{
                .data = self.data,
                .ndim = shape.len,
                .shape = new_shape,
                .strides = new_strides,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        /// Returns a sub-array view extracted from specified ranges.
        ///
        /// Trailing unspecified dimensions are automatically treated as full
        /// slices (`[:]`). Passing an integer index range (`array.Range.index`)
        /// collapses that dimension, reducing the resulting array's `ndim` by
        /// 1.
        ///
        /// ## Arguments
        /// * `self` (`array.Dense(N)`): The array to slice.
        /// * `ranges` (`[]const array.Range`): A slice of ranges for each
        ///   dimension.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The sliced array view.
        ///
        /// ## Errors
        /// * `array.Error.DimensionMismatch`: If `ranges.len > self.ndim`.
        /// * `array.Error.RangeOutOfBounds`: If any resolved slice exceeds
        ///   dimension boundaries.
        pub fn sliceView(self: array.Dense(N), ranges: []const array.Range) !array.Dense(N) {
            if (ranges.len > self.ndim)
                return array.Error.DimensionMismatch;

            var new_ndim: usize = 0;
            var new_shape: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
            var new_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;

            var ptr_offset: isize = 0;

            var i: usize = 0;
            while (i < self.ndim) : (i += 1) {
                const old_stride = self.strides[i];
                const old_dim = self.shape[i];

                const rng = if (i < ranges.len) ranges[i] else array.Range.all;

                const bounds = try rng.resolve(old_dim);

                ptr_offset += numeric.cast(isize, bounds.start) * old_stride;

                if (rng.is_index) {
                    if (bounds.start >= old_dim)
                        return array.Error.RangeOutOfBounds;

                    continue;
                }

                new_shape[new_ndim] = bounds.len;
                new_strides[new_ndim] = old_stride * rng.step;
                new_ndim += 1;
            }

            const final_ptr = if (ptr_offset >= 0)
                self.data + numeric.cast(usize, ptr_offset)
            else
                self.data - numeric.cast(usize, -ptr_offset);

            return .{
                .data = final_ptr,
                .ndim = new_ndim,
                .shape = new_shape,
                .strides = new_strides,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        pub fn isContiguous(self: array.Dense(N), order: array.Order) bool {
            var expected_stride: isize = 1;
            var i: usize = 0;

            while (i < self.ndim) : (i += 1) {
                const idx: usize = if (order == .f) i else self.ndim - 1 - i;

                if (self.shape[idx] <= 1)
                    continue;

                if (self.strides[idx] != expected_stride) return false;

                expected_stride *= numeric.cast(isize, self.shape[idx]);
            }

            return true;
        }

        pub fn _size(self: array.Dense(N)) usize {
            var size: usize = 1;
            for (self.shape[0..self.ndim]) |dim| {
                size *= dim;
            }
            return size;
        }

        fn _index(self: array.Dense(N), index: []const usize) isize {
            var idx: isize = 0;
            var i: usize = 0;
            while (i < index.len) : (i += 1) {
                idx += numeric.cast(isize, index[i]) * self.strides[i];
            }

            return idx;
        }

        pub fn format(self: *const array.Dense(N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const array.Dense(N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .array = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                array: *const array.Dense(N),

                pub fn format(self: @This(), writer: *std.Io.Writer) !void {
                    const shape = self.array.shape[0..self.array.ndim];

                    try writer.print("zsl.array.Dense({s}) (", .{@typeName(N)});
                    for (0..self.array.ndim) |i| {
                        if (i > 0) try writer.writeAll(" × ");
                        try writer.print("{d}", .{self.array.shape[i]});
                    }
                    try writer.writeAll("):\n\n");

                    return arrutils.format(self, num_fmt, shape, writer);
                }
            };
        }
    };
}
