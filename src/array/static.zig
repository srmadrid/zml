const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

const arrutils = @import("utils.zig");

/// Static `n`-dimensional array type, represented as a contiguous array of
/// elements of type `N`, with shape and stride order fixed at compile time.
pub fn Static(comptime shape_: []const usize, N: type, comptime order: array.Order) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.array.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    if (shape_.len == 0)
        @compileError("zsl.array.Static: shape must not be empty\n");

    if (shape_.len > array.max_dimensions)
        @compileError("zsl.array.Static: shape.len must not exceed array.max_dimensions, got \n\tshape.len = " ++ std.fmt.comptimePrint("{d}", .{shape_.len}) ++ "\n\tarray.max_dimensions = " ++ std.fmt.comptimePrint("{d}", .{array.max_dimensions}) ++ "\n");

    for (shape_) |dim| {
        if (dim == 0)
            @compileError("zsl.array.Static: shape dimensions must not be zero, got \n\tshape = " ++ std.fmt.comptimePrint("{any}", .{shape_}) ++ "\n");
    }

    const ndim_ = shape_.len;

    comptime var shape_arr: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    comptime var strides_arr: [array.max_dimensions]usize = .{0} ** array.max_dimensions;
    comptime var size_: usize = 1;
    for (0..ndim_) |i| {
        const idx: usize = if (order == .f) i else ndim_ - i - 1;

        if (shape_[idx] == 1)
            strides_arr[idx] = 0
        else
            strides_arr[idx] = size_;

        size_ *= shape_[idx];
        shape_arr[i] = shape_[i];
    }

    return struct {
        data: [size]N,

        // Type signatures
        pub const is_array = true;
        pub const is_static = true;
        pub const ndim = ndim_;
        pub const shape = shape_arr;
        pub const strides = strides_arr;
        pub const storage_order = order;
        pub const size = size;

        // Numeric type
        pub const Numeric = N;

        pub const empty: array.Static(shape, N, order) = .{
            .data = undefined,
        };

        pub const init = empty;

        /// Initializes a new `array.Static(shape, N, order)`, with all elements
        /// set to the specified value.
        ///
        /// ## Arguments
        /// * `value` (`N`): The value to fill the array with.
        ///
        /// ## Returns
        /// `array.Static(shape, N, order)`: The newly initialized array.
        pub fn initValue(value: N) array.Static(shape, N, order) {
            return .{ .data = @splat(value) };
        }

        /// Initializes a new `array.Static(shape, N, order)`, with all elements
        /// set from the specified flat values, laid out according to `order`.
        ///
        /// ## Arguments
        /// * `values` (`[size]N`): The values to fill the array with.
        ///
        /// ## Returns
        /// `array.Static(shape, N, order)`: The newly initialized array.
        pub fn initArray(values: [size]N) array.Static(shape, N, order) {
            return .{ .data = values };
        }

        /// Initializes a new `array.Static(shape, N, order)`, with all elements
        /// set by calling the specified function with the given arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the array.
        /// * `args` (`anytype`): A tuple of the arguments to call the
        ///   function with.
        ///
        /// ## Returns
        /// `array.Static(shape, N, order)`: The newly initialized array.
        ///
        /// ## Errors
        /// Propagates any error returned by `@"fn"`.
        pub fn initFn(comptime @"fn": anytype, args: anytype) !array.Static(shape, N, order) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.array.Static(shape, N, order).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var arr: array.Static(shape, N, order) = .init;

            inline for (0..size) |i| {
                arr.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return arr;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`array.Static(shape, N, order)`): The array to get the
        ///   element from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn get(self: array.Static(shape, N, order), index: []const usize) !N {
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            return self.getAssumeInBounds(index);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`array.Statuc(shape, N, order)`): The array to get the
        ///   element from.
        /// * `index` (`[]const usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: array.Static(shape, N, order), index: []const usize) N {
            return self.data[_index(index)];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*array.Static(shape, N, order)`): A pointer to the array
        ///   to set the element in.
        /// * `index` (`[]const usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `array.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn set(self: *array.Static(shape, N, order), index: []const usize, value: N) !void {
            try arrutils.checkIndex(self.shape[0..self.ndim], index);

            self.setAssumeInBounds(index, value);
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*array.Static(shape, N, order)`): A pointer to the array
        ///   to set the element in.
        /// * `index` (`[]const usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *array.Static(shape, N, order), index: []const usize, value: N) void {
            self.data[_index(index)] = value;
        }

        /// Returns an `array.Dense(N)` view aliasing this array's data.
        ///
        /// ## Arguments
        /// * `self` (`*array.Static(shape, N, order)`): A pointer to the
        ///   array to get the view of.
        ///
        /// ## Returns
        /// `array.Dense(N)`: The dense view of the array.
        pub fn denseView(self: *array.Static(shape, N, order)) array.Dense(N) {
            var dense_strides: [array.max_dimensions]isize = .{0} ** array.max_dimensions;
            for (0..ndim_) |i| {
                dense_strides[i] = numeric.cast(isize, strides_arr[i]);
            }

            return .{
                .data = &self.data,
                .ndim = ndim,
                .shape = shape,
                .strides = dense_strides,
                .flags = .{ .owns_data = false },
            };
        }

        fn _index(index: []const usize) usize {
            var idx: usize = 0;
            for (0..index.len) |i| {
                idx += index[i] * strides_arr[i];
            }

            return idx;
        }
    };
}
