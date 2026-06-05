const meta = @import("../meta.zig");

const vector = @import("../vector.zig");

/// Static vector type, represented as a contiguous array of elements of type
/// `N`.
pub fn Static(len_: comptime_int, N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.vector.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [len]N,

        // Type signatures
        pub const is_vector = true;
        pub const is_static = true;
        pub const len = len_;

        // Numeric type
        pub const Numeric = N;

        pub const empty = vector.Static(len, N){
            .data = undefined,
        };

        pub const init = empty;

        /// Initializes a new `vector.Static(len, N)`, with all elements set to
        /// the specified value.
        ///
        /// ## Arguments
        /// * `value` (`N`): The value to fill the vector with.
        ///
        /// ## Returns
        /// `vector.Static(len, N)`: The newly initialized vector.
        pub fn initValue(value: N) vector.Static(len, N) {
            return .{ .data = @splat(value) };
        }

        /// Initializes a new `vector.Static(len, N)`, with all elements set by
        /// calling the specified function with the given arguments.
        ///
        /// ## Arguments
        /// * `@"fn"` (`anytype`): The function to call to fill the vector.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `vector.Static(len, N)`: The newly initialized vector.
        ///
        /// ## Errors
        pub fn initFn(comptime @"fn": anytype, args: anytype) !vector.Static(len, N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.vector.Static(len, N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var vec: vector.Static(len, N) = .init;

            inline for (0..len) |i| {
                vec.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return vec;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`vector.Static(len, N)`): The vector to get the element
        ///   from.
        /// * `index` (`usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn get(self: vector.Static(len, N), index: usize) !N {
            if (index >= len)
                return vector.Error.PositionOutOfBounds;

            return self.data[index];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`vector.Static(len, N)`): The vector to get the element
        ///   from.
        /// * `index` (`usize`): The index of the element to get. Assumed to be
        ///   within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: vector.Static(len, N), index: usize) N {
            return self.data[index];
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Static(len, N)`): A pointer to the vector to set
        ///   the element in.
        /// * `index` (`usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn set(self: *vector.Static(len, N), index: usize, value: N) !void {
            if (index >= len)
                return vector.Error.PositionOutOfBounds;

            self.data[index] = value;
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Static(len, N)`): A pointer to the vector to set
        ///   the element in.
        /// * `index` (`usize`): The index of the element to set. Assumed to be
        ///   within bounds.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *vector.Static(len, N), index: usize, value: N) void {
            self.data[index] = value;
        }

        /// Sets all elements of the vector from an array.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Static(len, N)`): A pointer to the vector to set
        ///   the element in.
        /// * `values` (`[len]N`): The values to set the elements to.
        ///
        /// ## Returns
        /// `void`
        pub fn setArray(self: *vector.Static(len, N), values: [len]N) void {
            self.data = values;
        }
    };
}
