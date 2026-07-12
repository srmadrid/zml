const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");

const int = @import("../int.zig");

const vecutils = @import("utils.zig");

/// Dense vector type, represented as a contiguous array of elements of type
/// `N` and a stride.
pub fn Dense(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.vector.Dense: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        len: usize,
        inc: isize,
        flags: vector.Flags,

        // Type signatures
        pub const is_vector = true;
        pub const is_dense = true;

        // Numeric type
        pub const Numeric = N;

        pub const empty = vector.Dense(N){
            .data = &.{},
            .len = 0,
            .inc = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `vector.Dense(N)` with the specified length.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `len` (`usize`): The length of the vector.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.ZeroLength`: If `len` is zero.
        pub fn init(allocator: std.mem.Allocator, len: usize) !vector.Dense(N) {
            if (len == 0)
                return vector.Error.ZeroLength;

            return .{
                .data = (try allocator.alloc(N, len)).ptr,
                .len = len,
                .inc = 1,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new `vector.Dense(N)` with the given buffer.
        ///
        /// ## Arguments
        /// * `buffer` (`[]N`): The buffer.
        /// * `inc` (`isize`): The step size between elements in memory.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `vector.Error.ZeroLength`: If if the length of the buffer is zero.
        pub fn initBuffer(buffer: []N, inc: isize) !vector.Dense(N) {
            if (buffer.len == 0)
                return vector.Error.ZeroLength;

            return .{
                .data = buffer.ptr,
                .len = int.div(buffer.len, numeric.cast(usize, int.abs(inc))),
                .inc = inc,
                .flags = .{ .owns_data = false },
            };
        }

        /// Initializes a new `vector.Dense(N)` with the specified length, with
        /// all elements set to the specified value.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `len` (`usize`): The length of the vector.
        /// * `value` (`N`): The value to fill the vector with.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.ZeroLength`: If `len` is zero.
        pub fn initValue(allocator: std.mem.Allocator, len: usize, value: N) !vector.Dense(N) {
            const vec: vector.Dense(N) = try .init(allocator, len);

            var i: usize = 0;
            while (i < len) : (i += 1) {
                vec.data[i] = value;
            }

            return vec;
        }

        /// Initializes a new `vector.Dense(N)` with the specified length, with
        /// all elements set by calling the specified function with the given
        /// arguments.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `len` (`usize`): The length of the vector.
        /// * `@"fn"` (`anytype`): The function to call to fill the vector.
        /// * `args` (`anytype`): A tuple of the arguments to call the function
        ///   with.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The newly initialized vector.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `vector.Error.ZeroLength`: If `len` is zero.
        pub fn initFn(allocator: std.mem.Allocator, len: usize, comptime @"fn": anytype, args: anytype) !vector.Dense(N) {
            const Fn = @TypeOf(@"fn");
            const Args = @TypeOf(args);

            const fn_info = @typeInfo(Fn);
            const args_info = @typeInfo(Args);

            comptime if (fn_info != .@"fn" or args_info != .@"struct")
                @compileError("zsl.vector.Dense(N).initFn: @\"fn\" must be a function and args must be a struct, got \n\t@\"fn\": " ++ @typeName(Fn) ++ "\n\targs: " ++ @typeName(Args) ++ "\n");

            var vec: vector.Dense(N) = try .init(allocator, len);
            errdefer vec.deinit(allocator);

            var i: usize = 0;
            while (i < len) : (i += 1) {
                vec.data[i] = if (comptime @typeInfo(meta.ReturnTypeFromInputs(@"fn", &meta.structToArrayOfTypes(Args))) == .error_union)
                    try @call(.auto, @"fn", args)
                else
                    @call(.auto, @"fn", args);
            }

            return vec;
        }

        /// Deinitializes the vector, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Dense(N)`): A pointer to the vector to
        ///   deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *vector.Dense(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0 .. self.len * numeric.cast(usize, int.abs(self.inc))]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`vector.Dense(N)`): The vector to get the element from.
        /// * `index` (`usize`): The index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn get(self: vector.Dense(N), index: usize) !N {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            return self.getAssumeInBounds(index);
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`vector.Dense(N)`): The vector to get the element from.
        /// * `index` (`usize`): The index of the element to get. Assumed to be
        ///   within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        pub fn getAssumeInBounds(self: vector.Dense(N), index: usize) N {
            return if (self.flags.noconj)
                self.data[self._index(index)]
            else
                numeric.conj(self.data[self._index(index)]);
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Dense(N)`): A pointer to the vector to set the
        ///   element in.
        /// * `index` (`usize`): The index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `vector.Error.PositionOutOfBounds`: If `index` is out of bounds.
        pub fn set(self: *vector.Dense(N), index: usize, value: N) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            return self.setAssumeInBounds(index, value);
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`*vector.Dense(N)`): A pointer to the vector to set the
        ///   element in.
        /// * `index` (`usize`): The index of the element to set. Assumed to be
        ///   within bounds.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        pub fn setAssumeInBounds(self: *vector.Dense(N), index: usize, value: N) void {
            self.data[self._index(index)] = if (self.flags.noconj) value else numeric.conj(value);
        }

        /// Creates a copy of the vector.
        ///
        /// ## Arguments
        /// * `self` (`vector.Dense(N)`): The vector to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `vector.Dense(N)`: The copied vector.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn clone(self: vector.Dense(N), allocator: std.mem.Allocator) !vector.Dense(N) {
            var vec: vector.Dense(N) = try .init(allocator, self.len);

            var i: usize = 0;
            while (i < vec.len) : (i += 1) {
                vec.data[i] = self.getAssumeInBounds(i);
            }

            return vec;
        }

        /// Views the vector as a diagonal matrix.
        ///
        /// ## Arguments
        /// * `self` (`vector.Dense(N)`): A pointer to the vector to view as a
        ///   diagonal matrix.
        /// * `rows` (`usize`): The number of rows of the resulting matrix.
        /// * `cols` (`usize`): The number of columns of the resulting matrix.
        ///
        /// ## Returns
        /// `matrix.Diagonal(N)`: The diagonal matrix view of the vector.
        ///
        /// ## Errors
        /// * `matrix.Error.ZeroDimension`: If `rows` or `cols` is zero.
        /// * `vector.Error.DimensionMismatch`: If both `rows` and `cols` are
        ///   greater than the length of the vector.
        ///*  `vector.Error.NonContiguousData`: If the vector data is not
        ///   contiguous (`inc != 1`).
        pub fn asDiagonal(self: vector.Dense(N), rows: usize, cols: usize) !matrix.diagonal.Sparse(N) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (rows > self.len and cols > self.len)
                return vector.Error.DimensionMismatch;

            if (self.inc != 1)
                return vector.Error.NonContiguousData;

            return .{
                .data = self.data,
                .rows = rows,
                .cols = cols,
                .flags = .{ .owns_data = false, .noconj = self.flags.noconj },
            };
        }

        pub fn _index(self: Dense(N), index: usize) usize {
            return if (self.inc > 0)
                index * numeric.cast(usize, self.inc)
            else
                numeric.cast(usize, (numeric.cast(isize, index) - numeric.cast(isize, self.len) + 1) * self.inc);
        }

        pub fn format(self: vector.Dense(N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const vector.Dense(N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .vector = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                vector: *const vector.Dense(N),

                pub fn format(self: vector.Dense(N).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const len = self.vector.len;

                    try writer.print("zsl.vector.Dense({s}) ({d}):\n\n", .{ @typeName(N), len });

                    return vecutils.format(self, num_fmt, len, writer);
                }
            };
        }
    };
}
