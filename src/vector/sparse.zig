const std = @import("std");

const types = @import("../types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const ops = @import("../ops.zig");
const constants = @import("../constants.zig");
const int = @import("../int.zig");

const vector = @import("../vector.zig");
const Flags = vector.Flags;

/// Sparse vector type, represented as a contiguous array of non-zero elements
/// of type `T` along with their corresponding indices, in ascending order.
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

        /// Initializes a new sparse vector with the specified length and an
        /// initial capacity for non-zero elements.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `len` (`u32`):
        /// The length of the vector.
        ///
        /// `nnz` (`u32`):
        /// The initial capacity for non-zero elements.
        ///
        /// Returns
        /// -------
        /// `vector.Sparse(T)`:
        /// The newly initialized vector.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `vector.Error.ZeroLength`:
        /// If `len` is zero.
        ///
        /// `vector.Error.DimensionMismatch`:
        /// If `nnz` is zero or greater than `len`.
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

        /// Deinitializes the vector, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Sparse(T)`):
        /// A pointer to the vector to deinitialize.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory deallocation. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, `cleanup` must be
        /// called before `deinit` to properly deinitialize the elements.
        pub fn deinit(self: *Sparse(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.idx[0..self._ilen]);
            }

            self.* = undefined;
        }

        /// Reserves space for at least `new_nnz` non-zero elements.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Sparse(T)`):
        /// A pointer to the vector.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `new_nnz` (`u32`):
        /// The new capacity for non-zero elements.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `vector.Error.DimensionMismatch`:
        /// If `new_nnz` is greater than the length of the vector.
        ///
        /// Notes
        /// -----
        /// If `self` does not own its data or if `new_nnz` is less than or
        /// equal to the current capacity, this function does nothing.
        pub fn reserve(self: *Sparse(T), allocator: std.mem.Allocator, new_nnz: u32) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._ilen)
                return;

            if (new_nnz > self.len)
                return vector.Error.DimensionMismatch;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._ilen) {
                self.idx = (try allocator.realloc(self.idx[0..self._ilen], new_nnz)).ptr;
                self._ilen = new_nnz;
            }
        }

        /// Gets the element at the specified index.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const vector.Sparse(T)`):
        /// A pointer to the vector to get the element from.
        ///
        /// `index` (`u32`):
        /// The index of the element to get.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified index.
        ///
        /// Errors
        /// ------
        /// `vector.Error.PositionOutOfBounds`:
        /// If `index` is out of bounds.
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

        /// Gets the element at the specified index without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const vector.Sparse(T)`):
        /// A pointer to the vector to get the element from.
        ///
        /// `index` (`u32`):
        /// The index of the element to get. Assumed to be valid.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified index.
        pub fn at(self: *Sparse(T), index: u32) T {
            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index)
                    return self.data[i]
                else if (self.idx[i] > index)
                    break;
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        /// Sets the element at the specified index, inserting it if it does not
        /// already exist and shifting elements as necessary to maintain index
        /// order.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Sparse(T)`):
        /// A pointer to the vector to set the element in.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `index` (`u32`):
        /// The index of the element to set.
        ///
        /// `value` (`T`):
        /// The value to set the element to.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails when inserting a new element.
        ///
        /// `vector.Error.PositionOutOfBounds`:
        /// If `index` is out of bounds.
        ///
        /// `vector.Error.DataNotOwned`:
        /// If the vector does not own its data and a resize is required.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type and an existing
        /// element is being overwritten at `index`, the existing is not
        /// deinitialized. The user must ensure that no memory leaks occur.
        ///
        /// If the elements are of arbitrary precision type, the vector takes
        /// ownership of `value`.
        pub fn set(self: *Sparse(T), allocator: std.mem.Allocator, index: u32, value: T) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

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
                if (!self.flags.owns_data)
                    return vector.Error.DataNotOwned;

                // Need more space
                var new_nnz = if (self.nnz * 2 > self.len) self.len else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

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

        /// Sets the element at the specified index without bounds or space
        /// checking, inserting it if it does not already exist and shifting
        /// elements as necessary to maintain index order.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Sparse(T)`):
        /// A pointer to the vector to set the element in. If needed, assumed
        /// to have enough space to insert a new element.
        ///
        /// `index` (`u32`):
        /// The index of the element to set. Assumed to be valid.
        ///
        /// `value` (`T`):
        /// The value to set the element to.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// Attempting to set an index that does not exist or when there is no
        /// space will result in undefined behavior.
        ///
        /// If the elements are of arbitrary precision type and an existing
        /// element is being overwritten at `index`, the existing is not
        /// deinitialized. The user must ensure that no memory leaks occur.
        ///
        /// If the elements are of arbitrary precision type, the vector takes
        /// ownership of `value`.
        pub fn put(self: *Sparse(T), index: u32, value: T) void {
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

        /// Accumulates the specified value at the given index, adding to the
        /// existing value if it exists, or inserting it if it does not.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Sparse(T)`):
        /// A pointer to the vector to set the element in.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations. Must be the same
        /// allocator used to initialize `self`.
        ///
        /// `index` (`u32`):
        /// The index of the element to accumulate at.
        ///
        /// `value` (`T`):
        /// The value to accumulate at the specified index.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If  the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails when inserting a new element.
        ///
        /// `vector.Error.PositionOutOfBounds`:
        /// If `index` is out of bounds.
        ///
        /// `vector.Error.DataNotOwned`:
        /// If the vector does not own its data and a resize is required.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type and the element is
        /// being inserted, the vector takes ownership of `value`.
        ///
        /// When `T` is of arbitrary precision, the context may provide an
        /// optional pre-allocated buffer to store intermediate results of the
        /// addition, avoiding repeated allocations in scenarios where
        /// `accumulate` is called multiple times. If no buffer is provided, the
        /// operation will allocate a temporary buffer internally, using the
        /// allocator specified in the context.
        pub fn accumulate(self: *Sparse(T), allocator: std.mem.Allocator, index: u32, value: T, ctx: anytype) !void {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                            .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                            .buffer = .{ .type = ?*T, .required = false, .default = null },
                        },
                    );
                },
            };

            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            var i: u32 = 0;
            while (i < self.nnz) : (i += 1) {
                if (self.idx[i] == index) {
                    try ops.add_(
                        &self.data[i],
                        self.data[i],
                        value,
                        types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                    );

                    return;
                } else if (self.idx[i] > index) {
                    break;
                }
            }

            if (self.nnz == self._dlen or self.nnz == self._ilen) {
                if (!self.flags.owns_data)
                    return;

                // Need more space
                var new_nnz = if (self.nnz * 2 > self.len) self.len else self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

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

        /// Cleans up the elements of the vector, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Sparse(T)`):
        /// A pointer to the vector to clean up.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T`. If the
        /// context is missing required fields or contains unnecessary or
        /// wrongly typed fields, the compiler will emit a detailed error
        /// message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// This function must be called before `deinit` if the elements are of
        /// arbitrary precision type to properly deinitialize them.
        pub fn cleanup(self: *Sparse(T), ctx: anytype) void {
            return _cleanup(self, self.nnz, ctx);
        }

        pub fn _cleanup(self: *Sparse(T), num_elems: u32, ctx: anytype) void {
            switch (comptime types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    comptime types.validateContext(@TypeOf(ctx), .{});

                    // No cleanup needed for fixed precision types.
                },
                .integer, .rational, .real, .complex => {
                    comptime types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );

                    var i: u32 = 0;
                    while (i < num_elems) : (i += 1) {
                        ops.deinit(
                            &self.data[i],
                            types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                        );
                    }
                },
            }
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
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    if (comptime !types.isSparseVector(@TypeOf(x))) {
        var result: Sparse(R) = try .init(allocator, y.len, y.nnz);
        errdefer result.deinit(allocator);

        var i: u32 = 0;

        errdefer result.cleanup(
            types.renameStructFields(
                types.keepStructFields(
                    ctx,
                    &.{"allocator"},
                ),
                .{ .allocator = "element_allocator" },
            ),
        );

        const opinfo = @typeInfo(@TypeOf(op));
        while (i < y.nnz) : (i += 1) {
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
        var result: Sparse(R) = try .init(allocator, x.len, x.nnz);
        errdefer result.deinit(allocator);

        var i: u32 = 0;

        errdefer result.cleanup(
            types.renameStructFields(
                types.keepStructFields(
                    ctx,
                    &.{"allocator"},
                ),
                .{ .allocator = "element_allocator" },
            ),
        );

        const opinfo = @typeInfo(@TypeOf(op));
        while (i < x.nnz) : (i += 1) {
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

    var result: Sparse(R) = try .init(allocator, x.len, int.min(x.nnz + y.nnz, x.len));
    errdefer result.deinit(allocator);

    var i: u32 = 0;
    var j: u32 = 0;

    errdefer result.cleanup(
        types.renameStructFields(
            types.keepStructFields(
                ctx,
                &.{"allocator"},
            ),
            .{ .allocator = "element_allocator" },
        ),
    );

    const opinfo = @typeInfo(@TypeOf(op));
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
