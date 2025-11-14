const std = @import("std");

const types = @import("../types.zig");
const ReturnType2 = types.ReturnType2;
const Numeric = types.Numeric;
const ops = @import("../ops.zig");
const int = @import("../int.zig");

const vector = @import("../vector.zig");
const Flags = vector.Flags;
const matrix = @import("../matrix.zig");

/// A dense vector type.
pub fn Dense(T: type) type {
    if (!types.isNumeric(T))
        @compileError("vector.Dense requires a numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        len: u32,
        inc: i32,
        flags: Flags = .{},

        pub const empty = Dense(T){
            .data = &.{},
            .len = 0,
            .inc = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `vector.Dense(T)` with the specified length.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `len` (`u32`):
        /// The length of the vector.
        ///
        /// Returns
        /// -------
        /// `vector.Dense(T)`:
        /// The newly initialized `vector.Dense(T)`.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `vector.Error.ZeroLength`:
        /// If `len` is zero.
        pub fn init(allocator: std.mem.Allocator, len: u32) !Dense(T) {
            if (len == 0)
                return vector.Error.ZeroLength;

            return .{
                .data = (try allocator.alloc(T, len)).ptr,
                .len = len,
                .inc = 1,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new `vector.Dense(T)` with the specified length,
        /// filled with the specified value.
        ///
        /// Parameters
        /// ----------
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory allocations.
        ///
        /// `len` (`u32`):
        /// The length of the vector.
        ///
        /// `value` (`anytype`):
        /// The value to fill the vector with.
        ///
        /// `ctx` (`anytype`):
        /// A context struct providing necessary resources and configuration for
        /// the operation. The required fields depend on the type `T` and the
        /// type of `value`. If  the context is missing required fields or
        /// contains unnecessary or wrongly typed fields, the compiler will emit
        /// a detailed error message describing the expected structure.
        ///
        /// Returns
        /// -------
        /// `vector.Dense(T)`:
        /// The newly initialized `vector.Dense(T)`.
        ///
        /// Errors
        /// ------
        /// `std.mem.Allocator.Error.OutOfMemory`:
        /// If memory allocation fails.
        ///
        /// `vector.Error.ZeroLength`:
        /// If `len` is zero.
        ///
        /// Notes
        /// -----
        /// The vector does not take ownership of `value` if it is an arbitrary
        /// precision type.
        pub fn full(
            allocator: std.mem.Allocator,
            len: u32,
            value: anytype,
            ctx: anytype,
        ) !Dense(T) {
            comptime switch (types.numericType(T)) {
                .bool, .int, .float, .cfloat => {
                    types.validateContext(@TypeOf(ctx), .{});
                },
                .integer, .rational, .real, .complex => {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                },
            };

            var vec: Dense(T) = try .init(allocator, len);
            errdefer vec.deinit(allocator);

            const value_casted: T = types.scast(T, value);

            var i: u32 = 0;

            errdefer vec._cleanup(i, ctx);

            while (i < len) : (i += 1) {
                vec.data[i] = value_casted;
            }

            return vec;
        }

        /// Deinitializes the vector, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Dense(T)`):
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
        pub fn deinit(self: *Dense(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.len]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const vector.Dense(T)`):
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
        pub fn get(self: *const Dense(T), index: u32) !T {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))];
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const vector.Dense(T)`):
        /// A pointer to the vector to get the element from.
        ///
        /// `index` (`u32`):
        /// The index of the element to get. Assumed to be valid.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified index.
        pub inline fn at(self: *const Dense(T), index: u32) T {
            return if (self.inc > 0)
                self.data[index * types.scast(u32, self.inc)]
            else
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))];
        }

        /// Sets the element at the specified index.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Dense(T)`):
        /// A pointer to the vector to set the element in.
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
        /// `vector.Error.PositionOutOfBounds`:
        /// If `index` is out of bounds.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at `index` is not deinitialized. The user must ensure that
        /// no memory leaks occur. Additionally, the vector takes ownership of
        /// `value`.
        pub fn set(self: *Dense(T), index: u32, value: T) !void {
            if (index >= self.len)
                return vector.Error.PositionOutOfBounds;

            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        /// Sets the element at the specified index without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Dense(T)`):
        /// A pointer to the vector to set the element in.
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
        /// If the elements are of arbitrary precision type, the existing
        /// element at `index` is not deinitialized. The user must ensure that
        /// no memory leaks occur. Additionally, the vector takes ownership of
        /// `value`.
        pub inline fn put(self: *Dense(T), index: u32, value: T) void {
            if (self.inc > 0) {
                self.data[index * types.scast(u32, self.inc)] = value;
            } else {
                self.data[types.scast(u32, (-types.scast(i32, self.len) + 1) * self.inc) - index * types.scast(u32, int.abs(self.inc))] = value;
            }
        }

        /// Views the vector as a diagonal matrix.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const vector.Dense(T)`):
        /// A pointer to the vector to view as a diagonal matrix.
        ///
        /// `rows` (`u32`):
        /// The number of rows of the resulting matrix.
        ///
        /// `cols` (`u32`):
        /// The number of columns of the resulting matrix.
        ///
        /// Returns
        /// -------
        /// `matrix.Diagonal(T)`:
        /// The diagonal matrix view of the vector.
        ///
        /// Errors
        /// ------
        /// `vector.Error.ZeroDimension`:
        /// If `rows` or `cols` is zero.
        ///
        /// `vector.Error.DimensionMismatch`:
        /// If both `rows` and `cols` are greater than the length of the vector.
        ///
        /// `vector.Error.NonContiguousData`:
        /// If the vector data is not contiguous (inc != 1).
        pub fn asDiagonal(self: *const Dense(T), rows: u32, cols: u32) !matrix.Diagonal(T) {
            if (rows == 0 or cols == 0)
                return vector.Error.ZeroDimension;

            if (rows > self.len and cols > self.len)
                return vector.Error.DimensionMismatch;

            if (self.inc != 1)
                return vector.Error.NonContiguousData;

            return .{
                .data = self.data,
                .rows = rows,
                .cols = cols,
                .flags = .{ .owns_data = false },
            };
        }

        /// Cleans up the elements of the vector, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*vector.Dense(T)`):
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
        pub fn cleanup(self: *Dense(T), ctx: anytype) void {
            return _cleanup(self, self.len, ctx);
        }

        pub fn _cleanup(self: *Dense(T), num_elems: u32, ctx: anytype) void {
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

                    if (self.inc == 1) {
                        var i: u32 = 0;
                        while (i < num_elems) : (i += 1) {
                            ops.deinit(
                                &self.data[i],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );
                        }
                    } else {
                        var is: i32 = if (self.inc < 0) (-types.scast(i32, self.len) + 1) * self.inc else 0;
                        var i: u32 = 0;
                        while (i < num_elems) : (i += 1) {
                            ops.deinit(
                                &self.data[types.scast(u32, is)],
                                types.renameStructFields(ctx, .{ .element_allocator = "allocator" }),
                            );

                            is += self.inc;
                        }
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
) !Dense(ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const R: type = ReturnType2(op, Numeric(X), Numeric(Y));

    if (comptime !types.isDenseVector(@TypeOf(x))) {
        var result: Dense(R) = try .init(allocator, y.len);
        errdefer result.deinit(allocator);

        var i: u32 = 0;

        errdefer result._cleanup(
            i,
            types.renameStructFields(
                types.keepStructFields(
                    ctx,
                    &.{"allocator"},
                ),
                .{ .allocator = "element_allocator" },
            ),
        );

        const opinfo = @typeInfo(@TypeOf(op));
        if (y.inc == 1) {
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[i]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[i], ctx);
                }
            }
        } else {
            var iy: i32 = if (y.inc < 0) (-types.scast(i32, y.len) + 1) * y.inc else 0;
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x, y.data[types.scast(u32, iy)]);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x, y.data[types.scast(u32, iy)], ctx);
                }

                iy += y.inc;
            }
        }

        return result;
    } else if (comptime !types.isDenseVector(@TypeOf(y))) {
        var result: Dense(R) = try .init(allocator, x.len);
        errdefer result.deinit(allocator);

        var i: u32 = 0;

        errdefer result._cleanup(
            i,
            types.renameStructFields(
                types.keepStructFields(
                    ctx,
                    &.{"allocator"},
                ),
                .{ .allocator = "element_allocator" },
            ),
        );

        const opinfo = @typeInfo(@TypeOf(op));
        if (x.inc == 1) {
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[i], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[i], y, ctx);
                }
            }
        } else {
            var ix: i32 = if (x.inc < 0) (-types.scast(i32, x.len) + 1) * x.inc else 0;
            while (i < result.len) : (i += 1) {
                if (comptime opinfo.@"fn".params.len == 2) {
                    result.data[i] = op(x.data[types.scast(u32, ix)], y);
                } else if (comptime opinfo.@"fn".params.len == 3) {
                    result.data[i] = try op(x.data[types.scast(u32, ix)], y, ctx);
                }

                ix += x.inc;
            }
        }

        return result;
    }

    if (x.len != y.len)
        return vector.Error.DimensionMismatch;

    var result: Dense(R) = try .init(allocator, x.len);
    errdefer result.deinit(allocator);

    var i: u32 = 0;

    errdefer result._cleanup(
        i,
        types.renameStructFields(
            types.keepStructFields(
                ctx,
                &.{"allocator"},
            ),
            .{ .allocator = "element_allocator" },
        ),
    );

    const opinfo = @typeInfo(@TypeOf(op));
    if (x.inc == 1 and y.inc == 1) {
        while (i < result.len) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[i], y.data[i]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[i], y.data[i], ctx);
            }
        }
    } else {
        var ix: i32 = if (x.inc < 0) (-types.scast(i32, x.len) + 1) * x.inc else 0;
        var iy: i32 = if (y.inc < 0) (-types.scast(i32, y.len) + 1) * y.inc else 0;
        while (i < result.len) : (i += 1) {
            if (comptime opinfo.@"fn".params.len == 2) {
                result.data[i] = op(x.data[types.scast(u32, ix)], y.data[types.scast(u32, iy)]);
            } else if (comptime opinfo.@"fn".params.len == 3) {
                result.data[i] = try op(x.data[types.scast(u32, ix)], y.data[types.scast(u32, iy)], ctx);
            }

            ix += x.inc;
            iy += y.inc;
        }
    }

    return result;
}
