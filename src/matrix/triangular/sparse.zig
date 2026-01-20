const std = @import("std");

const types = @import("../../types.zig");
const Layout = types.Layout;
const Uplo = types.Uplo;
const Diag = types.Diag;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

/// Sparse triangular matrix type, represented in either CSC or CSR format,
/// depending on if `order` is column-major or row-major, respectively. Only the
/// upper or lower triangular part of the matrix is stored, depending on the
/// `uplo` parameter, and the diagonal can be either unit, meaning all diagonal
/// elements are assumed to be 1 and not stored, or non-unit, meaning the
/// diagonal elements are stored normally.
pub fn Sparse(T: type, uplo: Uplo, diag: Diag, layout: Layout) type {
    if (!types.isNumeric(T))
        @compileError("T must be a numeric type");

    return struct {
        data: [*]T,
        idx: [*]u32,
        ptr: [*]u32,
        nnz: u32,
        rows: u32,
        cols: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_matrix = {};
        pub const is_sparse = {};
        pub const is_triangular = {};
        pub const storage_layout = layout;
        pub const storage_uplo = uplo;
        pub const storage_diag = diag;

        /// Numeric type
        pub const Numeric = T;

        pub const empty = Sparse(T, uplo, diag, layout){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangua.Sparse(T, uplo, diag, order)`):
        /// A pointer to the matrix to deinitialize.
        ///
        /// `allocator` (`std.mem.Allocator`):
        /// The allocator to use for memory deallocation. Must be the same
        /// allocator used to compile `self`.
        ///
        /// Returns
        /// -------
        /// `void`
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, `cleanup` must be
        /// called before `deinit` to properly deinitialize the elements.
        pub fn deinit(self: *Sparse(T, uplo, diag, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.nnz]);
                allocator.free(self.idx[0..self.nnz]);
                allocator.free(self.ptr[0..(if (comptime layout == .col_major) self.cols + 1 else self.rows + 1)]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Sparse(T, uplo, diag, order)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        ///
        /// Errors
        /// ------
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` or `c` is out of bounds.
        pub fn get(self: *const Sparse(T, uplo, diag, layout), r: u32, c: u32) !T {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return constants.zero(T, .{}) catch unreachable;
            } else {
                if (r < c)
                    return constants.zero(T, .{}) catch unreachable;
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return constants.one(T, .{}) catch unreachable;
            }

            if (comptime layout == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r)
                        return self.data[i]
                    else if (self.idx[i] > r)
                        break;
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c)
                        return self.data[j]
                    else if (self.idx[j] > c)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.triangular.Sparse(T, uplo, diag, order)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds, on
        /// the correct triangular part, and outside the diagonal if `diag` is
        /// unit.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds,
        /// on the correct triangular part, and outside the diagonal if `diag`
        /// is unit.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub fn at(self: *Sparse(T, uplo, diag, layout), r: u32, c: u32) T {
            if (comptime layout == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r)
                        return self.data[i]
                    else if (self.idx[i] > r)
                        break;
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c)
                        return self.data[j]
                    else if (self.idx[j] > c)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        /// Sets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangular.Sparse(T, uplo, diag, order)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set.
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
        /// `matrix.Error.PositionOutOfBounds`:
        /// If `r` or `c` is out of bounds.
        ///
        /// `matrix.Error.BreaksStructure`:
        /// If no existing element is present at the specified position, if
        /// `r == c` and `diag` is unit, or if the position `(r, c)` is outside
        /// the correct triangular part of the matrix.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Sparse(T, uplo, diag, layout), r: u32, c: u32, value: T) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime uplo == .upper) {
                if (r > c)
                    return matrix.Error.BreaksStructure;
            } else {
                if (r < c)
                    return matrix.Error.BreaksStructure;
            }

            if (comptime diag == .unit) {
                if (r == c)
                    return matrix.Error.BreaksStructure;
            }

            // Find the position to update. If the position does not exist,
            // we return an error.
            if (comptime layout == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r) {
                        self.data[i] = value;
                        return;
                    } else if (self.idx[i] > r) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c) {
                        self.data[j] = value;
                        return;
                    } else if (self.idx[j] > c) {
                        break;
                    }
                }
            }

            return matrix.Error.BreaksStructure;
        }

        /// Sets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangular.Sparse(T, uplo, diag, order)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds, on
        /// the correct triangular part, outside the diagonal if `diag` is unit,
        /// and that an existing element is present at the position.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds,
        /// on the correct triangular part, outside the diagonal if `diag` is
        /// unit, and that an existing element is present at the position.
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
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn put(self: *Sparse(T, uplo, diag, layout), r: u32, c: u32, value: T) void {
            if (comptime layout == .col_major) {
                const col_start = self.ptr[c];
                const col_end = self.ptr[c + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == r) {
                        self.data[i] = value;
                        return;
                    } else if (self.idx[i] > r) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[r];
                const row_end = self.ptr[r + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == c) {
                        self.data[j] = value;
                        return;
                    } else if (self.idx[j] > c) {
                        break;
                    }
                }
            }

            return;
        }

        /// Cleans up the elements of the matrix, deinitializing them if
        /// necessary.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.triangular.Sparse(T, uplo, diag, order)`):
        /// A pointer to the matrix to clean up.
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
        pub fn cleanup(self: *Sparse(T, uplo, diag, layout), ctx: anytype) void {
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
                    while (i < self.nnz) : (i += 1) {
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
