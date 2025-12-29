const std = @import("std");

const types = @import("../../types.zig");
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");

const matrix = @import("../../matrix.zig");
const Flags = matrix.Flags;

/// Sparse hermitian matrix type, represented in either CSC or CSR format,
/// depending on if `order` is column-major or row-major, respectively. Only the
/// upper or lower triangular part of the matrix is stored, depending on the
/// `uplo` parameter.
pub fn Sparse(T: type, uplo: Uplo, order: Order) type {
    if (!types.isNumeric(T) or !types.isComplex(T))
        @compileError("matrix.hermitian.Sparse requires a complex numeric type, got " ++ @typeName(T));

    return struct {
        data: [*]T,
        idx: [*]u32,
        ptr: [*]u32,
        nnz: u32,
        size: u32,
        flags: Flags = .{},

        /// Type signatures
        pub const is_matrix = {};
        pub const is_sparse = {};
        pub const is_hermitian = {};

        /// Numeric type
        pub const Numeric = T;

        pub const empty = Sparse(T, uplo, order){
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .size = 0,
            .flags = .{ .owns_data = false },
        };

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*matrix.hermitian.Sparse(T, uplo, order)`):
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
        pub fn deinit(self: *Sparse(T, uplo, order), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self.nnz]);
                allocator.free(self.idx[0..self.nnz]);
                allocator.free(self.ptr[0 .. self.size + 1]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified position.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.hermitian.Sparse(T, uplo, order)`):
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
        pub fn get(self: *const Sparse(T, uplo, order), r: u32, c: u32) !T {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var rr: u32 = r;
            var cc: u32 = c;
            var noconj: bool = true;
            if (comptime uplo == .upper) {
                if (rr > cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    noconj = false;
                }
            } else {
                if (rr < cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    noconj = false;
                }
            }

            if (comptime order == .col_major) {
                const col_start = self.ptr[cc];
                const col_end = self.ptr[cc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == rr)
                        return if (noconj)
                            self.data[i]
                        else
                            ops.conj(self.data[i], .{}) catch unreachable
                    else if (self.idx[i] > rr)
                        break;
                }
            } else {
                const row_start = self.ptr[rr];
                const row_end = self.ptr[rr + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == cc)
                        return if (noconj)
                            self.data[j]
                        else
                            ops.conj(self.data[j], .{}) catch unreachable
                    else if (self.idx[j] > cc)
                        break;
                }
            }

            return constants.zero(T, .{}) catch unreachable;
        }

        /// Gets the element at the specified position without bounds checking.
        ///
        /// Parameters
        /// ----------
        /// `self` (`*const matrix.hermitian.Sparse(T, uplo, order)`):
        /// A pointer to the matrix to get the element from.
        ///
        /// `r` (`u32`):
        /// The row index of the element to get. Assumed to be within bounds and
        /// on the correct triangular part.
        ///
        /// `c` (`u32`):
        /// The column index of the element to get. Assumed to be within bounds
        /// and on the correct triangular part.
        ///
        /// Returns
        /// -------
        /// `T`:
        /// The element at the specified position.
        pub fn at(self: *Sparse(T, uplo, order), r: u32, c: u32) T {
            if (comptime order == .col_major) {
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
        /// `self` (`*matrix.hermitian.Sparse(T, uplo, order)`):
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
        /// If no existing element is present at the specified position.
        ///
        /// Notes
        /// -----
        /// If the elements are of arbitrary precision type, the existing
        /// element at the position is not deinitialized. The user must ensure
        /// that no memory leaks occur. Additionally, the matrix takes ownership
        /// of `value`.
        pub fn set(self: *Sparse(T, uplo, order), r: u32, c: u32, value: T) !void {
            if (r >= self.size or c >= self.size)
                return matrix.Error.PositionOutOfBounds;

            var rr: u32 = r;
            var cc: u32 = c;
            var conj: bool = false;
            if (comptime uplo == .upper) {
                if (rr > cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    conj = true;
                }
            } else {
                if (rr < cc) {
                    const temp: u32 = rr;
                    rr = cc;
                    cc = temp;
                    conj = true;
                }
            }

            // Find the position to update. If the position does not exist,
            // we return an error.
            if (comptime order == .col_major) {
                const col_start = self.ptr[cc];
                const col_end = self.ptr[cc + 1];

                var i: u32 = col_start;
                while (i < col_end) : (i += 1) {
                    if (self.idx[i] == rr) {
                        self.data[i] = value;
                        if (conj) {
                            ops.conj_(
                                &self.data[i],
                                self.data[i],
                                .{},
                            ) catch unreachable;
                        }
                        return;
                    } else if (self.idx[i] > rr) {
                        break;
                    }
                }
            } else {
                const row_start = self.ptr[rr];
                const row_end = self.ptr[rr + 1];

                var j: u32 = row_start;
                while (j < row_end) : (j += 1) {
                    if (self.idx[j] == cc) {
                        self.data[j] = value;
                        if (conj) {
                            ops.conj_(
                                &self.data[j],
                                self.data[j],
                                .{},
                            ) catch unreachable;
                        }
                        return;
                    } else if (self.idx[j] > cc) {
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
        /// `self` (`*matrix.hermitian.Sparse(T, uplo, order)`):
        /// A pointer to the matrix to set the element in.
        ///
        /// `r` (`u32`):
        /// The row index of the element to set. Assumed to be within bounds and
        /// on the correct triangular part, and that an existing element is
        /// present at the position.
        ///
        /// `c` (`u32`):
        /// The column index of the element to set. Assumed to be within bounds
        /// and on the correct triangular part, and that an existing element is
        /// present at the position.
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
        pub fn put(self: *Sparse(T, uplo, order), r: u32, c: u32, value: T) void {
            if (comptime order == .col_major) {
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
        /// `self` (`*matrix.hermitian.Sparse(T, uplo, order)`):
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
        pub fn cleanup(self: *Sparse(T, uplo, order), ctx: anytype) void {
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
