const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const matrix = @import("../../matrix.zig");

/// Sparse symmetric matrix type, represented in either CSC or CSR format,
/// depending on if `layout` is column-major or row-major, respectively. Only
/// the upper or lower triangular part of the matrix is stored, depending on the
/// `uplo` parameter.
pub fn Sparse(N: type, uplo: matrix.Uplo, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.symmetric.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        idx: [*]usize,
        ptr: [*]usize,
        nnz: usize,
        rows: usize,
        cols: usize,
        _dlen: usize,
        _ilen: usize,
        flags: matrix.Flags,

        // Type signatures
        pub const is_matrix = true;
        pub const is_sparse = true;
        pub const is_symmetric = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = uplo;
        pub const storage_diag: ?matrix.Diag = null;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.symmetric.Sparse(N, uplo, layout) = .{
            .data = &.{},
            .idx = &.{},
            .ptr = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            ._dlen = 0,
            ._ilen = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.symmetric.Sparse(N, uplo, layout)` with
        /// the specified size, and a capacity for non-zero elements.
        ///
        /// Initialization with this function is only meant for pre-allocating
        /// an empty matrix to serve as the output destination for in-place
        /// mathematical functions.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the matrix.
        /// * `nnz` (`usize`): The capacity for non-zero elements.
        ///
        /// ## Returns
        /// `matrix.symmetric.Sparse(N, uplo, layout)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        /// * `matrix.Error.DimensionMismatch`: If `nnz` is greater than
        ///   `size * size`.
        pub fn init(allocator: std.mem.Allocator, size: usize, nnz: usize) !matrix.symmetric.Sparse(N, uplo, layout) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            if (nnz > size * size)
                return matrix.Error.DimensionMismatch;

            const data: []N = try allocator.alloc(N, nnz);
            errdefer allocator.free(data);

            const idx: []usize = try allocator.alloc(usize, nnz);
            errdefer allocator.free(idx);

            return .{
                .data = data.ptr,
                .idx = idx.ptr,
                .ptr = (try allocator.alloc(usize, size + 1)).ptr,
                .nnz = 0,
                .rows = size,
                .cols = size,
                ._dlen = nnz,
                ._ilen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Sparse(N, uplo, layout)`): A pointer to
        ///   the matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to compile `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.symmetric.Sparse(N, uplo, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.idx[0..self._ilen]);
                allocator.free(self.ptr[0 .. self.rows + 1]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.symmetric.Sparse(N, uplo, layout)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        pub fn get(self: matrix.symmetric.Sparse(N, uplo, layout), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: usize = r;
            var j: usize = c;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            }

            const major = if (comptime layout == .col_major) j else i;
            const minor = if (comptime layout == .col_major) i else j;

            var left: usize = self.ptr[major];
            var right: usize = self.ptr[major + 1];
            while (left < right) {
                const mid = left + (right - left) / 2;
                const mid_idx = self.idx[mid];

                if (mid_idx == minor) {
                    return self.data[mid];
                } else if (mid_idx < minor) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }

            return numeric.zero(N);
        }

        /// Sets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.symmetric.Sparse(N, uplo, layout)`): A pointer to
        ///   the matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set.
        /// * `c` (`usize`): The column index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        /// * `matrix.Error.BreaksStructure`: If no existing element is present
        ///   at the specified index.
        pub fn set(self: *matrix.symmetric.Sparse(N, uplo, layout), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            var i: usize = r;
            var j: usize = c;
            if (comptime uplo == .upper) {
                if (i > j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            } else {
                if (i < j) {
                    const temp: usize = i;
                    i = j;
                    j = temp;
                }
            }

            const major = if (comptime layout == .col_major) j else i;
            const minor = if (comptime layout == .col_major) i else j;

            var left: usize = self.ptr[major];
            var right: usize = self.ptr[major + 1];
            while (left < right) {
                const mid = left + (right - left) / 2;
                const mid_idx = self.idx[mid];

                if (mid_idx == minor) {
                    self.data[mid] = value;
                    return;
                } else if (mid_idx < minor) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }

            return matrix.Error.BreaksStructure;
        }
    };
}
