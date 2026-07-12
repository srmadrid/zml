const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");

const matrix = @import("../../matrix.zig");

const matutils = @import("../utils.zig");

/// Sparse general matrix type, represented in either CSC or CSR format,
/// depending on if `layout` is column-major or row-major, respectively.
pub fn Sparse(N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.general.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

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
        pub const is_general = true;
        pub const storage_layout: ?matrix.Layout = layout;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.general.Sparse(N, layout) = .{
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

        /// Initializes a new `matrix.general.Sparse(N, layout)` with the
        /// specified rows and columns, and a capacity for non-zero elements.
        ///
        /// Initialization with this function is only meant for pre-allocating
        /// an empty matrix to serve as the output destination for in-place
        /// mathematical functions.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the matrix.
        /// * `cols` (`usize`): The columns of the matrix.
        /// * `nnz` (`usize`): The capacity for non-zero elements.
        ///
        /// ## Returns
        /// `matrix.general.Sparse(N, layout)`: The newly initialized matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        /// * `matrix.Error.DimensionMismatch`: If `nnz` is greater than
        ///   `rows * cols`.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize, nnz: usize) !matrix.general.Sparse(N, layout) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (nnz > rows * cols)
                return matrix.Error.DimensionMismatch;

            const data: []N = try allocator.alloc(N, nnz);
            errdefer allocator.free(data);

            const idx: []usize = try allocator.alloc(usize, nnz);
            errdefer allocator.free(idx);

            return .{
                .data = data.ptr,
                .idx = idx.ptr,
                .ptr = (try allocator.alloc(usize, if (comptime layout == .col_major) cols + 1 else rows + 1)).ptr,
                .nnz = 0,
                .rows = rows,
                .cols = cols,
                ._dlen = nnz,
                ._ilen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.general.Sparse(N, layout)`): A pointer to the
        ///   matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to compile `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.general.Sparse(N, layout), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.idx[0..self._ilen]);
                allocator.free(self.ptr[0..(if (comptime layout == .col_major) self.cols + 1 else self.rows + 1)]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.general.Sparse(N, layout)`): The matrix to get the
        ///   element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.general.Sparse(N, layout), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            const major = if (comptime layout == .col_major) c else r;
            const minor = if (comptime layout == .col_major) r else c;

            var left: usize = self.ptr[major];
            var right: usize = self.ptr[major + 1];
            while (left < right) {
                const mid = left + (right - left) / 2;
                const mid_idx = self.idx[mid];

                if (mid_idx == minor) {
                    return if (self.flags.noconj) self.data[mid] else numeric.conj(self.data[mid]);
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
        /// * `self` (`*matrix.general.Sparse(N, layout)`): A pointer to the
        ///   matrix to set the element in.
        /// * `r` (`usize`): The row index of the element to set.
        /// * `c` (`usize`): The column index of the element to set.
        /// * `value` (`N`): The value to set the element to.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        /// * `matrix.Error.BreaksStructure`: If no existing element is present
        ///   at the specified index.
        pub fn set(self: *matrix.general.Sparse(N, layout), r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            const major = if (comptime layout == .col_major) c else r;
            const minor = if (comptime layout == .col_major) r else c;

            var left: usize = self.ptr[major];
            var right: usize = self.ptr[major + 1];
            while (left < right) {
                const mid = left + (right - left) / 2;
                const mid_idx = self.idx[mid];

                if (mid_idx == minor) {
                    self.data[mid] = if (self.flags.noconj) value else numeric.conj(value);
                    return;
                } else if (mid_idx < minor) {
                    left = mid + 1;
                } else {
                    right = mid;
                }
            }

            return matrix.Error.BreaksStructure;
        }

        pub fn format(self: matrix.general.Sparse(N, layout), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const matrix.general.Sparse(N, layout), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .matrix = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                matrix: *const matrix.general.Sparse(N, layout),

                pub fn format(self: matrix.general.Sparse(N, layout).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const rows = self.matrix.rows;
                    const cols = self.matrix.cols;

                    try writer.print("zsl.matrix.general.Sparse({s}, .{s}) ({d} × {d}):\n\n", .{ @typeName(N), @tagName(layout), rows, cols });

                    return matutils.format(self, num_fmt, rows, cols, writer);
                }
            };
        }
    };
}
