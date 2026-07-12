const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

const matutils = @import("../utils.zig");

/// Sparse permutation matrix type, represented as a contiguous array of `size`
/// elements of type `usize` holding a permutation of `0 .. size`. If
/// `direction` is forward, the element at index `i` indicates the column
/// index of the 1 in row `i`, i.e., if `data[i] = j`, then the element at
/// row `i` and column `j` is 1, and all other elements in row `i` are 0. If
/// `direction` is backward, the same applies but for columns, i.e., if
/// `data[j] = i`, then the element at row `i` and column `j` is 1, and all
/// other elements in column `j` are 0.
pub fn Sparse(N: type, direction_: matrix.permutation.Direction) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.permutation.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        idx: [*]usize,
        rows: usize,
        cols: usize,
        flags: matrix.Flags,

        // Type signatures
        pub const is_matrix = true;
        pub const is_sparse = true;
        pub const is_permutation = true;
        pub const storage_layout: ?matrix.Layout = null;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;
        pub const direction = direction_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.permutation.Sparse(N, direction) = .{
            .idx = &.{},
            .rows = 0,
            .cols = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.permutation.Sparse(N, direction)` with the
        /// specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.permutation.Sparse(N, direction)`: The newly initialized
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn init(allocator: std.mem.Allocator, size: usize) !matrix.permutation.Sparse(N, direction) {
            if (size == 0)
                return matrix.Error.ZeroDimension;

            return .{
                .idx = (try allocator.alloc(usize, size)).ptr,
                .rows = size,
                .cols = size,
                .flags = .{ .owns_data = true },
            };
        }

        /// Initializes a new identity `matrix.permutation.Sparse(N, direction)`
        /// of the specified size.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `size` (`usize`): The size of the (square) matrix.
        ///
        /// ## Returns
        /// `matrix.permutation.Sparse(N, direction)`: The newly initialized
        /// identity matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If `size` is zero.
        pub fn initIdentity(allocator: std.mem.Allocator, size: usize) !matrix.permutation.Sparse(N, direction) {
            const mat: matrix.permutation.Sparse(N, direction) = try .init(allocator, size);

            var i: usize = 0;
            while (i < size) : (i += 1) {
                mat.idx[i] = i;
            }

            return mat;
        }

        /// Deinitializes the matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.permutation.Sparse(N, direction)`): A pointer to
        ///   the matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///  `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.permutation.Sparse(N, direction), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.idx[0..self.rows]);
            }

            self.* = undefined;
        }

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Sparse(N, direction)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.permutation.Sparse(N, direction), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime direction == .forward) {
                if (self.idx[r] == c) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            } else {
                if (self.idx[c] == r) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            }
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Sparse(N, direction)`): The matrix to
        ///   get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified position.
        pub fn getAssumeInBounds(self: matrix.permutation.Sparse(N, direction), r: usize, c: usize) N {
            if (comptime direction == .forward) {
                if (self.idx[r] == c) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            } else {
                if (self.idx[c] == r) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            }
        }

        // pub fn set(self: *Permutation(T), row: usize, col: usize, value: usize) !void {
        //     if (row >= self.rows or col >= self.cols)
        //         return matrix.Error.PositionOutOfBounds;

        //     if (value != 0 and value != 1)
        //         return matrix.Error.BreaksStructure;
        // }

        // pub  fn put(self: *Permutation(T), row: usize, col: usize, value: usize) void {
        //     // Unchecked version of set. Assumes row and col are valid and
        //     // in banded range.
        //     if (value == 1) {
        //         self.data[row] = col;
        //     }
        // }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Sparse(N, direction)`): The matrix to
        ///   transpose.
        ///
        /// ## Returns
        /// `matrix.permutation.Sparse(N, direction.invert())`: The transposed
        /// matrix.
        pub fn transposeView(self: matrix.permutation.Sparse(N, direction)) matrix.permutation.Sparse(N, direction.invert()) {
            return .{
                .idx = self.idx,
                .rows = self.cols,
                .cols = self.rows,
                .flags = self.flags,
            };
        }

        pub fn format(self: matrix.permutation.Sparse(N), writer: *std.Io.Writer) !void {
            try self.formatter("{e}").format(writer);
        }

        pub fn formatter(self: *const matrix.permutation.Sparse(N), comptime num_fmt: []const u8) Formatter(num_fmt) {
            return .{ .matrix = self };
        }

        pub fn Formatter(comptime num_fmt: []const u8) type {
            return struct {
                matrix: *const matrix.permutation.Sparse(N),

                pub fn format(self: matrix.permutation.Sparse(N).Formatter(num_fmt), writer: *std.Io.Writer) !void {
                    const rows = self.matrix.rows;
                    const cols = self.matrix.cols;

                    try writer.print("zsl.matrix.permutation.Sparse({s}) ({d} × {d}):\n\n", .{ @typeName(N), rows, cols });

                    return matutils.format(self, num_fmt, rows, cols, writer);
                }
            };
        }
    };
}
