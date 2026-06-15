const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

/// Static permutation matrix type, represented as a contiguous array of `size`
/// elements of type `usize` holding a permutation of `0 .. size`. If
/// `direction` is forward, the element at index `i` indicates the column
/// index of the 1 in row `i`, i.e., if `data[i] = j`, then the element at
/// row `i` and column `j` is 1, and all other elements in row `i` are 0. If
/// `direction` is backward, the same applies but for columns, i.e., if
/// `data[j] = i`, then the element at row `i` and column `j` is 1, and all
/// other elements in column `j` are 0.
pub fn Static(size_: comptime_int, N: type, direction_: matrix.permutation.Direction) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.permutation.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [size]usize,

        // Type signatures
        pub const is_matrix = true;
        pub const is_static = true;
        pub const is_permutation = true;
        pub const storage_layout: ?matrix.Layout = null;
        pub const storage_uplo: ?matrix.Uplo = null;
        pub const storage_diag: ?matrix.Diag = null;
        pub const direction = direction_;

        pub const size = size_;
        pub const rows = size_;
        pub const cols = size_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: matrix.permutation.Static(size, N, direction) = .{
            .data = undefined,
        };

        pub const init = undefined;

        /// Gets the element at the specified index.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Static(size, N, direction)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get.
        /// * `c` (`usize`): The column index of the element to get.
        ///
        /// ## Returns
        /// `N`: The element at the specified index.
        ///
        /// ## Errors
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of
        ///   bounds.
        pub fn get(self: matrix.permutation.Static(size, N, direction), r: usize, c: usize) !N {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (comptime direction == .forward) {
                if (self.data[r] == c) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            } else {
                if (self.data[c] == r) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            }
        }

        /// Gets the element at the specified index without bounds checking.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Static(size, N, direction)`): The
        ///   matrix to get the element from.
        /// * `r` (`usize`): The row index of the element to get. Assumed to be
        ///   within bounds.
        /// * `c` (`usize`): The column index of the element to get. Assumed to
        ///   be within bounds.
        ///
        /// ## Returns
        /// `N`: The element at the specified position.
        pub fn getAssumeInBounds(self: matrix.permutation.Static(size, N, direction), r: usize, c: usize) N {
            if (comptime direction == .forward) {
                if (self.data[r] == c) {
                    return numeric.one(N);
                } else {
                    return numeric.zero(N);
                }
            } else {
                if (self.data[c] == r) {
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

        /// Returns a transposed copy of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Static(size, N, direction)`): The
        ///   matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.permutation.Static(size, N, direction.invert())`: The
        /// transposed matrix.
        pub fn transposeCopy(self: matrix.permutation.Static(size, N, direction)) matrix.permutation.Static(size, N, direction.invert()) {
            return .{ .data = self.data };
        }

        /// Returns a transposed view of the matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.permutation.Static(size, N, direction)`): The
        ///   matrix to transpose.
        ///
        /// ## Returns
        /// `matrix.permutation.Sparse(N, direction.invert())`: The transposed
        /// matrix.
        pub fn transposeView(self: matrix.permutation.Static(size, N, direction)) matrix.permutation.Sparse(N, direction.invert()) {
            return .{
                .data = &self.data,
                .rows = cols,
                .cols = rows,
                .flags = .{ .owns_data = false },
            };
        }
    };
}
