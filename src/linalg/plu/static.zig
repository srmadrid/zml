const std = @import("std");

const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");
const matrix = @import("../../matrix.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn Static(rows_: comptime_int, cols_: comptime_int, N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.linalg.plu.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [rows * cols]N,
        pivots: [int.min(rows, cols)]usize,

        // Type signatures
        pub const is_plu = true;
        pub const is_static = true;
        pub const storage_layout: matrix.Layout = layout;
        pub const rows = rows_;
        pub const cols = cols_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: linalg.plu.Static(rows, cols, N, layout) = .{
            .data = undefined,
            .pivots = undefined,
        };

        pub const init = empty;

        pub fn p(self: linalg.plu.Static(rows, cols, N, layout)) matrix.permutation.Static(rows, N, .forward) {
            var _p: matrix.permutation.Static(rows, N, .forward) = .init;

            inline for (0..rows) |i| {
                _p.idx[i] = i;
            }

            comptime var i: usize = int.min(rows, cols);
            inline while (i > 0) {
                i -= 1;

                const temp = _p.idx[i];
                _p.idx[i] = _p.idx[self.pivots[i]];
                _p.idx[self.pivots[i]] = temp;
            }

            return _p;
        }

        pub fn l(self: linalg.plu.Static(rows, cols, N, layout)) matrix.triangular.Static(rows, int.min(rows, cols), N, .lower, .unit, layout) {
            var _l: matrix.triangular.Static(rows, int.min(rows, cols), N, .lower, .unit, layout) = .init;

            if (comptime layout == .col_major) {
                inline for (0..int.min(rows, cols)) |j| {
                    inline for (j + 1..rows) |i| {
                        numeric.set(&_l.data[i + j * rows], self.data[i + j * rows]);
                    }
                }
            } else {
                inline for (0..rows) |i| {
                    inline for (0..i) |j| {
                        numeric.set(&_l.data[i * int.min(rows, cols) + j], self.data[i * int.min(rows, cols) + j]);
                    }
                }
            }

            return _l;
        }

        pub fn u(self: linalg.plu.Static(rows, cols, N, layout)) matrix.triangular.Static(int.min(rows, cols), cols, N, .upper, .non_unit, layout) {
            var _u: matrix.triangular.Static(int.min(rows, cols), cols, N, .upper, .non_unit, layout) = .init;

            if (comptime layout == .col_major) {
                inline for (0..cols) |j| {
                    inline for (0..j + 1) |i| {
                        numeric.set(&_u.data[i + j * int.min(rows, cols)], self.data[i + j * int.min(rows, cols)]);
                    }
                }
            } else {
                inline for (0..int.min(rows, cols)) |i| {
                    inline for (i..cols) |j| {
                        numeric.set(&_u.data[i * cols + j], self.data[i * cols + j]);
                    }
                }
            }

            return _u;
        }
    };
}
