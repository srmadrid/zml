const std = @import("std");

const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");
const matrix = @import("../../matrix.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn Static(size_: comptime_int, N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.linalg.llt.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [rows * cols]N,

        // Type signatures
        pub const is_llt = true;
        pub const is_static = true;
        pub const storage_layout: matrix.Layout = layout;
        pub const size = size_;
        pub const rows = size_;
        pub const cols = size_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: linalg.llt.Static(size, N, layout) = .{
            .data = undefined,
        };

        pub const init = empty;

        pub fn l(self: linalg.llt.Static(size, N, layout)) matrix.triangular.Static(size, size, N, .lower, .non_unit, layout) {
            var _l: matrix.triangular.Static(size, size, N, .lower, .non_unit, layout) = .init;

            if (comptime layout == .col_major) {
                inline for (0..size) |j| {
                    inline for (j..size) |i| {
                        numeric.set(&_l.data[i + j * size], self.data[i + j * size]);
                    }
                }
            } else {
                inline for (0..size) |i| {
                    inline for (0..i + 1) |j| {
                        numeric.set(&_l.data[i * size + j], self.data[i * size + j]);
                    }
                }
            }

            return _l;
        }
    };
}
