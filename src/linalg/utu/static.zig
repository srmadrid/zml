const std = @import("std");

const int = @import("../../int.zig");
const linalg = @import("../../linalg.zig");
const matrix = @import("../../matrix.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn Static(size_: comptime_int, N: type, layout: matrix.Layout) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.linalg.utu.Static: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [rows * cols]N,

        // Type signatures
        pub const is_utu = true;
        pub const is_static = true;
        pub const storage_layout: matrix.Layout = layout;
        pub const size = size_;
        pub const rows = size_;
        pub const cols = size_;

        // Numeric type
        pub const Numeric = N;

        pub const empty: linalg.utu.Static(size, N, layout) = .{
            .data = undefined,
        };

        pub const init = empty;

        pub fn u(self: linalg.utu.Static(size, N, layout)) matrix.triangular.Static(size, size, N, .upper, .non_unit, layout) {
            var _u: matrix.triangular.Static(size, size, N, .upper, .non_unit, layout) = .init;

            if (comptime layout == .col_major) {
                inline for (0..size) |j| {
                    inline for (0..j + 1) |i| {
                        numeric.set(&_u.data[i + j * size], self.data[i + j * size]);
                    }
                }
            } else {
                inline for (0..size) |i| {
                    inline for (i..size) |j| {
                        numeric.set(&_u.data[i * size + j], self.data[i * size + j]);
                    }
                }
            }

            return _u;
        }
    };
}
