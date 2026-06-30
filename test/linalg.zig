const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("zsl.zig");

pub fn correctMatmul(comptime O: type, allocator: std.mem.Allocator, rows: usize, K: usize, cols: usize, A: anytype, B: anytype) !if (zsl.meta.isVector(@TypeOf(A)) or zsl.meta.isVector(@TypeOf(B))) zsl.vector.Dense(O) else zsl.matrix.general.Dense(O, .col_major) {
    if (comptime zsl.meta.isVector(@TypeOf(A)) or zsl.meta.isVector(@TypeOf(B))) {
        if (comptime zsl.meta.isVector(@TypeOf(A))) {
            const result: zsl.vector.Dense(O) = try .init(allocator, cols);

            var j: usize = 0;
            while (j < cols) : (j += 1) {
                var sum = zsl.numeric.zero(zsl.meta.Accumulator(O));

                var k: usize = 0;
                while (k < K) : (k += 1) {
                    zsl.numeric.fmaInto(
                        &sum,
                        A.get(k) catch unreachable,
                        B.get(k, j) catch unreachable,
                        sum,
                    );
                }

                zsl.numeric.set(&result.data[result._index(j)], sum);
            }

            return result;
        } else {
            const result: zsl.vector.Dense(O) = try .init(allocator, rows);

            var i: usize = 0;
            while (i < rows) : (i += 1) {
                var sum = zsl.numeric.zero(zsl.meta.Accumulator(O));

                var k: usize = 0;
                while (k < K) : (k += 1) {
                    zsl.numeric.fmaInto(
                        &sum,
                        A.get(i, k) catch unreachable,
                        B.get(k) catch unreachable,
                        sum,
                    );
                }

                zsl.numeric.set(&result.data[result._index(i)], sum);
            }

            return result;
        }
    } else {
        const result: zsl.matrix.general.Dense(O, .col_major) = try .init(allocator, rows, cols);

        var j: usize = 0;
        while (j < result.cols) : (j += 1) {
            var i: usize = 0;
            while (i < result.rows) : (i += 1) {
                var sum = zsl.numeric.zero(zsl.meta.Accumulator(O));

                var k: usize = 0;
                while (k < K) : (k += 1) {
                    zsl.numeric.fmaInto(
                        &sum,
                        A.get(i, k) catch unreachable,
                        B.get(k, j) catch unreachable,
                        sum,
                    );
                }

                zsl.numeric.set(&result.data[result._index(i, j)], sum);
            }
        }

        return result;
    }
}

test {
    const test_blas = true;
    const test_matmul = false;

    if (test_blas) {
        _ = @import("linalg/blas.zig");
    }

    if (test_matmul) {
        _ = @import("linalg/matmul.zig");
        _ = @import("linalg/matmulAlloc.zig");
    }
}
