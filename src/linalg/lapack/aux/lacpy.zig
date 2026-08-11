const int = @import("../../../int.zig");
const matrix = @import("../../../matrix.zig");
const numeric = @import("../../../numeric.zig");

pub fn lacpy(
    layout: matrix.Layout,
    uplo: union(enum) { uplo: matrix.Uplo, full: void },
    m: usize,
    n: usize,
    a: anytype,
    lda: usize,
    b: anytype,
    ldb: usize,
) !void {
    switch (uplo) {
        .uplo => |_uplo| {
            if (layout == .col_major) {
                if (_uplo == .upper) {
                    for (0..n) |j| {
                        for (0..int.min(j + 1, m)) |i| {
                            numeric.set(&b[i + j * ldb], a[i + j * lda]);
                        }
                    }
                } else {
                    for (0..n) |j| {
                        for (j..m) |i| {
                            numeric.set(&b[i + j * ldb], a[i + j * lda]);
                        }
                    }
                }
            } else {
                if (_uplo == .upper) {
                    for (0..m) |i| {
                        for (i..n) |j| {
                            numeric.set(&b[i * ldb + j], a[i * lda + j]);
                        }
                    }
                } else {
                    for (0..m) |i| {
                        for (0..int.min(i + 1, n)) |j| {
                            numeric.set(&b[i * ldb + j], a[i * lda + j]);
                        }
                    }
                }
            }
        },
        .full => {
            if (layout == .col_major) {
                for (0..n) |j| {
                    for (0..m) |i| {
                        numeric.set(&b[i + j * ldb], a[i + j * lda]);
                    }
                }
            } else {
                for (0..m) |i| {
                    for (0..n) |j| {
                        numeric.set(&b[i * ldb + j], a[i * lda + j]);
                    }
                }
            }
        },
    }
}
