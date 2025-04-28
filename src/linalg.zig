const std = @import("std");

pub const blas = @import("linalg/blas.zig");
//pub const lapack = @import("linalg/lapack.zig");

test {
    const test_blas = true;

    if (test_blas) {
        _ = @import("linalg/blas.zig");
    }
}
