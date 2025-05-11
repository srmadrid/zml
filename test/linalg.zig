test {
    const test_blas = true;

    if (test_blas) {
        _ = @import("linalg/blas.zig");
    }
}
