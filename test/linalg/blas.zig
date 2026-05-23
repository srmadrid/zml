test {
    const test_level_1 = true;
    const test_level_2 = true;
    const test_level_3 = true;

    if (test_level_1) {
        _ = @import("blas/asum.zig");
        _ = @import("blas/axpy.zig");
        _ = @import("blas/copy.zig");
        _ = @import("blas/dot.zig");
        _ = @import("blas/dotc.zig");
        _ = @import("blas/nrm2.zig");
        _ = @import("blas/scal.zig");
        _ = @import("blas/swap.zig");
    }

    if (test_level_2) {}

    if (test_level_3) {}
}
