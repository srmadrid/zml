test {
    // Override test flags
    const test_all = false;

    // Individual test flags
    const test_level_1 = false;
    const test_level_2 = true;
    const test_level_3 = false;

    if (test_all or test_level_1) {
        _ = @import("blas/level1/asum.zig");
        _ = @import("blas/level1/axpy.zig");
        _ = @import("blas/level1/copy.zig");
        _ = @import("blas/level1/dot.zig");
        _ = @import("blas/level1/dotc.zig");
        _ = @import("blas/level1/nrm2.zig");
        // _ = @import("blas/level1/rot.zig");
        // _ = @import("blas/level1/rotg.zig");
        // _ = @import("blas/level1/rotm.zig");
        // _ = @import("blas/level1/rotmg.zig");
        _ = @import("blas/level1/scal.zig");
        _ = @import("blas/level1/swap.zig");
        _ = @import("blas/level1/iamax.zig");
        _ = @import("blas/level1/iamin.zig");
    }

    if (test_all or test_level_2) {
        _ = @import("blas/level2/gemv.zig");
        _ = @import("blas/level2/ger.zig");
        _ = @import("blas/level2/gerc.zig");
        _ = @import("blas/level2/hemv.zig");
        _ = @import("blas/level2/her.zig");
        _ = @import("blas/level2/her2.zig");
        _ = @import("blas/level2/symv.zig");
        _ = @import("blas/level2/syr.zig");
        _ = @import("blas/level2/syr2.zig");
        // _ = @import("blas/level2/trmv.zig");
        // _ = @import("blas/level2/trsv.zig");
    }

    if (test_all or test_level_3) {
        // _ = @import("blas/level3/gemm.zig");
        // _ = @import("blas/level3/gemmtr.zig");
        // _ = @import("blas/level3/hemm.zig");
        // _ = @import("blas/level3/herk.zig");
        // _ = @import("blas/level3/her2k.zig");
        // _ = @import("blas/level3/symm.zig");
        // _ = @import("blas/level3/syrk.zig");
        // _ = @import("blas/level3/syr2k.zig");
        // _ = @import("blas/level3/trmm.zig");
        // _ = @import("blas/level3/trsm.zig");
    }
}
