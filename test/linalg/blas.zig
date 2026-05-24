test {
    const test_level_1 = true;
    const test_level_2 = true;
    const test_level_3 = true;

    if (test_level_1) {
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

    if (test_level_2) {
        // _ = @import("blas/level2/gbmv.zig");
        // _ = @import("blas/level2/gemv.zig");
        // _ = @import("blas/level2/ger.zig");
        // _ = @import("blas/level2/gerc.zig");
        // _ = @import("blas/level2/hbmv.zig");
        // _ = @import("blas/level2/hemv.zig");
        // _ = @import("blas/level2/her.zig");
        // _ = @import("blas/level2/her2.zig");
        // _ = @import("blas/level2/hpmv.zig");
        // _ = @import("blas/level2/hpr.zig");
        // _ = @import("blas/level2/hpr2.zig");
        // _ = @import("blas/level2/sbmv.zig");
        // _ = @import("blas/level2/spmv.zig");
        // _ = @import("blas/level2/spr.zig");
        // _ = @import("blas/level2/spr2.zig");
        // _ = @import("blas/level2/symv.zig");
        // _ = @import("blas/level2/syr.zig");
        // _ = @import("blas/level2/syr2.zig");
        // _ = @import("blas/level2/tbmv.zig");
        // _ = @import("blas/level2/tbsv.zig");
        // _ = @import("blas/level2/tpmv.zig");
        // _ = @import("blas/level2/tpsv.zig");
        // _ = @import("blas/level2/trmv.zig");
        // _ = @import("blas/level2/trsv.zig");
    }

    if (test_level_3) {
        _ = @import("blas/level3/gemm.zig");
        _ = @import("blas/level3/gemmtr.zig");
        _ = @import("blas/level3/hemm.zig");
        _ = @import("blas/level3/herk.zig");
        _ = @import("blas/level3/her2k.zig");
        _ = @import("blas/level3/symm.zig");
        _ = @import("blas/level3/syrk.zig");
        _ = @import("blas/level3/syr2k.zig");
        _ = @import("blas/level3/trmm.zig");
        _ = @import("blas/level3/trsm.zig");
    }
}
