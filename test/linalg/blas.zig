test {
    const test_level_1 = true;
    const test_level_2 = false;
    const test_level_3 = false;

    if (test_level_1) {
        _ = @import("blas/asum.zig");
        _ = @import("blas/axpy.zig");
        _ = @import("blas/copy.zig");
        _ = @import("blas/dot.zig");
        _ = @import("blas/dotc.zig");
        _ = @import("blas/dotc_sub.zig");
        _ = @import("blas/dotu.zig");
        _ = @import("blas/dotu_sub.zig");
        _ = @import("blas/nrm2.zig");
        _ = @import("blas/rot.zig");
        _ = @import("blas/rotg.zig");
        _ = @import("blas/rotm.zig");
        _ = @import("blas/rotmg.zig");
        _ = @import("blas/scal.zig");
        _ = @import("blas/swap.zig");
        _ = @import("blas/iamax.zig");
        _ = @import("blas/iamin.zig");
    }

    if (test_level_2) {
        _ = @import("blas/gbmv.zig");
        _ = @import("blas/gemv.zig");
        _ = @import("blas/ger.zig");
        _ = @import("blas/gerc.zig");
        _ = @import("blas/geru.zig");
        _ = @import("blas/hbmv.zig");
        _ = @import("blas/hemv.zig");
        _ = @import("blas/her.zig");
        _ = @import("blas/her2.zig");
        _ = @import("blas/hpmv.zig");
        _ = @import("blas/hpr.zig");
        _ = @import("blas/hpr2.zig");
        _ = @import("blas/sbmv.zig");
        _ = @import("blas/spmv.zig");
        _ = @import("blas/spr.zig");
        _ = @import("blas/spr2.zig");
        _ = @import("blas/symv.zig");
        _ = @import("blas/syr.zig");
        _ = @import("blas/syr2.zig");
        _ = @import("blas/tbmv.zig");
        _ = @import("blas/tbsv.zig");
        _ = @import("blas/tpmv.zig");
        _ = @import("blas/tpsv.zig");
        _ = @import("blas/trmv.zig");
        _ = @import("blas/trsv.zig");
    }

    if (test_level_3) {
        _ = @import("blas/gemm.zig");
        _ = @import("blas/hemm.zig");
        _ = @import("blas/herk.zig");
        _ = @import("blas/her2k.zig");
        _ = @import("blas/symm.zig");
        _ = @import("blas/syrk.zig");
        _ = @import("blas/syr2k.zig");
        _ = @import("blas/trmm.zig");
        _ = @import("blas/trsm.zig");
    }
}
