// Level 1
pub const Asum = @import("blas/level1/asum.zig").Asum;
pub const asum = @import("blas/level1/asum.zig").asum;
pub const axpy = @import("blas/level1/axpy.zig").axpy;
pub const copy = @import("blas/level1/copy.zig").copy;
pub const Dot = @import("blas/level1/dot.zig").Dot;
pub const dot = @import("blas/level1/dot.zig").dot;
pub const Dotc = @import("blas/level1/dotc.zig").Dotc;
pub const dotc = @import("blas/level1/dotc.zig").dotc;
pub const Nrm2 = @import("blas/level1/nrm2.zig").Nrm2;
pub const nrm2 = @import("blas/level1/nrm2.zig").nrm2;
pub const rot = @import("blas/level1/rot.zig").rot; // Untested and unbenchmarked
// pub const rotg = @import("blas/level1/rotg.zig").rotg;
// pub const rotm = @import("blas/level1/rotm.zig").rotm;
// pub const rotmg = @import("blas/level1/rotmg.zig").rotmg;
pub const scal = @import("blas/level1/scal.zig").scal;
pub const swap = @import("blas/level1/swap.zig").swap;
pub const iamax = @import("blas/level1/iamax.zig").iamax;
pub const iamin = @import("blas/level1/iamin.zig").iamin;

// Level 2
pub const gemv = @import("blas/level2/gemv.zig").gemv;
pub const ger = @import("blas/level2/ger.zig").ger;
pub const gerc = @import("blas/level2/gerc.zig").gerc;
pub const hemv = @import("blas/level2/hemv.zig").hemv;
// pub const her = @import("blas/level2/her.zig").her;
// pub const her2 = @import("blas/level2/her2.zig").her2;
pub const symv = @import("blas/level2/symv.zig").symv;
// pub const syr = @import("blas/level2/syr.zig").syr;
// pub const syr2 = @import("blas/level2/syr2.zig").syr2;
// pub const tpsv = @import("blas/level2/tpsv.zig").tpsv;
// pub const trmv = @import("blas/level2/trmv.zig").trmv;
// pub const trsv = @import("blas/level2/trsv.zig").trsv;

// Level 3
// pub const gemm = @import("blas/level3/gemm.zig").gemm;
// pub const gemmtr = @import("blas/level3/gemmtr.zig").gemmtr;
// pub const hemm = @import("blas/level3/hemm.zig").hemm;
// pub const herk = @import("blas/level3/herk.zig").herk;
// pub const her2k = @import("blas/level3/her2k.zig").her2k;
// pub const symm = @import("blas/level3/symm.zig").symm;
// pub const syrk = @import("blas/level3/syrk.zig").syrk;
// pub const syr2k = @import("blas/level3/syr2k.zig").syr2k;
// pub const trmm = @import("blas/level3/trmm.zig").trmm;
// pub const trsm = @import("blas/level3/trsm.zig").trsm;

pub const Error = error{
    InvalidArgument,
};
