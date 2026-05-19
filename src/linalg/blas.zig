// Level 1
pub const Asum = @import("blas/asum.zig").Asum;
pub const asum = @import("blas/asum.zig").asum;
pub const axpy = @import("blas/axpy.zig").axpy;
pub const copy = @import("blas/copy.zig").copy;
pub const Dot = @import("blas/dot.zig").Dot;
pub const dot = @import("blas/dot.zig").dot;
pub const Dotc = @import("blas/dotc.zig").Dotc;
pub const dotc = @import("blas/dotc.zig").dotc;
pub const Nrm2 = @import("blas/nrm2.zig").Nrm2;
pub const nrm2 = @import("blas/nrm2.zig").nrm2;
pub const rot = @import("blas/rot.zig").rot;
// pub const rotg = @import("blas/rotg.zig").rotg;
// pub const rotm = @import("blas/rotm.zig").rotm;
// pub const rotmg = @import("blas/rotmg.zig").rotmg;
pub const scal = @import("blas/scal.zig").scal;
pub const swap = @import("blas/swap.zig").swap;
pub const iamax = @import("blas/iamax.zig").iamax;
pub const iamin = @import("blas/iamin.zig").iamin;

// Level 2
pub const gbmv = @import("blas/gbmv.zig").gbmv;
pub const gemv = @import("blas/gemv.zig").gemv;
pub const ger = @import("blas/ger.zig").ger;
// pub const gerc = @import("blas/gerc.zig").gerc;
pub const hbmv = @import("blas/hbmv.zig").hbmv;
// pub const hemv = @import("blas/hemv.zig").hemv;
// pub const her = @import("blas/her.zig").her;
// pub const her2 = @import("blas/her2.zig").her2;
pub const hpmv = @import("blas/hpmv.zig").hpmv;
pub const hpr = @import("blas/hpr.zig").hpr;
pub const hpr2 = @import("blas/hpr2.zig").hpr2;
pub const sbmv = @import("blas/sbmv.zig").sbmv;
pub const spmv = @import("blas/spmv.zig").spmv;
pub const spr = @import("blas/spr.zig").spr;
pub const spr2 = @import("blas/spr2.zig").spr2;
// pub const symv = @import("blas/symv.zig").symv;
// pub const syr = @import("blas/syr.zig").syr;
// pub const syr2 = @import("blas/syr2.zig").syr2;
pub const tbmv = @import("blas/tbmv.zig").tbmv;
pub const tbsv = @import("blas/tbsv.zig").tbsv;
pub const tpmv = @import("blas/tpmv.zig").tpmv;
// pub const tpsv = @import("blas/tpsv.zig").tpsv;
// pub const trmv = @import("blas/trmv.zig").trmv;
// pub const trsv = @import("blas/trsv.zig").trsv;

// Level 3
// pub const gemm = @import("blas/gemm.zig").gemm;
// pub const gemmtr = @import("blas/gemmtr.zig").gemmtr;
// pub const hemm = @import("blas/hemm.zig").hemm;
// pub const herk = @import("blas/herk.zig").herk;
// pub const her2k = @import("blas/her2k.zig").her2k;
// pub const symm = @import("blas/symm.zig").symm;
// pub const syrk = @import("blas/syrk.zig").syrk;
// pub const syr2k = @import("blas/syr2k.zig").syr2k;
// pub const trmm = @import("blas/trmm.zig").trmm;
// pub const trsm = @import("blas/trsm.zig").trsm;

pub const Error = error{
    InvalidArgument,
};
