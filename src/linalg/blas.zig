const meta = @import("../meta.zig");

pub const Transpose = enum(u2) {
    no_trans,
    trans,
    conj_trans,
    conj_no_trans,

    pub fn toInt(self: Transpose, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.linalg.Transpose.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .no_trans => if (comptime Int == u8) 'N' else 111,
            .trans => if (comptime Int == u8) 'T' else 112,
            .conj_trans => if (comptime Int == u8) 'C' else 113,
            .conj_no_trans => if (comptime Int == u8) 0 else 114,
        };
    }

    pub fn invert(self: Transpose) Transpose {
        return switch (self) {
            .no_trans => .trans,
            .trans => .no_trans,
            .conj_no_trans => .conj_trans,
            .conj_trans => .conj_no_trans,
        };
    }

    pub fn reverse(self: Transpose) Transpose {
        return switch (self) {
            .no_trans => .conj_trans,
            .trans => .conj_no_trans,
            .conj_no_trans => .trans,
            .conj_trans => .no_trans,
        };
    }
};

pub const Side = enum(u1) {
    left,
    right,

    pub fn toInt(self: Side, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.linalg.Side.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .left => if (comptime Int == u8) 'L' else 141,
            .right => if (comptime Int == u8) 'R' else 142,
        };
    }

    pub fn invert(self: Side) Side {
        return switch (self) {
            .left => .right,
            .right => .left,
        };
    }
};

// Level 1
pub const Asum = @import("blas/level1/asum.zig").Asum;
pub const asum = @import("blas/level1/asum.zig").asum;
pub const asumParallel = @import("blas/level1/asum.zig").asumParallel;
pub const axpy = @import("blas/level1/axpy.zig").axpy;
pub const axpyParallel = @import("blas/level1/axpy.zig").axpyParallel;
pub const copy = @import("blas/level1/copy.zig").copy;
pub const copyParallel = @import("blas/level1/copy.zig").copyParallel;
pub const Dot = @import("blas/level1/dot.zig").Dot;
pub const dot = @import("blas/level1/dot.zig").dot;
pub const dotParallel = @import("blas/level1/dot.zig").dotParallel;
pub const Dotc = @import("blas/level1/dotc.zig").Dotc;
pub const dotc = @import("blas/level1/dotc.zig").dotc;
pub const dotcParallel = @import("blas/level1/dotc.zig").dotcParallel;
pub const Nrm2 = @import("blas/level1/nrm2.zig").Nrm2;
pub const nrm2 = @import("blas/level1/nrm2.zig").nrm2;
pub const nrm2Parallel = @import("blas/level1/nrm2.zig").nrm2Parallel;
pub const rot = @import("blas/level1/rot.zig").rot;
pub const rotParallel = @import("blas/level1/rot.zig").rotParallel;
pub const scal = @import("blas/level1/scal.zig").scal;
pub const scalParallel = @import("blas/level1/scal.zig").scalParallel;
pub const swap = @import("blas/level1/swap.zig").swap;
pub const swapParallel = @import("blas/level1/swap.zig").swapParallel;
pub const iamax = @import("blas/level1/iamax.zig").iamax;
pub const iamaxParallel = @import("blas/level1/iamax.zig").iamaxParallel;
pub const iamin = @import("blas/level1/iamin.zig").iamin;
pub const iaminParallel = @import("blas/level1/iamin.zig").iaminParallel;

// Level 2
pub const gemv = @import("blas/level2/gemv.zig").gemv;
pub const gemvParallel = @import("blas/level2/gemv.zig").gemvParallel;
pub const ger = @import("blas/level2/ger.zig").ger;
pub const gerParallel = @import("blas/level2/ger.zig").gerParallel;
pub const gerc = @import("blas/level2/gerc.zig").gerc;
pub const gercParallel = @import("blas/level2/gerc.zig").gercParallel;
pub const hemv = @import("blas/level2/hemv.zig").hemv; // Unparallelized
pub const her = @import("blas/level2/her.zig").her; // Unparallelized
pub const her2 = @import("blas/level2/her2.zig").her2; // Unparallelized
pub const symv = @import("blas/level2/symv.zig").symv; // Unparallelized
pub const syr = @import("blas/level2/syr.zig").syr; // Unparallelized
pub const syr2 = @import("blas/level2/syr2.zig").syr2; // Unparallelized
pub const trmv = @import("blas/level2/trmv.zig").trmv;
pub const trmvParallel = @import("blas/level2/trmv.zig").trmvParallel;
pub const trsv = @import("blas/level2/trsv.zig").trsv; // Unparallelized

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
