const meta = @import("meta.zig");

const int = @import("int.zig");

pub const cblas = @import("linalg/cblas.zig");
pub const blas = @import("linalg/blas.zig");
// pub const lapacke = @import("linalg/lapacke.zig");
// pub const lapack = @import("linalg/lapack.zig");

// Vector operations
pub const Dot = @import("linalg/dot.zig").Dot;
pub const dot = @import("linalg/dot.zig").dot;

// outer

pub const Cross = @import("linalg/cross.zig").Cross;
pub const cross = @import("linalg/cross.zig").cross;

pub const Normalize = @import("linalg/normalize.zig").Normalize;
pub const normalize = @import("linalg/normalize.zig").normalize;
pub const normalizeAlloc = @import("linalg/normalize.zig").normalizeAlloc;
pub const normalizeInto = @import("linalg/normalize.zig").normalizeInto;

// Vector/matrix operations
pub const Norm = @import("linalg/norm.zig").Norm;
pub const NormOrder = @import("linalg/norm.zig").NormOrder;
pub const norm = @import("linalg/norm.zig").norm;
pub const normAlloc = @import("linalg/norm.zig").normAlloc;
// pub const normInto = @import("linalg/norm.zig").normInto;

// pub const matmul = @import("linalg/matmul.zig").matmul;

// solve (with SolveMethod, to choose decomposition, and optional out parameter to save it)
// Allow vector and matrix rhs

// leastSquares
// Allow vector and matrix rhs

// Matrix operations
// trace

// det

// inv

// pinv

// const _lu = @import("linalg/lu.zig");
// pub const LU = _lu.LU;
// pub const lu = _lu.lu;
// pub const PLU = _lu.PLU;
// pub const plu = _lu.plu;
// pub const PLUQ = _lu.PLUQ;
// pub const pluq = _lu.pluq;

// const _cholesky = @import("linalg/cholesky.zig");
// pub const LLT = _cholesky.LLT;
// pub const llt = _cholesky.llt;
// pub const UTU = _cholesky.UTU;
// pub const utu = _cholesky.utu;
// pub const cholesky = _cholesky.cholesky;

// const _bunchkaufman = @import("linalg/bunchkaufman.zig");
// pub const LDLT = _bunchkaufman.LDLT;
// pub const ldlt = _bunchkaufman.ldlt;
// pub const UDUT = _bunchkaufman.UDUT;
// pub const udut = _bunchkaufman.udut;
// pub const bunchkaufman = _bunchkaufman.bunchkaufman;

// const qr_ = @import("linalg/qr.zig");
// pub const QR = qr_.QR;
// pub const qr = qr_.qr;
// pub const QRP = qr_.QRP;
// pub const qrp = qr_.qrp;

// svd

// eig

// schur

// Array (tensor) operations
// contract (generalization of the dot product, C[i, j, k] = sum(l, A[i, j, l] * B[l, k])

// einsum

// tensorProduct (generalization of the outer product, A (shape i, j) ⊗ B (shape k, l) gives (shape i, j, k, l))

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

pub const Error = error{
    DimensionMismatch,
    FactorizationFailed,
    SingularMatrix,
    IndefiniteMatrix,
    NotImplemented,
};
