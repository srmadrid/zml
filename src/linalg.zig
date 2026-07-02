const meta = @import("meta.zig");

const int = @import("int.zig");

pub const cblas = @import("linalg/cblas.zig");
pub const blas = @import("linalg/blas.zig");
pub const lapacke = @import("linalg/lapacke.zig");
// pub const lapack = @import("linalg/lapack.zig");

// Vector operations
pub const Dot = @import("linalg/dot.zig").Dot;
pub const dot = @import("linalg/dot.zig").dot;
pub const dotUnchecked = @import("linalg/dot.zig").dotUnchecked;

// outer

pub const Cross = @import("linalg/cross.zig").Cross;
pub const cross = @import("linalg/cross.zig").cross;

pub const Normalize = @import("linalg/normalize.zig").Normalize;
pub const normalize = @import("linalg/normalize.zig").normalize;
pub const normalizeAlloc = @import("linalg/normalize.zig").normalizeAlloc;
pub const normalizeInto = @import("linalg/normalize.zig").normalizeInto;

// distance (with norm order)

// project(x, y): ((x ⋅ y)/(y ⋅ y)) * y

// reject(x, y): x - project(x, y)

// Vector/matrix operations
pub const Norm = @import("linalg/norm.zig").Norm;
pub const NormOrder = @import("linalg/norm.zig").NormOrder;
pub const norm = @import("linalg/norm.zig").norm;
pub const normAlloc = @import("linalg/norm.zig").normAlloc;

pub const Matmul = @import("linalg/matmul.zig").Matmul;
pub const matmul = @import("linalg/matmul.zig").matmul;
pub const matmulUnchecked = @import("linalg/matmul.zig").matmulUnchecked;
pub const matmulAlloc = @import("linalg/matmul.zig").matmulAlloc;
pub const matmulInto = @import("linalg/matmul.zig").matmulInto;
pub const matmulIntoAlloc = @import("linalg/matmul.zig").matmulIntoAlloc;
pub const matmulIntoUnchecked = @import("linalg/matmul.zig").matmulIntoUnchecked;

// hadamard (also for vectors)

// solve (with SolveMethod, to choose decomposition, and optional out parameter to save it)
// Allow vector and matrix rhs

// lstsq
// Allow vector and matrix rhs

// Matrix operations
// trace

// det

// cond (condition number)

// gram (gramian matrix, Xᴴ X)

// cong (matrix congruence, A X Aᴴ)

// kron (A ⊗ B; (m × n) ⊗ (p × q) -> (mp × nq))

// expm

// sqrtm

// inv

// pinv

// pub const LU = @import("linalg/lu.zig").LU;
// pub const lu = @import("linalg/lu.zig").lu;
// pub const PLU = @import("linalg/lu.zig").PLU;
// pub const plu = @import("linalg/lu.zig").plu;
// pub const PLUQ = @import("linalg/lu.zig").PLUQ;
// pub const pluq = @import("linalg/lu.zig").pluq;

// ilu

// pub const LLT = @import("linalg/cholesky.zig").LLT;
// pub const llt = @import("linalg/cholesky.zig").llt;
// pub const UTU = @import("linalg/cholesky.zig").UTU;
// pub const utu = @import("linalg/cholesky.zig").utu;
// pub const cholesky = @import("linalg/cholesky.zig").cholesky;

// pub const LDLT = @import("linalg/bunchkaufman.zig").LDLT;
// pub const ldlt = @import("linalg/bunchkaufman.zig").ldlt;
// pub const UDUT = @import("linalg/bunchkaufman.zig").UDUT;
// pub const udut = @import("linalg/bunchkaufman.zig").udut;
// pub const bunchkaufman = @import("linalg/bunchkaufman.zig").bunchkaufman;

// pub const QR = @import("linalg/qr.zig").QR;
// pub const qr = @import("linalg/qr.zig").qr;
// pub const QRP = @import("linalg/qr.zig").QRP;
// pub const qrp = @import("linalg/qr.zig").qrp;

// rq

// ql

// lq

// rz

// svd

// eig

// schur

// Array (tensor) operations
// contract (generalization of the dot product, C[i, j, k] = sum(l, A[i, j, l] * B[l, k])

// einsum

// tensor (tensor product, generalization of the outer product, A (shape i, j) ⊗ B (shape k, l) gives (shape i, j, k, l))

pub const Error = error{
    DimensionMismatch,
    InsufficientSpace,
    FactorizationFailed,
    SingularMatrix,
    IndefiniteMatrix,
    NotImplemented,
};
