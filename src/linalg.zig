pub const cblas = @import("linalg/cblas.zig");
pub const blas = @import("linalg/blas.zig");
pub const lapacke = @import("linalg/lapacke.zig");
pub const lapack = @import("linalg/lapack.zig");

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
pub const normalizeIntoUnchecked = @import("linalg/normalize.zig").normalizeIntoUnchecked;

// distance (with norm order) without realizing y - x, i.e., computing it on the go

// project(x, y): ((x ⋅ y)/(y ⋅ y)) * y

// reject(x, y): x - project(x, y)

// Vector/matrix operations
pub const Norm = @import("linalg/norm.zig").Norm;
pub const NormOrder = @import("linalg/norm.zig").NormOrder;
pub const norm = @import("linalg/norm.zig").norm;
pub const normAlloc = @import("linalg/norm.zig").normAlloc;

// hadamard (also for vectors)

pub const Matmul = @import("linalg/matmul.zig").Matmul;
pub const matmul = @import("linalg/matmul.zig").matmul;
pub const matmulUnchecked = @import("linalg/matmul.zig").matmulUnchecked;
pub const matmulAlloc = @import("linalg/matmul.zig").matmulAlloc;
pub const matmulInto = @import("linalg/matmul.zig").matmulInto;
pub const matmulIntoAlloc = @import("linalg/matmul.zig").matmulIntoAlloc;
pub const matmulIntoUnchecked = @import("linalg/matmul.zig").matmulIntoUnchecked;

// pub const Solve = @import("linalg/solve.zig").Solve;
pub const SolveMethod = @import("linalg/solve.zig").SolveMethod;
// solve (with SolveMethod, to choose decomposition, and optional out parameter to save it)
// Allow vector and matrix rhs

// pub const Lstsq = @import("linalg/lstsq.zig").Lstsq;
pub const LstsqMethod = @import("linalg/lstsq.zig").LstsqMethod;
// lstsq
// Allow vector and matrix rhs

// Matrix operations
pub const Trace = @import("linalg/trace.zig").Trace;
pub const trace = @import("linalg/trace.zig").trace;
pub const traceUnchecked = @import("linalg/trace.zig").traceUnchecked;

// det

// cond (condition number)

// gram (gramian matrix, Xᴴ X)

// cong (matrix congruence, A X Aᴴ)

// kron (A ⊗ B; (m × n) ⊗ (p × q) -> (mp × nq))

// expm

// sqrtm

// inv

// pinv

// pub const lu = @import("linalg/lu.zig");
pub const plu = @import("linalg/plu.zig");
// pub const pluq = @import("linalg/pluq.zig");

// ilu

pub const llt = @import("linalg/llt.zig");
pub const utu = @import("linalg/utu.zig");

// pub const ldlt = @import("linalg/ldlt.zig");
// pub const udut = @import("linalg/udut.zig");

// pub const qr = @import("linalg/qr.zig");
// pub const qrp = @import("linalg/qrp.zig");

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
    ZeroDimension,
};
