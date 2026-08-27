const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const linalg = @import("../linalg.zig");

pub fn SolveMethod(A: type) type {
    if (!meta.isMatrix(A) or meta.isBuilderMatrix(A))
        @compileError("zsl.linalg.SolveMethod: A must be a non-builder matrix type, got \n\tA = " ++ @typeName(A) ++ "\n");

    return union(enum) {
        // Structural solvers
        /// Triangular solver. For triangular matrices.
        tri,
        /// Diagonal solver. For diagonal matrices.
        dia,
        /// Permutation solver. For permutation matrices.
        perm,

        // Direct solvers
        /// LU decomposition, no pivoting. For general matrices.
        lu: struct {
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.lu.Factor(meta.Numeric(A)),
        },
        /// LU decomposition, row pivoting. For general matrices.
        plu: struct {
            /// Pointer to save the factorization of the matrix during the solve.
            save: ?*linalg.plu.Factor(meta.Numeric(A)),
        },
        /// LU decomposition, row and column pivoting. For general matrices.
        pluq: struct {
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.pluq.Factor(meta.Numeric(A)),
        },
        /// Lower Cholesky decomposition. For lower positive-definite Hermitian matrices.
        llt: struct {
            /// Pointer to save the factorization of the matrix during the solve.
            save: ?*linalg.llt.Factor(meta.Numeric(A)),
        },
        /// Upper Cholesky decomposition. For upper positive-definite Hermitian matrices.
        utu: struct {
            /// Pointer to save the factorization of the matrix during the solve.
            save: ?*linalg.utu.Factor(meta.Numeric(A)),
        },
        /// Lower Bunch-Kaufman decomposition. For lower Hermitian matrices.
        ldlt: struct {
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.ldlt.Factor(meta.Numeric(A)),
        },
        /// Upper Bunch-Kaufman decomposition. For upper Hermitian matrices.
        udut: struct {
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.udut.Factor(meta.Numeric(A)),
        },
        /// QR decomposition, no pivoting. For general matrices.
        qr: struct {
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.qr.Factor(meta.Numeric(A)),
            // maybe also workspace arrays (optional)
        },
        /// QR decomposition, column pivoting. For general matrices.
        qrp: struct {
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.qrp.Factor(meta.Numeric(A)),
            // maybe also workspace arrays (optional)
        },
        /// SVD decomposition. For general matrices.
        svd: struct {
            /// Singular values smaller than rcond * maxᵢ σᵢ treated as 0.
            rcond: meta.Numeric(A) = numeric.cast(meta.Numeric(A), 1e-15),
            // /// Pointer to save the factorization of the matrix during the solve.
            // save: ?*linalg.svd.Factor(meta.Numeric(A)),
            // maybe also workspace arrays (optional)
        },

        // Krylov subspace iterative solvers
        /// Conjugate gradient method. For positive-definite Hermitian matrices.
        cg: struct {
            /// Relative residual tolerance required for convergence.
            tol: meta.Numeric(A) = numeric.cast(meta.Numeric(A), 1e-6),
            /// Maximum number of iterations allowed before returning an error.
            max_iter: usize = 1_000,
            // preconditioner
        },
        /// Biconjugate gradient stabilized method. For general matrices.
        bicgstab: struct {
            /// Relative residual tolerance required for convergence.
            tol: meta.Numeric(A) = numeric.cast(meta.Numeric(A), 1e-6),
            /// Maximum number of iterations allowed before returning an error.
            max_iter: usize = 1_000,
            // preconditioner
        },
        /// Minimal residual method. For Hermitian matrices.
        minres: struct {
            /// Relative residual tolerance required for convergence.
            tol: meta.Numeric(A) = numeric.cast(meta.Numeric(A), 1e-6),
            /// Maximum number of iterations allowed before returning an error.
            max_iter: usize = 1_000,
            // preconditioner
        },
        /// Generalized minimal residual method. For general matrices.
        gmres: struct {
            /// Relative residual tolerance required for convergence.
            tol: meta.Numeric(A) = numeric.cast(meta.Numeric(A), 1e-6),
            /// Maximum number of iterations allowed before returning an error.
            max_iter: usize = 1_000,
            restart: usize = 30,
            // preconditioner
        },

        pub const default: linalg.SolveMethod(A) = switch (meta.matrixType(A)) {
            .general_static, .general_dense, .general_sparse => .{ .plu = .{} },
            .symmetric_static, .symmetric_dense, .symmetric_sparse => if (meta.uploOf(A) == .upper) .{ .udut = .{} } else .{ .ldlt = .{} },
            .hermitian_static, .hermitian_dense, .hermitian_sparse => if (meta.uploOf(A) == .upper) .{ .udut = .{} } else .{ .ldlt = .{} },
            .triangular_static, .triangular_dense, .triangular_sparse => .tri,
            .diagonal_static, .diagonal_sparse => .dia,
            .permutation_static, .permutation_sparse => .perm,
            else => unreachable,
        };
    };
}
