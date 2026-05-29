const meta = @import("meta.zig");

const int = @import("int.zig");

pub const cblas = @import("linalg/cblas.zig");
pub const blas = @import("linalg/blas.zig");
// pub const lapacke = @import("linalg/lapacke.zig");
// pub const lapack = @import("linalg/lapack.zig");

// pub fn dot(x: anytype, y: anytype, ctx: anytype) !Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y))) {
//     // Dot(X, Y) = begin
//     //     const T = numeric.Mul(meta.Numeric(X), meta.Numeric(Y))
//     //     return numeric.Add(T, T)
//     // end
//     const X: type = @TypeOf(x);
//     const Y: type = @TypeOf(y);
//     const C: type = Coerce(Numeric(X), Numeric(Y));

//     comptime if (!meta.isVector(X) or !meta.isVector(Y))
//         @compileError("dot: both arguments must be vectors, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

//     comptime if (meta.isArbitraryPrecision(C)) {
//         @compileError("zsl.linalg.blas.dotc not implemented for arbitrary precision types yet");
//     } else {
//         meta.validateContext(@TypeOf(ctx), .{});
//     };

//     if (x.len != y.len)
//         return Error.DimensionMismatch;

//     if (comptime meta.isDenseVector(X)) {
//         if (comptime meta.isDenseVector(Y)) {
//             return blas.dotc(meta.scast(i32, x.len), x.data, x.inc, y.data, y.inc, ctx);
//         } else {
//             return Error.NotImplemented;
//         }
//     } else {
//         if (comptime meta.isDenseVector(Y)) {
//             return Error.NotImplemented;
//         } else {
//             return Error.NotImplemented;
//         }
//     }
// }

// pub const matmul = @import("linalg/matmul.zig").matmul;

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
