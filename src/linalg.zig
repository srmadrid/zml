const std = @import("std");

const types = @import("types.zig");
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const MulCoerce = types.MulCoerce;
const int = @import("int.zig");
const ops = @import("ops.zig");

pub const blas = @import("linalg/blas.zig");
pub const lapack = @import("linalg/lapack.zig");

pub inline fn dot(x: anytype, y: anytype, ctx: anytype) !Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y))) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const C: type = Coerce(Numeric(X), Numeric(Y));

    comptime if (!types.isVector(X) or !types.isVector(Y))
        @compileError("dot: both arguments must be vectors, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.dotc not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (x.len != y.len)
        return Error.DimensionMismatch;

    if (comptime types.isDenseVector(X)) {
        if (comptime types.isDenseVector(Y)) {
            if (comptime types.isComplex(C)) {
                return blas.dotc(types.scast(i32, x.len), x.data, x.inc, y.data, y.inc, ctx);
            } else {
                return blas.dot(types.scast(i32, x.len), x.data, x.inc, y.data, y.inc, ctx);
            }
        } else {
            return Error.NotImplemented;
        }
    } else {
        if (comptime types.isDenseVector(Y)) {
            return Error.NotImplemented;
        } else {
            return Error.NotImplemented;
        }
    }
}

pub const matmul = @import("linalg/matmul.zig").matmul;

// const _lu = @import("linalg/lu.zig");
// pub const LU = _lu.LU;
// pub const lu = _lu.lu;
// pub const PLU = _lu.PLU;
// pub const plu = _lu.plu;
// pub const PLUQ = _lu.PLUQ;
// pub const pluq = _lu.pluq;

const _cholesky = @import("linalg/cholesky.zig");
pub const LLT = _cholesky.LLT;
pub const llt = _cholesky.llt;
pub const UTU = _cholesky.UTU;
pub const utu = _cholesky.utu;
pub const cholesky = _cholesky.cholesky;

const _bunchkaufman = @import("linalg/bunchkaufman.zig");
pub const LDLT = _bunchkaufman.LDLT;
pub const ldlt = _bunchkaufman.ldlt;
pub const UDUT = _bunchkaufman.UDUT;
pub const udut = _bunchkaufman.udut;
pub const bunchkaufman = _bunchkaufman.bunchkaufman;

const qr_ = @import("linalg/qr.zig");
pub const QR = qr_.QR;
pub const qr = qr_.qr;
pub const QRP = qr_.QRP;
pub const qrp = qr_.qrp;

pub const Transpose = enum(u2) {
    no_trans,
    trans,
    conj_trans,
    conj_no_trans,

    pub inline fn toCUInt(self: Transpose) c_uint {
        return switch (self) {
            .no_trans => 111,
            .trans => 112,
            .conj_trans => 113,
            .conj_no_trans => 114,
        };
    }

    pub inline fn toCInt(self: Transpose) c_int {
        return switch (self) {
            .no_trans => 111,
            .trans => 112,
            .conj_trans => 113,
            .conj_no_trans => 114,
        };
    }

    pub inline fn toChar(self: Transpose) u8 {
        return switch (self) {
            .no_trans => 'N',
            .trans => 'T',
            .conj_no_trans => 0,
            .conj_trans => 'C',
        };
    }

    pub inline fn invert(self: Transpose) Transpose {
        return switch (self) {
            .no_trans => .trans,
            .trans => .no_trans,
            .conj_no_trans => .conj_trans,
            .conj_trans => .conj_no_trans,
        };
    }

    pub inline fn reverse(self: Transpose) Transpose {
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

    pub inline fn toCUInt(self: Side) c_uint {
        return switch (self) {
            .left => 141,
            .right => 142,
        };
    }

    pub inline fn toCInt(self: Side) c_int {
        return switch (self) {
            .left => 141,
            .right => 142,
        };
    }

    pub inline fn toChar(self: Side) u8 {
        return switch (self) {
            .left => 'L',
            .right => 'R',
        };
    }

    pub inline fn invert(self: Side) Side {
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
