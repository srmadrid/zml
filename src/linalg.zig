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
    const C: type = Coerce(Numeric(@TypeOf(x)), Numeric(@TypeOf(y)));

    comptime if (!types.isVector(@TypeOf(x)) or !types.isVector(@TypeOf(y)))
        @compileError("dot: both arguments must be vectors, got " ++ @typeName(@TypeOf(x)) ++ " and " ++ @typeName(@TypeOf(y)));

    comptime if (types.isArbitraryPrecision(C)) {
        @compileError("zml.linalg.blas.dotc not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    if (x.len != y.len)
        return Error.DimensionMismatch;

    if (comptime types.isComplex(C)) {
        return blas.dotc(types.scast(i32, x.len), x.data, x.inc, y.data, y.inc, ctx);
    } else {
        return blas.dot(types.scast(i32, x.len), x.data, x.inc, y.data, y.inc, ctx);
    }
}

pub const matmul = @import("linalg/matmul.zig").matmul;

const lu_ = @import("linalg/lu.zig");
pub const LU = lu_.LU;
pub const lu = lu_.lu;
pub const PLU = lu_.PLU;
pub const plu = lu_.plu;
pub const PLUQ = lu_.PLUQ;
pub const pluq = lu_.pluq;

const cholesky_ = @import("linalg/cholesky.zig");
pub const LLT = cholesky_.LLT;
pub const llt = cholesky_.llt;
pub const UTU = cholesky_.UTU;
pub const utu = cholesky_.utu;
pub const cholesky = cholesky_.cholesky;

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
    SingularMatrix,
    IndefiniteMatrix,
    NotImplemented,
};
