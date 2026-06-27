//! Namespace for matrix types and operations.

const meta = @import("meta.zig");

pub const general = @import("matrix/general.zig");
pub const symmetric = @import("matrix/symmetric.zig");
pub const hermitian = @import("matrix/hermitian.zig");
pub const triangular = @import("matrix/triangular.zig");
pub const diagonal = @import("matrix/diagonal.zig");
pub const permutation = @import("matrix/permutation.zig");

pub const builder = @import("matrix/builder.zig");

const matops = @import("matrix/ops.zig");
pub const Add = matops.Add;
pub const add = matops.add;
pub const addAlloc = matops.addAlloc;
pub const addInto = matops.addInto;
pub const addIntoUnchecked = matops.addIntoUnchecked;
pub const Sub = matops.Sub;
pub const sub = matops.sub;
pub const subAlloc = matops.subAlloc;
pub const subInto = matops.subInto;
pub const subIntoUnchecked = matops.subIntoUnchecked;
pub const Mul = matops.Mul;
pub const mul = matops.mul;
pub const mulAlloc = matops.mulAlloc;
pub const mulInto = matops.mulInto;
pub const mulIntoUnchecked = matops.mulIntoUnchecked;
pub const Div = matops.Div;
pub const div = matops.div;
pub const divAlloc = matops.divAlloc;
pub const divInto = matops.divInto;
pub const divIntoUnchecked = matops.divIntoUnchecked;

pub const Layout = enum(u1) {
    row_major,
    col_major,

    pub const default: Layout = .col_major;

    pub fn toInt(self: Layout, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.matrix.Layout.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .row_major => 101,
            .col_major => 102,
        };
    }

    pub fn invert(self: Layout) Layout {
        return switch (self) {
            .row_major => .col_major,
            .col_major => .row_major,
        };
    }
};

pub const Uplo = enum(u1) {
    upper,
    lower,

    pub const default: Uplo = .upper;

    pub fn toInt(self: Uplo, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.matrix.Uplo.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .upper => if (comptime Int == u8) 'U' else 121,
            .lower => if (comptime Int == u8) 'L' else 122,
        };
    }

    pub fn invert(self: Uplo) Uplo {
        return switch (self) {
            .upper => .lower,
            .lower => .upper,
        };
    }
};

pub const Diag = enum(u1) {
    non_unit,
    unit,

    pub const default: Diag = .non_unit;

    pub fn toInt(self: Diag, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.matrix.Diag.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .non_unit => if (comptime Int == u8) 'N' else 131,
            .unit => if (comptime Int == u8) 'U' else 132,
        };
    }

    pub fn invert(self: Diag) Diag {
        return switch (self) {
            .non_unit => .unit,
            .unit => .non_unit,
        };
    }
};

pub const Error = error{
    ZeroDimension,
    PositionOutOfBounds,
    BreaksStructure,
    InvalidRange,
    DimensionMismatch,
    InvalidBandwidth,
    NotSquare,
    DataNotOwned,
    InsufficientSpace,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    noconj: bool = true,
};
