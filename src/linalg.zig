pub const blas = @import("linalg/blas.zig");
pub const lapack = @import("linalg/lapack.zig");

pub const Order = enum(c_uint) {
    row_major = 101,
    col_major = 102,
};

pub const Transpose = enum(c_uint) {
    no_trans = 111,
    trans = 112,
    conj_trans = 113,
    conj_no_trans = 114,

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

pub const Uplo = enum(c_uint) {
    upper = 121,
    lower = 122,

    pub inline fn toChar(self: Uplo) u8 {
        return switch (self) {
            .upper => 'U',
            .lower => 'L',
        };
    }

    pub inline fn invert(self: Uplo) Uplo {
        return switch (self) {
            .upper => .lower,
            .lower => .upper,
        };
    }
};

pub const Diag = enum(c_uint) {
    non_unit = 131,
    unit = 132,

    pub inline fn invert(self: Diag) Diag {
        return switch (self) {
            .non_unit => .unit,
            .unit => .non_unit,
        };
    }
};

pub const Side = enum(c_uint) {
    left = 141,
    right = 142,

    pub inline fn invert(self: Side) Side {
        return switch (self) {
            .left => .right,
            .right => .left,
        };
    }
};
