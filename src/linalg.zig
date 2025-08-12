pub const blas = @import("linalg/blas.zig");
pub const lapack = @import("linalg/lapack.zig");

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
