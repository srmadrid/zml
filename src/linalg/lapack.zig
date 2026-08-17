const meta = @import("../meta.zig");

pub const Direction = enum {
    forward,
    backward,

    pub fn toInt(self: Direction, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.linalg.lapack.Direction.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .forward => 'F',
            .backward => 'B',
        };
    }
};

pub const Storage = enum {
    columnwise,
    rowwise,

    pub fn toInt(self: Storage, comptime Int: type) Int {
        comptime if (!meta.isNumeric(Int) or meta.numericType(Int) != .int)
            @compileError("zsl.linalg.lapack.Storage.toInt: Int must be an int type, got:\n\tInt = " ++ @typeName(Int) ++ "\n");

        return switch (self) {
            .columnwise => 'C',
            .rowwise => 'R',
        };
    }
};

// LAMCH cheatsheet:
// key  (u8) -> value (value when N is integral)
// eps   (E) -> numeric.eps (0)
// sfmin (S) -> numeric.smallest (1)
// base  (B) -> 2
// prec  (P) -> numeric.eps * 2
// t     (N) -> *not used*
// rnd   (R) -> *not used*
// emin  (M) -> *not used*
// rmin  (U) -> numeric.smallest (1)
// emax  (L) -> *not used*
// rmax  (O) -> numeric.highest (highest)

// LU
pub const getrf2 = @import("lapack/lu/getrf2.zig").getrf2; // Unoptimized
pub const getrf = @import("lapack/lu/getrf.zig").getrf; // Unoptimized
pub const getrs = @import("lapack/lu/getrs.zig").getrs; // Unoptimized
pub const gesv = @import("lapack/lu/gesv.zig").gesv; // Unoptimized

// Cholesky
pub const potrf2 = @import("lapack/cholesky/potrf2.zig").potrf2; // Unoptimized
pub const potrf = @import("lapack/cholesky/potrf.zig").potrf; // Unoptimized
pub const potrs = @import("lapack/cholesky/potrs.zig").potrs; // Unoptimized
pub const posv = @import("lapack/cholesky/posv.zig").posv; // Unoptimized

// Bunch-Kaufman
// pub const sytf2 = @import("lapack/bunchkaufman/sytf2.zig").sytf2;
// pub const sytrf = @import("lapack/bunchkaufman/sytrf.zig").sytrf;
// pub const hetf2 = @import("lapack/bunchkaufman/hetf2.zig").hetf2;
// pub const hetrf = @import("lapack/bunchkaufman/hetrf.zig").hetrf;

// QR
pub const geqr2 = @import("lapack/qr/geqr2.zig").geqr2; // Unoptimized
// pub const geqrf = @import("lapack/qr/geqrf.zig").geqrf;
// pub const org2r = @import("lapack/qr/org2r.zig").org2r;
// pub const orgqr = @import("lapack/qr/orgqr.zig").orgqr;
// pub const ung2r = @import("lapack/qr/ung2r.zig").ung2r;
// pub const ungqr = @import("lapack/qr/ungqr.zig").ungqr;

// Auxiliary functions
pub const lacgv = @import("lapack/aux/lacgv.zig").lacgv; // Unoptimized
pub const ilalc = @import("lapack/aux/ilalc.zig").ilalc; // Unoptimized
pub const ilalr = @import("lapack/aux/ilalr.zig").ilalr; // Unoptimized
pub const lacpy = @import("lapack/aux/lacpy.zig").lacpy; // Unoptimized
pub const laswp = @import("lapack/aux/laswp.zig").laswp; // Unoptimized
// pub const lasyf = @import("lapack/aux/lasyf.zig").lasyf;
// pub const lahef = @import("lapack/aux/lahef.zig").lahef;
pub const larfg = @import("lapack/aux/larfg.zig").larfg; // Unoptimized
pub const larf1f = @import("lapack/aux/larf1f.zig").larf1f; // Unoptimized
pub const larft = @import("lapack/aux/larft.zig").larft; // Unoptimized
// pub const larfb = @import("lapack/aux/larfb.zig").larfb; // Unoptimized

pub const Error = error{
    InvalidArgument,
};
