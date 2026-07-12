pub const meta = @import("meta.zig");

pub const int = @import("int.zig");
pub const float = @import("float.zig");
pub const dyadic = @import("dyadic.zig");
pub const Dyadic = dyadic.Dyadic;
pub const complex = @import("complex.zig");
pub const Complex = complex.Complex;
pub const cf16 = complex.cf16;
pub const cf32 = complex.cf32;
pub const cf64 = complex.cf64;
pub const cf80 = complex.cf80;
pub const cf128 = complex.cf128;
pub const comptime_complex = complex.comptime_complex;

// Domain namespaces
pub const numeric = @import("numeric.zig");
pub const vector = @import("vector.zig");
pub const matrix = @import("matrix.zig");
pub const array = @import("array.zig");
pub const poly = @import("poly.zig");

// Module namespaces
pub const stats = @import("stats.zig");
pub const linalg = @import("linalg.zig");
pub const autodiff = @import("autodiff.zig");
