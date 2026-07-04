//! Namespace for automatic differentiation types and operations.

// Forward mode
pub const dual = @import("autodiff/dual.zig");
pub const Dual = dual.Dual;

// Backward mode
pub const Op = enum {
    // Leaf
    @"var",

    // Basic operations
    abs,
    abs1,
    abs2,
    neg,
    re,
    im,
    conj,
    sign,

    // Arithmetic operations
    add,
    mul,
    sub,
    div,

    // Comparison operations
    min,
    max,

    // Exponential functions
    exp,
    ln,
};

pub const Tape = @import("autodiff/tape.zig").Tape;
pub const @"var" = @import("autodiff/var.zig");
pub const Var = @"var".Var;
