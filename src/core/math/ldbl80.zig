const ieee_f80_shape = packed struct {
    lsw: u32,
    msw: u32,
    sign_exponent: i16,
    empty: u16,
};

/// Get three 32 bit ints from a double.
pub inline fn getWords(exp: anytype, ix0: anytype, ix1: anytype, d: f80) void {
    const ew_u: ieee_f80_shape = @bitCast(d);
    exp.* = @bitCast(ew_u.sign_exponent);
    ix0.* = @bitCast(ew_u.msw);
    ix1.* = @bitCast(ew_u.lsw);
}

/// Set a double from two 32 bit ints.
pub inline fn setWords(d: *f80, exp: anytype, ix0: anytype, ix1: anytype) void {
    const iw_u: ieee_f80_shape = .{ .sign_exponent = @bitCast(exp), .msw = @bitCast(ix0), .lsw = @bitCast(ix1) };
    d.* = @bitCast(iw_u);
}

/// Get the more significant 32 bits of a long double mantissa.
pub inline fn getMsw(v: anytype, d: f80) void {
    const sh_u: ieee_f80_shape = @bitCast(d);
    v.* = @bitCast(sh_u.msw);
}

/// Set the more significant 32 bits of a long double mantissa from an int.
pub inline fn setMsw(d: *f80, v: anytype) void {
    var sh_u: ieee_f80_shape = @bitCast(d.*);
    sh_u.msw = @bitCast(v);
    d.* = @bitCast(sh_u);
}

/// Get int from the exponent of a long double.
pub inline fn getExp(exp: anytype, d: f80) void {
    const ge_u: ieee_f80_shape = @bitCast(d);
    exp.* = @bitCast(ge_u.sign_exponent);
}

/// Set exponent of a long double from an int.
pub inline fn setExp(d: *f80, exp: anytype) void {
    var se_u: ieee_f80_shape = @bitCast(d.*);
    se_u.sign_exponent = @bitCast(exp);
    d.* = @bitCast(se_u);
}
