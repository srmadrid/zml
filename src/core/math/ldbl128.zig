pub const ieee_f128_shape64 = packed struct {
    lsw: u64,
    msw: u64,
};

pub const ieee_f128_shape32 = packed struct {
    w3: u32,
    w2: u32,
    w1: u32,
    w0: u32,
};

/// Get two 64 bit ints from a long double.
pub inline fn getWords(ix0: anytype, ix1: anytype, d: f128) void {
    const qw_u: ieee_f128_shape64 = @bitCast(d);
    ix0.* = @bitCast(qw_u.msw);
    ix1.* = @bitCast(qw_u.lsw);
}

/// Set a long double from two 64 bit ints.
pub inline fn setWords(d: *f128, ix0: anytype, ix1: anytype) void {
    const qw_u: ieee_f128_shape64 = .{ .msw = @bitCast(ix0), .lsw = @bitCast(ix1) };
    d.* = @bitCast(qw_u);
}

/// Get the more significant 64 bits of a long double mantissa.
pub inline fn getMsw(v: anytype, d: f128) void {
    const sh_u: ieee_f128_shape64 = @bitCast(d);
    v.* = @bitCast(sh_u.msw);
}

/// Set the more significant 64 bits of a long double mantissa from an int.
pub inline fn setMsw(d: *f128, v: anytype) void {
    var sh_u: ieee_f128_shape64 = @bitCast(d.*);
    sh_u.msw = @bitCast(v);
    d.* = @bitCast(sh_u);
}

/// Get the least significant 64 bits of a long double mantissa.
pub inline fn getLsw(v: anytype, d: f128) void {
    const sh_u: ieee_f128_shape64 = @bitCast(d);
    v.* = @bitCast(sh_u.lsw);
}
