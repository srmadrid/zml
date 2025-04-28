pub const ieee_f32_shape = packed struct {
    man: u23,
    exp: u8,
    sign: u1,
};

/// Get a 32 bit int from a float.
pub inline fn getWord(i: anytype, d: f32) void {
    const gf_u: u32 = @bitCast(d);
    i.* = @bitCast(gf_u);
}

/// Set a float from a 32 bit int.
pub inline fn setWord(d: *f32, i: anytype) void {
    const sf_u: u32 = @bitCast(i);
    d.* = @bitCast(sf_u);
}
