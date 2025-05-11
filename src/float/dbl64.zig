const ieee_f64_shape = packed struct {
    lsw: u32,
    msw: u32,
};

/// Get two 32 bit ints from a double.
pub inline fn extractWords(ix0: anytype, ix1: anytype, d: f64) void {
    const ew_u: ieee_f64_shape = @bitCast(d);
    ix0.* = @bitCast(ew_u.msw);
    ix1.* = @bitCast(ew_u.lsw);
}

/// Get the more significant 32 bit int from a double.
pub inline fn getHighWord(i: anytype, d: f64) void {
    const gh_u: ieee_f64_shape = @bitCast(d);
    i.* = @bitCast(gh_u.msw);
}

/// Get the less significant 32 bit int from a double.
pub inline fn getLowWord(i: anytype, d: f64) void {
    const gl_u: ieee_f64_shape = @bitCast(d);
    i.* = @bitCast(gl_u.lsw);
}

/// Get all in one, efficient on 64-bit machines.
pub inline fn extractWords64(i: anytype, d: f64) void {
    const gh_u: u64 = @bitCast(d);
    i.* = @bitCast(gh_u);
}

/// Set a double from two 32 bit ints.
pub inline fn insertWords(d: *f64, ix0: anytype, ix1: anytype) void {
    const iw_u: ieee_f64_shape = .{ .msw = @bitCast(ix0), .lsw = @bitCast(ix1) };
    d.* = @bitCast(iw_u);
}

/// Get all in one, efficient on 64-bit machines.
pub inline fn insertWords64(d: *f64, i: anytype) void {
    const iw_u: u64 = @bitCast(i);
    d.* = @bitCast(iw_u);
}

/// Set the more significant 32 bits of a double from an int.
pub inline fn setHighWord(d: *f64, v: anytype) void {
    var sh_u: ieee_f64_shape = @bitCast(d.*);
    sh_u.msw = @bitCast(v);
    d.* = @bitCast(sh_u);
}

/// Set the less significant 32 bits of a double from an int.
pub inline fn setLowWord(d: *f64, v: anytype) void {
    var sl_u: ieee_f64_shape = @bitCast(d.*);
    sl_u.lsw = @bitCast(v);
    d.* = @bitCast(sl_u);
}
