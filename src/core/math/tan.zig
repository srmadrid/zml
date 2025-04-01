pub fn tan(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @tan(x),
        f32 => return @tan(x),
        f64 => return @tan(x),
        f80 => return @tan(x),
        f128 => return @tan(x),
        else => @compileError("x must be a float"),
    }
}
