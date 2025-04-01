pub fn exp(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @exp(x),
        f32 => return @exp(x),
        f64 => return @exp(x),
        f80 => return @exp(x),
        f128 => return @exp(x),
        else => @compileError("x must be a float"),
    }
}
