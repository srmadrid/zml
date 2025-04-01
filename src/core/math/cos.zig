pub fn cos(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @cos(x),
        f32 => return @cos(x),
        f64 => return @cos(x),
        f80 => return @cos(x),
        f128 => return @cos(x),
        else => @compileError("x must be a float"),
    }
}
