pub fn log10(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @log10(x),
        f32 => return @log10(x),
        f64 => return @log10(x),
        f80 => return @log10(x),
        f128 => return @log10(x),
        else => @compileError("x must be a float"),
    }
}
