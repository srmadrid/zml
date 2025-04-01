pub fn log2(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @log2(x),
        f32 => return @log2(x),
        f64 => return @log2(x),
        f80 => return @log2(x),
        f128 => return @log2(x),
        else => @compileError("x must be a float"),
    }
}
