pub fn log(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @log(x),
        f32 => return @log(x),
        f64 => return @log(x),
        f80 => return @log(x),
        f128 => return @log(x),
        else => @compileError("x must be a float"),
    }
}
