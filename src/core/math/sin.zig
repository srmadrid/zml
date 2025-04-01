pub fn sin(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @sin(x),
        f32 => return @sin(x),
        f64 => return @sin(x),
        f80 => return @sin(x),
        f128 => return @sin(x),
        else => @compileError("x must be a float"),
    }
}
