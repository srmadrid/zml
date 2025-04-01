pub fn sqrt(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @sqrt(x),
        f32 => return @sqrt(x),
        f64 => return @sqrt(x),
        f80 => return @sqrt(x),
        f128 => return @sqrt(x),
        else => @compileError("x must be a float"),
    }
}
