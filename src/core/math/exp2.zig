pub fn exp2(x: anytype) @TypeOf(x) {
    switch (@TypeOf(x)) {
        f16 => return @exp2(x),
        f32 => return @exp2(x),
        f64 => return @exp2(x),
        f80 => return @exp2(x),
        f128 => return @exp2(x),
        else => @compileError("x must be a float"),
    }
}
