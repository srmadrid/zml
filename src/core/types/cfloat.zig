const numericType = @import("../types.zig").numericType;

pub fn cfloat(comptime T: type) type {
    if (numericType(T) != .float) @compileError("Unsupported type for cfloat: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,
    };
}

pub const cf16 = cfloat(f16);
pub const cf32 = cfloat(f32);
pub const cf64 = cfloat(f64);
pub const cf80 = cfloat(f80);
pub const cf128 = cfloat(f128);
pub const comptime_complex = cfloat(comptime_float);
