const std = @import("std");
const zml = @import("../zml.zig");

pub const cf16 = cf(f16);
pub const cf32 = cf(f32);
pub const cf64 = cf(f64);
pub const cf128 = cf(f128);
// Maybe keep like this, maybe do all four separate so the arguments of the
// functions are specific (fxx and not T).
fn cf(comptime T: type) type {
    return struct {
        const Self = @This();
        re: T,
        im: T,

        pub fn init(re: T, im: T) Self {
            return Self{
                .re = re,
                .im = im,
            };
        }
    };
}
