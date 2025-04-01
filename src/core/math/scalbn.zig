const std = @import("std");

pub fn scalbn(x: anytype, n: anytype) @TypeOf(x) {
    return std.math.scalbn(x, n);
}
