const std = @import("std");

pub const blas = @import("linalg/blas.zig");
//pub const lapack = @import("linalg/lapack.zig");

test {
    std.testing.refAllDeclsRecursive(@This());
}
