const std = @import("std");

pub const Integer = struct {
    limbs: []usize,
    metadata: usize,

    pub const empty: Integer = .{
        .limbs = &.{},
        .metadata = 0,
    };
};
