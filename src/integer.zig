const std = @import("std");

pub const Integer = struct {
    limbs: []usize,
    size: usize,
    positive: bool,

    pub const empty: Integer = .{
        .limbs = &.{},
        .size = 0,
        .positive = true,
    };
};
