const std = @import("std");

pub const Integer = struct {
    limbs: []u32,
    size: u32,
    positive: bool,

    pub const empty: Integer = .{
        .limbs = &.{},
        .size = 0,
        .positive = true,
    };
};
