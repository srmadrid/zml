const std = @import("std");

pub const supported = @import("supported.zig");

pub const CamelError = error{
    NullInput,
};
