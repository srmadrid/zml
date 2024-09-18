//! The symbol interface.

const std = @import("../std.zig");
const Symbol = @This();

// The type erased pointer to the symbol implementation.
ptr: *anyopaque,
vtable: *const VTable,

pub const VTable = struct {
    /// Some function.
    someFunc: *const fn (ctx: *anyopaque, some_data: usize) bool,
};
