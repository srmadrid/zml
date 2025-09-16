const types = @import("../types.zig");
const Order = types.Order;

pub inline fn index(order: Order, i: i32, j: i32, ld: i32) u32 {
    return types.scast(u32, if (order == .col_major) i + j * ld else i * ld + j);
}

// Comptime version
pub inline fn cindex(comptime order: Order, i: u32, j: u32, ld: u32) u32 {
    return if (comptime order == .col_major) i + j * ld else i * ld + j;
}

pub inline fn col_ld(order: Order, ld: i32) i32 {
    return if (order == .col_major) 1 else ld;
}

pub inline fn row_ld(order: Order, ld: i32) i32 {
    return if (order == .col_major) ld else 1;
}
