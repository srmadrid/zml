const numeric = @import("../../numeric.zig");

const matops = @import("../ops.zig");

pub fn mul_(o: anytype, x: anytype, y: anytype) !void {
    return matops.apply2_(o, x, y, numeric.mulInto);
}
