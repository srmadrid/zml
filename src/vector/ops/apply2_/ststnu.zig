pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    inline for (0..@TypeOf(x).len) |i| {
        op_(&o.data[i], x.data[i], y);
    }
}
