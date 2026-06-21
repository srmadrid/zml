pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    inline for (0..@TypeOf(y).len) |i| {
        opInto(&o.data[i], x, y.data[i]);
    }
}
