const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    inline for (0..(comptime int.min(O.rows, O.cols))) |i| {
        opInto(
            &o.data[i],
            if (comptime meta.isMatrix(X))
                x.get(i, i) catch unreachable
            else
                x,
            if (comptime meta.isMatrix(Y))
                y.get(i, i) catch unreachable
            else
                y,
        );
    }
}
