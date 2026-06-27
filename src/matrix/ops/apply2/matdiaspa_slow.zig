const int = @import("../../../int.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    var i: usize = 0;
    while (i < int.min(o.rows, o.cols)) : (i += 1) {
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

        if (!o.flags.noconj)
            numeric.conjInto(&o.data[i], o.data[i]);
    }
}
