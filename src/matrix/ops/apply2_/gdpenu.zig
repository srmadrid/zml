const meta = @import("../../../meta.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);

    o.setAll(numeric.zero(meta.Numeric(O)));

    var k: usize = 0;
    while (k < x.rows) : (k += 1) {
        const i = if (x.direction == .forward) k else x.data[k];
        const j = if (x.direction == .forward) x.data[k] else k;

        if (comptime op_ == numeric.mulInto)
            numeric.set(&o.data[o._index(i, j)], y)
        else
            op_(&o.data[o._index(i, j)], numeric.one(meta.Numeric(X)), y);
    }
}
