const meta = @import("../../../meta.zig");

const int = @import("../../../int.zig");

const numeric = @import("../../../numeric.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            if (comptime meta.uploOf(O) == .upper) {
                var i: usize = 0;
                while (i < int.min(j, o.rows)) : (i += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                }

                if (j < o.rows) {
                    if (comptime meta.diagOf(X) == .unit) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                o.data[o._index(j, j)] = numeric.two(meta.Numeric(O))
                            else
                                o.data[o._index(j, j)] = numeric.zero(meta.Numeric(O));
                        } else {
                            op_(&o.data[o._index(j, j)], numeric.one(meta.Numeric(X)), y.data[y._index(j, j)]);
                        }
                    } else {
                        if (comptime meta.diagOf(Y) == .unit)
                            op_(&o.data[o._index(j, j)], x.data[x._index(j, j)], numeric.one(meta.Numeric(Y)))
                        else
                            op_(&o.data[o._index(j, j)], x.data[x._index(j, j)], y.data[y._index(j, j)]);
                    }
                }
            } else {
                if (j < o.rows) {
                    if (comptime meta.diagOf(X) == .unit) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                o.data[o._index(j, j)] = numeric.two(meta.Numeric(O))
                            else
                                o.data[o._index(j, j)] = numeric.zero(meta.Numeric(O));
                        } else {
                            op_(&o.data[o._index(j, j)], numeric.one(meta.Numeric(X)), y.data[y._index(j, j)]);
                        }
                    } else {
                        if (comptime meta.diagOf(Y) == .unit)
                            op_(&o.data[o._index(j, j)], x.data[x._index(j, j)], numeric.one(meta.Numeric(Y)))
                        else
                            op_(&o.data[o._index(j, j)], x.data[x._index(j, j)], y.data[y._index(j, j)]);
                    }
                }

                var i: usize = int.min(j + 1, o.rows);
                while (i < o.rows) : (i += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            if (comptime meta.uploOf(O) == .lower) {
                var j: usize = 0;
                while (j < int.min(i, o.cols)) : (j += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                }

                if (i < o.cols) {
                    if (comptime meta.diagOf(X) == .unit) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                o.data[o._index(i, i)] = numeric.two(meta.Numeric(O))
                            else
                                o.data[o._index(i, i)] = numeric.zero(meta.Numeric(O));
                        } else {
                            op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y.data[y._index(i, i)]);
                        }
                    } else {
                        if (comptime meta.diagOf(Y) == .unit)
                            op_(&o.data[o._index(i, i)], x.data[x._index(i, i)], numeric.one(meta.Numeric(Y)))
                        else
                            op_(&o.data[o._index(i, i)], x.data[x._index(i, i)], y.data[y._index(i, i)]);
                    }
                }
            } else {
                if (i < o.cols) {
                    if (comptime meta.diagOf(X) == .unit) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                o.data[o._index(i, i)] = numeric.two(meta.Numeric(O))
                            else
                                o.data[o._index(i, i)] = numeric.zero(meta.Numeric(O));
                        } else {
                            op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y.data[y._index(i, i)]);
                        }
                    } else {
                        if (comptime meta.diagOf(Y) == .unit)
                            op_(&o.data[o._index(i, i)], x.data[x._index(i, i)], numeric.one(meta.Numeric(Y)))
                        else
                            op_(&o.data[o._index(i, i)], x.data[x._index(i, i)], y.data[y._index(i, i)]);
                    }
                }

                var j: usize = int.min(i + 1, o.cols);
                while (j < o.cols) : (j += 1) {
                    op_(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                }
            }
        }
    }
}
