const int = @import("../../../int.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        var j: usize = 0;
        while (j < o.cols) : (j += 1) {
            var i: usize = 0;
            if (comptime meta.uploOf(X) == .upper) {
                if (comptime meta.uploOf(Y) == .upper) {
                    while (i < int.min(j, o.rows)) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                    }
                } else {
                    while (i < int.min(j, o.rows)) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], numeric.zero(meta.Numeric(Y)));
                    }
                }
            } else {
                if (comptime meta.uploOf(Y) == .upper) {
                    while (i < int.min(j, o.rows)) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(i, j)]);
                    }
                } else {
                    while (i < int.min(j, o.rows)) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), numeric.zero(meta.Numeric(Y)));
                    }
                }
            }

            if (j < o.rows) {
                if (comptime meta.diagOf(X) == .unit) {
                    if (comptime meta.diagOf(Y) == .unit)
                        opInto(&o.data[o._index(j, j)], numeric.one(meta.Numeric(X)), numeric.one(meta.Numeric(Y)))
                    else
                        opInto(&o.data[o._index(j, j)], numeric.one(meta.Numeric(X)), y.data[y._index(j, j)]);
                } else {
                    if (comptime meta.diagOf(Y) == .unit)
                        opInto(&o.data[o._index(j, j)], x.data[x._index(j, j)], numeric.one(meta.Numeric(Y)))
                    else
                        opInto(&o.data[o._index(j, j)], x.data[x._index(j, j)], y.data[y._index(j, j)]);
                }
            }

            i = int.min(j + 1, o.rows);
            if (comptime meta.uploOf(X) == .lower) {
                if (comptime meta.uploOf(Y) == .lower) {
                    while (i < o.rows) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                    }
                } else {
                    while (i < o.rows) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], numeric.zero(meta.Numeric(Y)));
                    }
                }
            } else {
                if (comptime meta.uploOf(Y) == .lower) {
                    while (i < o.rows) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(i, j)]);
                    }
                } else {
                    while (i < o.rows) : (i += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), numeric.zero(meta.Numeric(Y)));
                    }
                }
            }
        }
    } else {
        var i: usize = 0;
        while (i < o.rows) : (i += 1) {
            var j: usize = 0;
            if (comptime meta.uploOf(X) == .lower) {
                if (comptime meta.uploOf(Y) == .lower) {
                    while (j < int.min(i, o.cols)) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                    }
                } else {
                    while (j < int.min(i, o.cols)) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], numeric.zero(meta.Numeric(Y)));
                    }
                }
            } else {
                if (comptime meta.uploOf(Y) == .lower) {
                    while (j < int.min(i, o.cols)) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(i, j)]);
                    }
                } else {
                    while (j < int.min(i, o.cols)) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), numeric.zero(meta.Numeric(Y)));
                    }
                }
            }

            if (i < o.cols) {
                if (comptime meta.diagOf(X) == .unit) {
                    if (comptime meta.diagOf(Y) == .unit)
                        opInto(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), numeric.one(meta.Numeric(Y)))
                    else
                        opInto(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y.data[y._index(i, i)]);
                } else {
                    if (comptime meta.diagOf(Y) == .unit)
                        opInto(&o.data[o._index(i, i)], x.data[x._index(i, i)], numeric.one(meta.Numeric(Y)))
                    else
                        opInto(&o.data[o._index(i, i)], x.data[x._index(i, i)], y.data[y._index(i, i)]);
                }
            }

            j = int.min(i + 1, o.cols);
            if (comptime meta.uploOf(X) == .upper) {
                if (comptime meta.uploOf(Y) == .upper) {
                    while (j < o.cols) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], y.data[y._index(i, j)]);
                    }
                } else {
                    while (j < o.cols) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], x.data[x._index(i, j)], numeric.zero(meta.Numeric(Y)));
                    }
                }
            } else {
                if (comptime meta.uploOf(Y) == .upper) {
                    while (j < o.cols) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), y.data[y._index(i, j)]);
                    }
                } else {
                    while (j < o.cols) : (j += 1) {
                        opInto(&o.data[o._index(i, j)], numeric.zero(meta.Numeric(X)), numeric.zero(meta.Numeric(Y)));
                    }
                }
            }
        }
    }
}
