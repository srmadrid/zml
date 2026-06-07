const std = @import("std");

const meta = @import("../../../meta.zig");
const int = @import("../../../int.zig");
const numeric = @import("../../../numeric.zig");
const matrix = @import("../../../matrix.zig");

pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const aliased = (comptime O == Y) and std.meta.eql(o.*, y);

    if ((comptime op_ == numeric.subInto) or !aliased) {
        if (comptime meta.layoutOf(O) == .col_major) {
            var j: usize = 0;
            while (j < o.cols) : (j += 1) {
                if (comptime meta.uploOf(O) == .upper) {
                    var i: usize = 0;
                    while (i < int.min(j, o.rows)) : (i += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }

                    if (j < o.rows) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(j, j)], numeric.one(meta.Numeric(O)))
                            else
                                numeric.set(&o.data[o._index(j, j)], numeric.neg(numeric.one(meta.Numeric(O))));
                        } else {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(j, j)], y.data[y._index(j, j)])
                            else
                                numeric.set(&o.data[o._index(j, j)], numeric.neg(y.data[y._index(j, j)]));
                        }
                    }
                } else {
                    if (j < o.rows) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(j, j)], numeric.one(meta.Numeric(O)))
                            else
                                numeric.set(&o.data[o._index(j, j)], numeric.neg(numeric.one(meta.Numeric(O))));
                        } else {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(j, j)], y.data[y._index(j, j)])
                            else
                                numeric.set(&o.data[o._index(j, j)], numeric.neg(y.data[y._index(j, j)]));
                        }
                    }

                    var i: usize = int.min(j + 1, o.rows);
                    while (i < o.rows) : (i += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }
                }
            }
        } else {
            var i: usize = 0;
            while (i < o.rows) : (i += 1) {
                if (comptime meta.uploOf(O) == .lower) {
                    var j: usize = 0;
                    while (j < int.min(i, o.cols)) : (j += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }

                    if (i < o.cols) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(i, i)], numeric.one(meta.Numeric(O)))
                            else
                                numeric.set(&o.data[o._index(i, i)], numeric.neg(numeric.one(meta.Numeric(O))));
                        } else {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(i, i)], y.data[y._index(i, i)])
                            else
                                numeric.set(&o.data[o._index(i, i)], numeric.neg(y.data[y._index(i, i)]));
                        }
                    }
                } else {
                    if (i < o.cols) {
                        if (comptime meta.diagOf(Y) == .unit) {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(i, i)], numeric.one(meta.Numeric(O)))
                            else
                                numeric.set(&o.data[o._index(i, i)], numeric.neg(numeric.one(meta.Numeric(O))));
                        } else {
                            if (comptime op_ == numeric.addInto)
                                numeric.set(&o.data[o._index(i, i)], y.data[y._index(i, i)])
                            else
                                numeric.set(&o.data[o._index(i, i)], numeric.neg(y.data[y._index(i, i)]));
                        }
                    }

                    var j: usize = int.min(i + 1, o.cols);
                    while (j < o.cols) : (j += 1) {
                        if (comptime op_ == numeric.addInto)
                            numeric.set(&o.data[o._index(i, j)], y.data[y._index(i, j)])
                        else
                            numeric.set(&o.data[o._index(i, j)], numeric.neg(y.data[y._index(i, j)]));
                    }
                }
            }
        }
    }

    if (comptime meta.diagOf(X) == .unit) {
        var i: usize = 0;
        while (i < int.min(o.rows, o.cols)) : (i += 1) {
            if (comptime meta.diagOf(Y) == .unit) {
                op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), numeric.one(meta.Numeric(Y)));
            } else {
                if ((comptime op_ == numeric.addInto) or !aliased)
                    op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), y.data[y._index(i, i)])
                else
                    op_(&o.data[o._index(i, i)], numeric.one(meta.Numeric(X)), numeric.neg(y.data[y._index(i, i)]));
            }
        }
    }

    if (comptime meta.layoutOf(X) == .col_major) {
        var j: usize = 0;
        while (j < x.cols) : (j += 1) {
            var p: usize = x.ptr[j];
            while (p < x.ptr[j + 1]) : (p += 1) {
                const ty = if (x.idx[p] == j)
                    (if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(j, j)])
                else
                    y.data[y._index(x.idx[p], j)];

                if ((comptime op_ == numeric.addInto) or !aliased)
                    op_(&o.data[o._index(x.idx[p], j)], x.data[p], ty)
                else
                    op_(&o.data[o._index(x.idx[p], j)], x.data[p], numeric.neg(ty));
            }
        }
    } else {
        var i: usize = 0;
        while (i < x.rows) : (i += 1) {
            var p: usize = x.ptr[i];
            while (p < x.ptr[i + 1]) : (p += 1) {
                const ty = if (i == x.idx[p])
                    (if (comptime meta.diagOf(Y) == .unit) numeric.one(meta.Numeric(Y)) else y.data[y._index(i, i)])
                else
                    y.data[y._index(i, x.idx[p])];

                if ((comptime op_ == numeric.addInto) or !aliased)
                    op_(&o.data[o._index(i, x.idx[p])], x.data[p], ty)
                else
                    op_(&o.data[o._index(i, x.idx[p])], x.data[p], numeric.neg(ty));
            }
        }
    }
}
