const meta = @import("../../../meta.zig");

pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    const O: type = meta.Child(@TypeOf(o));
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    if (comptime meta.layoutOf(O) == .col_major) {
        inline for (0..O.cols) |j| {
            if (comptime meta.uploOf(O) == .upper) {
                inline for (0..j) |i| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }

                opInto(
                    &o.data[O._index(j, j)],
                    if (comptime meta.isMatrix(X))
                        x.get(j, j) catch unreachable
                    else
                        x,
                    if (comptime meta.isMatrix(Y))
                        y.get(j, j) catch unreachable
                    else
                        y,
                );
            } else {
                opInto(
                    &o.data[O._index(j, j)],
                    if (comptime meta.isMatrix(X))
                        x.get(j, j) catch unreachable
                    else
                        x,
                    if (comptime meta.isMatrix(Y))
                        y.get(j, j) catch unreachable
                    else
                        y,
                );

                inline for (j + 1..O.rows) |i| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }
            }
        }
    } else {
        inline for (0..O.rows) |i| {
            if (comptime meta.uploOf(O) == .lower) {
                inline for (0..i) |j| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }

                opInto(
                    &o.data[O._index(i, i)],
                    if (comptime meta.isMatrix(X))
                        x.get(i, i) catch unreachable
                    else
                        x,
                    if (comptime meta.isMatrix(Y))
                        y.get(i, i) catch unreachable
                    else
                        y,
                );
            } else {
                opInto(
                    &o.data[O._index(i, i)],
                    if (comptime meta.isMatrix(X))
                        x.get(i, i) catch unreachable
                    else
                        x,
                    if (comptime meta.isMatrix(Y))
                        y.get(i, i) catch unreachable
                    else
                        y,
                );

                inline for (i + 1..O.cols) |j| {
                    opInto(
                        &o.data[O._index(i, j)],
                        if (comptime meta.isMatrix(X))
                            x.get(i, j) catch unreachable
                        else
                            x,
                        if (comptime meta.isMatrix(Y))
                            y.get(i, j) catch unreachable
                        else
                            y,
                    );
                }
            }
        }
    }
}
