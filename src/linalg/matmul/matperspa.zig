const meta = @import("../../meta.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    switch (comptime O.direction) {
        .forward => switch (comptime X.direction) {
            .forward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.idx[i] = y.idx[x.idx[i]];
                    }
                },
                .backward => {
                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        o.idx[y.idx[j]] = j;
                    }

                    const MSB: usize = @as(usize, 1) << (@bitSizeOf(usize) - 1);

                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        if ((o.idx[i] & MSB) != 0)
                            continue;

                        const temp = o.idx[i];
                        var curr = i;
                        var next = x.idx[curr];

                        while (next != i) {
                            o.idx[curr] = o.idx[next] | MSB;
                            curr = next;
                            next = x.idx[curr];
                        }

                        o.idx[curr] = temp | MSB;
                    }

                    i = 0;
                    while (i < o.rows) : (i += 1) {
                        o.idx[i] &= ~MSB;
                    }
                },
            },
            .backward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.idx[x.idx[i]] = y.idx[i];
                    }
                },
                .backward => {
                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        o.idx[x.idx[y.idx[j]]] = j;
                    }
                },
            },
        },
        .backward => switch (comptime X.direction) {
            .forward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.idx[y.idx[x.idx[i]]] = i;
                    }
                },
                .backward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.idx[x.idx[i]] = i;
                    }

                    const MSB: usize = @as(usize, 1) << (@bitSizeOf(usize) - 1);

                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        if ((o.idx[j] & MSB) != 0)
                            continue;

                        const temp = o.idx[j];
                        var curr = j;
                        var next = y.idx[curr];

                        while (next != j) {
                            o.idx[curr] = o.idx[next] | MSB;
                            curr = next;
                            next = y.idx[curr];
                        }

                        o.idx[curr] = temp | MSB;
                    }

                    j = 0;
                    while (j < o.rows) : (j += 1) {
                        o.idx[j] &= ~MSB;
                    }
                },
            },
            .backward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.idx[y.idx[i]] = x.idx[i];
                    }
                },
                .backward => {
                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        o.idx[j] = x.idx[y.idx[j]];
                    }
                },
            },
        },
    }
}
