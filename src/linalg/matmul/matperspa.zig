const meta = @import("../../meta.zig");

pub fn matmulInto(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    switch (comptime O.direction) {
        .forward => switch (comptime X.direction) {
            .forward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.data[i] = y.data[x.data[i]];
                    }
                },
                .backward => {
                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        o.data[y.data[j]] = j;
                    }

                    const MSB: usize = @as(usize, 1) << (@bitSizeOf(usize) - 1);

                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        if ((o.data[i] & MSB) != 0)
                            continue;

                        const temp = o.data[i];
                        var curr = i;
                        var next = x.data[curr];

                        while (next != i) {
                            o.data[curr] = o.data[next] | MSB;
                            curr = next;
                            next = x.data[curr];
                        }

                        o.data[curr] = temp | MSB;
                    }

                    i = 0;
                    while (i < o.rows) : (i += 1) {
                        o.data[i] &= ~MSB;
                    }
                },
            },
            .backward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.data[x.data[i]] = y.data[i];
                    }
                },
                .backward => {
                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        o.data[x.data[y.data[j]]] = j;
                    }
                },
            },
        },
        .backward => switch (comptime X.direction) {
            .forward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.data[y.data[x.data[i]]] = i;
                    }
                },
                .backward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.data[x.data[i]] = i;
                    }

                    const MSB: usize = @as(usize, 1) << (@bitSizeOf(usize) - 1);

                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        if ((o.data[j] & MSB) != 0)
                            continue;

                        const temp = o.data[j];
                        var curr = j;
                        var next = y.data[curr];

                        while (next != j) {
                            o.data[curr] = o.data[next] | MSB;
                            curr = next;
                            next = y.data[curr];
                        }

                        o.data[curr] = temp | MSB;
                    }

                    j = 0;
                    while (j < o.rows) : (j += 1) {
                        o.data[j] &= ~MSB;
                    }
                },
            },
            .backward => switch (comptime Y.direction) {
                .forward => {
                    var i: usize = 0;
                    while (i < o.rows) : (i += 1) {
                        o.data[y.data[i]] = x.data[i];
                    }
                },
                .backward => {
                    var j: usize = 0;
                    while (j < o.rows) : (j += 1) {
                        o.data[j] = x.data[y.data[j]];
                    }
                },
            },
        },
    }
}
