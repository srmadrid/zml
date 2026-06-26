const meta = @import("../../meta.zig");

pub fn matmulIntoUnchecked(o: anytype, x: anytype, y: anytype) void {
    const O = meta.Child(@TypeOf(o));
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    switch (comptime O.direction) {
        .forward => switch (comptime X.direction) {
            .forward => switch (comptime Y.direction) {
                .forward => {
                    inline for (0..O.rows) |i| {
                        o.idx[i] = y.idx[x.idx[i]];
                    }
                },
                .backward => {
                    inline for (0..O.rows) |j| {
                        o.idx[y.idx[j]] = j;
                    }

                    const MSB: usize = @as(usize, 1) << (@bitSizeOf(usize) - 1);

                    inline for (0..O.rows) |i| {
                        if ((o.idx[i] & MSB) == 0) {
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
                    }

                    inline for (0..O.rows) |i| {
                        o.idx[i] &= ~MSB;
                    }
                },
            },
            .backward => switch (comptime Y.direction) {
                .forward => {
                    inline for (0..O.rows) |i| {
                        o.idx[x.idx[i]] = y.idx[i];
                    }
                },
                .backward => {
                    inline for (0..O.rows) |j| {
                        o.idx[x.idx[y.idx[j]]] = j;
                    }
                },
            },
        },
        .backward => switch (comptime X.direction) {
            .forward => switch (comptime Y.direction) {
                .forward => {
                    inline for (0..O.rows) |i| {
                        o.idx[y.idx[x.idx[i]]] = i;
                    }
                },
                .backward => {
                    inline for (0..O.rows) |i| {
                        o.idx[x.idx[i]] = i;
                    }

                    const MSB: usize = @as(usize, 1) << (@bitSizeOf(usize) - 1);

                    inline for (0..O.rows) |j| {
                        if ((o.idx[j] & MSB) == 0) {
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
                    }

                    inline for (0..O.rows) |j| {
                        o.idx[j] &= ~MSB;
                    }
                },
            },
            .backward => switch (comptime Y.direction) {
                .forward => {
                    inline for (0..O.rows) |i| {
                        o.idx[y.idx[i]] = x.idx[i];
                    }
                },
                .backward => {
                    inline for (0..O.rows) |j| {
                        o.idx[j] = x.idx[y.idx[j]];
                    }
                },
            },
        },
    }
}
