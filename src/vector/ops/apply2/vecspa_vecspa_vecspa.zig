const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    return switch (o.flags.noconj) {
        true => switch (x.flags.noconj) {
            true => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, true, true, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, true, true, false),
            },
            false => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, true, false, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, true, false, false),
            },
        },
        false => switch (x.flags.noconj) {
            true => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, false, true, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, false, true, false),
            },
            false => switch (y.flags.noconj) {
                true => k_apply2IntoUnchecked(o, x, y, opInto, false, false, true),
                false => k_apply2IntoUnchecked(o, x, y, opInto, false, false, false),
            },
        },
    };
}

fn k_apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype, comptime noconj_o: bool, comptime noconj_x: bool, comptime noconj_y: bool) void {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    var ix: usize = 0;
    var iy: usize = 0;
    var io: usize = 0;
    while (ix < x.nnz and iy < y.nnz) {
        if (x.idx[ix] == y.idx[iy]) {
            opInto(
                &o.data[io],
                if (comptime noconj_o == noconj_x)
                    x.data[ix]
                else
                    numeric.conj(x.data[ix]),
                if (comptime noconj_o == noconj_y)
                    y.data[iy]
                else
                    numeric.conj(y.data[iy]),
            );

            o.idx[io] = x.idx[ix];
            ix += 1;
            iy += 1;
            io += 1;
        } else if (x.idx[ix] < y.idx[iy]) {
            opInto(
                &o.data[io],
                if (comptime noconj_o == noconj_x)
                    x.data[ix]
                else
                    numeric.conj(x.data[ix]),
                numeric.cast(meta.Numeric(Y), 0),
            );

            o.idx[io] = x.idx[ix];
            ix += 1;
            io += 1;
        } else {
            opInto(
                &o.data[io],
                numeric.cast(meta.Numeric(X), 0),
                if (comptime noconj_o == noconj_y)
                    y.data[iy]
                else
                    numeric.conj(y.data[iy]),
            );

            o.idx[io] = y.idx[iy];
            iy += 1;
            io += 1;
        }
    }

    while (ix < x.nnz) {
        opInto(
            &o.data[io],
            if (comptime noconj_o == noconj_x)
                x.data[ix]
            else
                numeric.conj(x.data[ix]),
            numeric.cast(meta.Numeric(Y), 0),
        );

        o.idx[io] = x.idx[ix];
        ix += 1;
        io += 1;
    }

    while (iy < y.nnz) {
        opInto(
            &o.data[io],
            numeric.cast(meta.Numeric(X), 0),
            if (comptime noconj_o == noconj_y)
                y.data[iy]
            else
                numeric.conj(y.data[iy]),
        );

        o.idx[io] = y.idx[iy];
        iy += 1;
        io += 1;
    }

    o.nnz = io;
}
