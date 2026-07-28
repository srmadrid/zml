const linalg = @import("../../linalg.zig");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

pub fn norm(x: anytype, comptime order: linalg.NormOrder(meta.Numeric(@TypeOf(x)))) linalg.Norm(@TypeOf(x), order) {
    const X = @TypeOf(x);
    const R = linalg.Norm(X, order);

    switch (comptime order) {
        .l1, .l2, .inf, .p => return numeric.one(R),
        .frobenius => return numeric.cast(R, numeric.sqrt(numeric.cast(R, x.rows))),
    }
}
