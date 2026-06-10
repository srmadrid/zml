const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const linalg = @import("../linalg.zig");

pub fn Dot(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isVector(X) or !meta.isVector(Y))
        @compileError("zsl.linalg.Dot: X and Y must be vector types, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    comptime if (meta.isStaticVector(X) and meta.isStaticVector(Y)) {
        if (X.len != Y.len)
            @compileError("zsl.linalg.Dot: static vector types X and Y must have equal lengths, got\n\tX = " ++
                @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");
    };

    return numeric.Add(
        numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
        numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
    );
}

pub fn dot(x: anytype, y: anytype) !linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;
    const y_len = if (comptime meta.isStaticVector(Y)) Y.len else y.len;

    if (comptime !meta.isStaticVector(X) or !meta.isStaticVector(Y)) {
        if (x_len != y_len)
            return linalg.Error.DimensionMismatch;
    }

    switch (comptime meta.vectorType(X)) {
        .static => switch (comptime meta.vectorType(Y)) {
            .static => return @import("dot/stst.zig").dot(x, y),
            .numeric => unreachable,
            else => @compileError("zsl.linalg.dot: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n"),
        },
        .numeric => unreachable,
        else => @compileError("zsl.linalg.dot: not implemented yet for \n\tX = " ++ @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n"),
    }
}
