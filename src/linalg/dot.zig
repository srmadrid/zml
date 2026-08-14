const linalg = @import("../linalg.zig");
const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

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

/// Computes the dot product of two vectors, `x ⋅ y = xᴴ y = ∑ᵢ x̄ᵢ yᵢ`.
///
/// ## Signature
/// ```zig
/// linalg.dot(x: X, y: Y) !linalg.Dot(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector.
/// * `y` (`anytype`): The right vector.
///
/// ## Returns
/// `linalg.Dot(@TypeOf(x), @TypeOf(y))`: The dot product `x ⋅ y`.
///
/// ## Errors
/// * `linalg.Error.DimensionMismatch`: If the inputs do not have equal length.
pub fn dot(x: anytype, y: anytype) !linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    const x_len = if (comptime meta.isStaticVector(X)) X.len else x.len;
    const y_len = if (comptime meta.isStaticVector(Y)) Y.len else y.len;

    if (comptime !meta.isStaticVector(X) or !meta.isStaticVector(Y)) {
        if (x_len != y_len)
            return linalg.Error.DimensionMismatch;
    }

    return dotUnchecked(x, y);
}

/// Computes the dot product of two vectors, `x ⋅ y = xᴴ y = ∑ᵢ x̄ᵢ yᵢ`, without
/// performing dimension checks.
///
/// ## Signature
/// ```zig
/// linalg.dotUnchecked(x: X, y: Y) linalg.Dot(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector.
/// * `y` (`anytype`): The right vector.
///
/// ## Returns
/// `linalg.Dot(@TypeOf(x), @TypeOf(y))`: The dot product `x ⋅ y`.
pub fn dotUnchecked(x: anytype, y: anytype) linalg.Dot(@TypeOf(x), @TypeOf(y)) {
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    return switch (comptime meta.vectorType(X)) {
        .static => switch (comptime meta.vectorType(Y)) {
            .static => @import("dot/vecsta_vecsta.zig").dotUnchecked(x, y),
            .dense => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecsta_vecden.zig").dotUnchecked(x, y),
            .sparse => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecsta_vecspa.zig").dotUnchecked(x, y),
            .numeric => unreachable,
        },
        .dense => switch (comptime meta.vectorType(Y)) {
            .static => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecden_vecsta.zig").dotUnchecked(x, y),
            .dense => @import("dot/vecden_vecden.zig").dotUnchecked(x, y),
            .sparse => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecden_vecspa.zig").dotUnchecked(x, y),
            .numeric => unreachable,
        },
        .sparse => switch (comptime meta.vectorType(Y)) {
            .static => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecspa_vecsta.zig").dotUnchecked(x, y),
            .dense => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecspa_vecden.zig").dotUnchecked(x, y),
            .sparse => @import("dot/slow.zig").dotUnchecked(x, y), // @import("dot/vecspa_vecspa.zig").dotUnchecked(x, y),
            .numeric => unreachable,
        },
        .numeric => unreachable,
    };
}
