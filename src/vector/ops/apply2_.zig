const meta = @import("../../meta.zig");

const vector = @import("../../vector.zig");

/// Applies a binary in-place operation elementwise between an output and two
/// input vectors, or between an output vector, an input vector and an input
/// numeric.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
///
/// For two static vectors the function cannot return an error unless op can.
///
/// For two input sparse vectors, or an input sparse vector and an input
/// numeric, the operation is only applied to the indices where at least one of
/// the vectors has a non-zero element.
///
/// ## Signature
/// ```zig
/// vector.apply2_(*O, x: X, y: Y, op_: Op) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
/// * `op_` (`comptime anytype`): An in-place binary numeric function to apply
///   elementwise to `o`, `x` and `y`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `vector.Error.DimensionMismatch`: If the vectors do not have the same
///   length.
pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Op: type = @TypeOf(op_);
    const opinfo = @typeInfo(Op);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        (!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.vector.apply2_: o must be a mutable one-itme pointer to a vector, at least one of x or y must be a vector, the other must be a vector or a numeric, and op_ must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

    O = meta.Child(O);

    const o_len = if (comptime meta.isStaticVector(O)) O.len else o.len;
    const x_len = if (comptime meta.isVector(X)) (if (comptime meta.isStaticVector(X)) X.len else x.len) else o_len;
    const y_len = if (comptime meta.isVector(Y)) (if (comptime meta.isStaticVector(Y)) Y.len else y.len) else o_len;

    if (comptime meta.isStaticVector(O) and
        (meta.isStaticVector(X) or meta.isNumeric(X)) and
        (meta.isStaticVector(Y) or meta.isNumeric(Y)))
    {
        if (comptime o_len != x_len or o_len != y_len)
            if (comptime meta.isVector(X) and meta.isVector(Y))
                @compileError("zsl.vector.apply2: static vectors o, x and y must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n")
            else if (comptime meta.isVector(X))
                @compileError("zsl.vector.apply2: static vectors o and x must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n")
            else
                @compileError("zsl.vector.apply2: static vectors o and y must have equal compile time lengths, got\n\to: *" ++
                    @typeName(O) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");
    } else {
        if (o_len != x_len or o_len != y_len)
            return vector.Error.DimensionMismatch;
    }

    switch (comptime meta.vectorType(O)) {
        .static => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/ststst.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/ststde.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/ststsp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/ststnu.zig").apply2_(o, x, y, op_),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/stdest.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/stdede.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/stdesp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/stdenu.zig").apply2_(o, x, y, op_),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/stspst.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/stspde.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/stspsp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/stspnu.zig").apply2_(o, x, y, op_),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/stnust.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/stnude.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/stnusp.zig").apply2_(o, x, y, op_),
                .numeric => unreachable,
            },
        },
        .dense => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/destst.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/destde.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/destsp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/destnu.zig").apply2_(o, x, y, op_),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/dedest.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/dedede.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/dedesp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/dedenu.zig").apply2_(o, x, y, op_),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/despst.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/despde.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/despsp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/despnu.zig").apply2_(o, x, y, op_),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2_/denust.zig").apply2_(o, x, y, op_),
                .dense => return @import("apply2_/denude.zig").apply2_(o, x, y, op_),
                .sparse => return @import("apply2_/denusp.zig").apply2_(o, x, y, op_),
                .numeric => unreachable,
            },
        },
        .sparse => switch (comptime meta.vectorType(X)) {
            .sparse => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2_/spspsp.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/spspnu.zig").apply2_(o, x, y, op_),
                else => @compileError("zsl.vector.apply2_: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2_/spnusp.zig").apply2_(o, x, y, op_),
                .numeric => unreachable,
                else => @compileError("zsl.vector.apply2_: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            else => @compileError("zsl.vector.apply2_: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
        },
        .numeric => unreachable,
    }
}
