const meta = @import("../../meta.zig");

const vector = @import("../../vector.zig");

/// Applies a binary into operation elementwise between an output and two input
/// vectors, or between an output vector, an input vector and an input numeric.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
///
/// For three static vectors, or a static output vector, an input static vector
/// and an input numeric, the function cannot return an error unless opInto can.
///
/// For two input sparse vectors, or an input sparse vector and an input
/// numeric, the operation is only applied to the indices where at least one of
/// the vectors has a non-zero element.
///
/// ## Signature
/// ```zig
/// vector.apply2Into(*O, x: X, y: Y, opInto: OpInto) !void
/// ```
///
/// ## Arguments
/// * `o` (`anytype`): The output operand.
/// * `x` (`anytype`): The left input operand.
/// * `y` (`anytype`): The right input operand.
/// * `opInto` (`comptime anytype`): An ito binary numeric function to apply
///   elementwise to `o`, `x` and `y`.
///
/// ## Returns
/// `void`
///
/// ## Errors
/// * `vector.Error.DimensionMismatch`: If the vectors do not have the same
///   length.
pub fn apply2Into(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isVector(meta.Child(O)) or
        (!meta.isVector(X) and !meta.isNumeric(X)) or (!meta.isVector(Y) and !meta.isNumeric(Y)) or
        (!meta.isVector(X) and !meta.isVector(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.vector.apply2Into: o must be a mutable one-itme pointer to a vector, at least one of x or y must be a vector, the other must be a vector or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

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
                .static => return @import("apply2Into/ststst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/ststde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/ststsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/ststnu.zig").apply2Into(o, x, y, opInto),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/stdest.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/stdede.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/stdesp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/stdenu.zig").apply2Into(o, x, y, opInto),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/stspst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/stspde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/stspsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/stspnu.zig").apply2Into(o, x, y, opInto),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/stnust.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/stnude.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/stnusp.zig").apply2Into(o, x, y, opInto),
                .numeric => unreachable,
            },
        },
        .dense => switch (comptime meta.vectorType(X)) {
            .static => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/destst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/destde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/destsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/destnu.zig").apply2Into(o, x, y, opInto),
            },
            .dense => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/dedest.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/dedede.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/dedesp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/dedenu.zig").apply2Into(o, x, y, opInto),
            },
            .sparse => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/despst.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/despde.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/despsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/despnu.zig").apply2Into(o, x, y, opInto),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .static => return @import("apply2Into/denust.zig").apply2Into(o, x, y, opInto),
                .dense => return @import("apply2Into/denude.zig").apply2Into(o, x, y, opInto),
                .sparse => return @import("apply2Into/denusp.zig").apply2Into(o, x, y, opInto),
                .numeric => unreachable,
            },
        },
        .sparse => switch (comptime meta.vectorType(X)) {
            .sparse => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2Into/spspsp.zig").apply2Into(o, x, y, opInto),
                .numeric => return @import("apply2Into/spspnu.zig").apply2Into(o, x, y, opInto),
                else => @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
            },
            .numeric => switch (comptime meta.vectorType(Y)) {
                .sparse => return @import("apply2Into/spnusp.zig").apply2Into(o, x, y, opInto),
                .numeric => unreachable,
                else => @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                    @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
            },
            else => @compileError("zsl.vector.apply2Into: o cannot point to a sparse vector if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n"),
        },
        .numeric => unreachable,
    }
}
