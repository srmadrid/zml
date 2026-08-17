const std = @import("std");

const meta = @import("../../meta.zig");

const numeric = @import("../../numeric.zig");
const array = @import("../../array.zig");

const arrutils = @import("../utils.zig");

pub fn apply2IntoUnchecked(o: anytype, x: anytype, y: anytype, comptime opInto: anytype) void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const OpInto: type = @TypeOf(opInto);
    const opinfo = @typeInfo(OpInto);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isArray(meta.Child(O)) or
        (!meta.isArray(X) and !meta.isNumeric(X)) or (!meta.isArray(Y) and !meta.isNumeric(Y)) or
        (!meta.isArray(X) and !meta.isArray(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.array.apply2IntoUnchecked: o must be a mutable one-item pointer to an array, at least one of x and y must be an array, the other must be an array or a numeric, and opInto must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\topInto: " ++ @typeName(OpInto) ++ "\n");

    comptime if (meta.isBuilderSparseArray(X) or meta.isBuilderSparseArray(Y))
        @compileError("zsl.array.apply2IntoUnchecked: builder array types are not allowed as inputs, got\n\tx: " ++
            @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n");

    O = meta.Child(O);

    switch (comptime meta.arrayType(O)) {
        .static => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
        .dense => switch (comptime meta.arrayType(X)) {
            .static => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
            .dense => switch (comptime meta.arrayType(Y)) {
                .static => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
                .dense => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
                .sparse => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
                    @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
                .numeric => return @import("apply2/arrden_arrden_num.zig").apply2IntoUnchecked(o, x, y, opInto),
                .builder_sparse => unreachable,
            },
            .sparse => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
            .numeric => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
            else => unreachable,
        },
        .sparse => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
        .builder_sparse => @compileError("zsl.array.apply2IntoUnchecked: not implemented yet for\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n"),
        .numeric => unreachable,
    }
}
