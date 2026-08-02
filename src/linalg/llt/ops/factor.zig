const std = @import("std");

const linalg = @import("../../../linalg.zig");
const meta = @import("../../../meta.zig");

pub fn Factor(comptime X: type) type {
    comptime if (!(meta.isSymmetricMatrix(X) and meta.isReal(meta.Numeric(X))) and
        !(meta.isHermitianMatrix(X) and meta.isComplex(meta.Numeric(X))))
        @compileError("zsl.linalg.llt.Factor: X must be a real symmetric or a complex Hermitian matrix type, got \n\tX = " ++ @typeName(X) ++ "\n");

    return switch (comptime meta.matrixStorage(X)) {
        .static => linalg.llt.Static(X.size, meta.Numeric(X), meta.layoutOf(X).?),
        .dense => linalg.llt.Dense(meta.Numeric(X), meta.layoutOf(X).?),
        .sparse => @compileError("zsl.linalg.llt.Factor: not implemented yet for sparse matrices"),
        .numeric => unreachable,
    };
}

pub fn factor(x: anytype) !linalg.llt.Factor(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = linalg.llt.Factor(X);

    if (comptime !@hasDecl(meta.Child(R), "is_static") or !meta.Child(R).is_static)
        @compileError("zsl.linalg.llt.factor: the result cannot be a heap-allocated llt type, i.e., x must be a static matrix, got\n\tx: " ++
            @typeName(X) ++ "\n\tresult: " ++ @typeName(R) ++ "\nFor this input use zsl.linalg.llt.factorAlloc instead.");

    var result = R.init;

    try factorIntoUnchecked(&result, x);

    return result;
}

pub fn factorAlloc(allocator: std.mem.Allocator, x: anytype) !linalg.llt.Factor(@TypeOf(x)) {
    const X: type = @TypeOf(x);
    const R: type = linalg.llt.Factor(X);

    const r_storage: meta.MatrixStorage = if (comptime @hasDecl(meta.Child(R), "is_static") and meta.Child(R).is_static)
        .static
    else if (comptime @hasDecl(meta.Child(R), "is_dense") and meta.Child(R).is_dense)
        .dense
    else if (comptime @hasDecl(meta.Child(R), "is_sparse") and meta.Child(R).is_sparse)
        .sparse
    else
        unreachable;

    var result = switch (comptime r_storage) {
        .static => R.init,
        .dense => try R.init(allocator, x.rows),
        // .sparse =>
        .numeric => unreachable,
        else => @compileError("zsl.linalg.llt.factorAlloc: not implemented yet for x: " ++ @typeName(X) ++ "\n"),
    };

    try factorIntoUnchecked(&result, x);

    return result;
}

pub fn factorInto(o: anytype, x: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or (!@hasDecl(meta.Child(O), "is_llt") or !meta.Child(O).is_llt) or
        (!(meta.isSymmetricMatrix(X) and meta.isReal(meta.Numeric(X))) and !(meta.isHermitianMatrix(X) and meta.isComplex(meta.Numeric(X)))))
        @compileError("zsl.linalg.llt.factorInto: o must be a mutable one-item pointer to a llt object, and x must be a real symmetric or a complex Hermitian matrix, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    const o_rows = if (comptime @hasDecl(meta.Child(O), "is_static") and meta.Child(O).is_static) O.rows else o.rows;
    const o_cols = if (comptime @hasDecl(meta.Child(O), "is_static") and meta.Child(O).is_static) O.cols else o.cols;
    const x_rows = if (comptime meta.isStaticMatrix(X)) X.rows else x.rows;
    const x_cols = if (comptime meta.isStaticMatrix(X)) X.cols else x.cols;

    if (comptime @hasDecl(meta.Child(O), "is_static") and meta.Child(O).is_static and meta.isStaticMatrix(X)) {
        if (comptime o_rows != x_rows or o_cols != x_cols)
            @compileError("zsl.linalg.llt.factorInto: static llt object o and static matrix xmust have equal compile time dimensions, got\n\to: *" ++
                @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");
    }

    if (comptime @hasDecl(meta.Child(O), "is_sparse") and meta.Child(O).is_sparse) {
        // nnz checks
    }

    return factorIntoUnchecked(o, x);
}

pub fn factorIntoUnchecked(o: anytype, x: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or (!@hasDecl(meta.Child(O), "is_llt") or !meta.Child(O).is_llt) or
        (!(meta.isSymmetricMatrix(X) and meta.isReal(meta.Numeric(X))) and !(meta.isHermitianMatrix(X) and meta.isComplex(meta.Numeric(X)))))
        @compileError("zsl.linalg.llt.factorIntoUnchecked: o must be a mutable one-item pointer to a llt object, and x must be a real symmetric or a complex Hermitian matrix, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n");

    O = meta.Child(O);

    const o_storage: meta.MatrixStorage = if (comptime @hasDecl(meta.Child(O), "is_static") and meta.Child(O).is_static)
        .static
    else if (comptime @hasDecl(meta.Child(O), "is_dense") and meta.Child(O).is_dense)
        .dense
    else if (comptime @hasDecl(meta.Child(O), "is_sparse") and meta.Child(O).is_sparse)
        .sparse
    else
        unreachable;

    switch (comptime o_storage) {
        .static => switch (comptime meta.matrixStorage(X)) {
            .static => try @import("factor/lltsta_mat___sta.zig").factorIntoUnchecked(o, x),
            .dense => try @import("factor/lltsta_mat___den.zig").factorIntoUnchecked(o, x),
            // .sparse => try @import("factor/lltsta_mat___spa.zig").factorIntoUnchecked(o, x),
            .numeric => unreachable,
            else => @compileError("zsl.linalg.llt.factorIntoUnchecked: not implemented yet for o: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n"),
        },
        .dense => switch (comptime meta.matrixStorage(X)) {
            .static => try @import("factor/lltden_mat___sta.zig").factorIntoUnchecked(o, x),
            .dense => try @import("factor/lltden_mat___den.zig").factorIntoUnchecked(o, x),
            // .sparse => try @import("factor/lltden_mat___spa.zig").factorIntoUnchecked(o, x),
            .numeric => unreachable,
            else => @compileError("zsl.linalg.llt.factorIntoUnchecked: not implemented yet for o: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n"),
        },
        .sparse => switch (comptime meta.matrixStorage(X)) {
            // .sparse => try @import("factor/lltspa_mat___spa.zig").factorIntoUnchecked(o, x),
            .sparse => @compileError("zsl.linalg.llt.factorIntoUnchecked: not implemented yet for o: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n"),
            else => @compileError("zsl.linalg.llt.factorIntoUnchecked: o cannot point to a sparse llt object if the result is static or dense, got\n\to: *" ++
                @typeName(O) ++ "x: " ++ @typeName(X) ++ "\n"),
        },
        .numeric => unreachable,
    }
}
