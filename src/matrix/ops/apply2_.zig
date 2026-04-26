const meta = @import("../../meta.zig");

const matrix = @import("../../matrix.zig");

/// Applies a binary in-place operation elementwise between an output and two
/// input matrices, or between an output matrix, an input matrix and an input
/// numeric.
///
/// Exact aliasing (in-place modification) between the output and an input is
/// permitted and often more efficient. Any other form of memory overlap might
/// yield incorrect results.
///
/// For two input sparse matrices, or an input sparse matrix and an input
/// numeric, if the output is also a sparse matrix, the operation is only
/// applied to the indices where at least one of the matrices has a non-zero
/// element.
///
/// ## Signature
/// ```zig
/// matrix.apply2_(*O, x: X, y: Y, op_: Op) !void
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
/// * `matrix.Error.DimensionMismatch`: If the matrices do not have the same
///   dimensions.
pub fn apply2_(o: anytype, x: anytype, y: anytype, comptime op_: anytype) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);
    const Op: type = @TypeOf(op_);
    const opinfo = @typeInfo(Op);

    comptime if (!meta.isPointer(O) or meta.isConstPointer(O) or !meta.isMatrix(meta.Child(O)) or
        (!meta.isMatrix(X) and !meta.isNumeric(X)) or (!meta.isMatrix(Y) and !meta.isNumeric(Y)) or
        (!meta.isMatrix(X) and !meta.isMatrix(Y)) or
        opinfo != .@"fn" or opinfo.@"fn".params.len != 3)
        @compileError("zsl.matrix.apply2_: o must be a mutable one-itme pointer to a matrix, at least one of x or y must be a matrix, the other must be a matrix or a numeric, and op_ must be a function of three arguments, got\n\to: " ++
            @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

    O = meta.Child(O);

    const x_rows = if (comptime meta.isMatrix(X)) x.rows else o.rows;
    const x_cols = if (comptime meta.isMatrix(X)) x.cols else o.cols;

    const y_rows = if (comptime meta.isMatrix(Y)) y.rows else o.rows;
    const y_cols = if (comptime meta.isMatrix(Y)) y.cols else o.cols;

    if (o.rows != x_rows or o.cols != x_cols or
        o.rows != y_rows or o.cols != y_cols)
        return matrix.Error.DimensionMismatch;

    switch (comptime meta.matrixType(O)) {
        .general_dense => switch (comptime meta.matrixType(X)) {
            .general_dense => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdgdgd.zig").apply2_(o, x, y, op_),
                .general_sparse => return @import("apply2_/gdgdgs.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdgdsd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/gdgdss.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdgdhd.zig").apply2_(o, x, y, op_),
                .hermitian_sparse => return @import("apply2_/gdgdhs.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdgdtd.zig").apply2_(o, x, y, op_),
                .triangular_sparse => return @import("apply2_/gdgdts.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gdgddi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gdgdpe.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/gdgdnu.zig").apply2_(o, x, y, op_),
            },
            .general_sparse => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdgsgd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdgssd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdgshd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdgstd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                // .diagonal => return @import("apply2_/gd_sdi.zig").apply2_(o, x, y, op_),
                // .permutation => return @import("apply2_/gd_spe.zig").apply2_(o, x, y, op_),
                // .numeric => return @import("apply2_/gd_snu.zig").apply2_(o, x, y, op_),
                else => @compileError("Not implemented yet"),
            },
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdsdgd.zig").apply2_(o, x, y, op_),
                .general_sparse => return @import("apply2_/gdsdgs.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdsdsd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/gdsdss.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdsdhd.zig").apply2_(o, x, y, op_),
                .hermitian_sparse => return @import("apply2_/gdsdhs.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdsdtd.zig").apply2_(o, x, y, op_),
                .triangular_sparse => return @import("apply2_/gdsdts.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gdsddi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gdsdpe.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/gdsdnu.zig").apply2_(o, x, y, op_),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdssgd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdsssd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdsshd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdsstd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                // .diagonal => return @import("apply2_/gd_sdi.zig").apply2_(o, x, y, op_),
                // .permutation => return @import("apply2_/gd_spe.zig").apply2_(o, x, y, op_),
                // .numeric => return @import("apply2_/gd_snu.zig").apply2_(o, x, y, op_),
                else => @compileError("Not implemented yet"),
            },
            .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdhdgd.zig").apply2_(o, x, y, op_),
                .general_sparse => return @import("apply2_/gdhdgs.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdhdsd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/gdhdss.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdhdhd.zig").apply2_(o, x, y, op_),
                .hermitian_sparse => return @import("apply2_/gdhdhs.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdhdtd.zig").apply2_(o, x, y, op_),
                .triangular_sparse => return @import("apply2_/gdhdts.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gdhddi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gdhdpe.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/gdhdnu.zig").apply2_(o, x, y, op_),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdhsgd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdhssd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdhshd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdhstd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                // .diagonal => return @import("apply2_/gd_sdi.zig").apply2_(o, x, y, op_),
                // .permutation => return @import("apply2_/gd_spe.zig").apply2_(o, x, y, op_),
                // .numeric => return @import("apply2_/gd_snu.zig").apply2_(o, x, y, op_),
                else => @compileError("Not implemented yet"),
            },
            .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdtdgd.zig").apply2_(o, x, y, op_),
                .general_sparse => return @import("apply2_/gdtdgs.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdtdsd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/gdtdss.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdtdhd.zig").apply2_(o, x, y, op_),
                .hermitian_sparse => return @import("apply2_/gdtdhs.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdtdtd.zig").apply2_(o, x, y, op_),
                .triangular_sparse => return @import("apply2_/gdtdts.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gdtddi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gdtdpe.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/gdtdnu.zig").apply2_(o, x, y, op_),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdtsgd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdtssd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdtshd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdtstd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/g__s_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                // .diagonal => return @import("apply2_/gd_sdi.zig").apply2_(o, x, y, op_),
                // .permutation => return @import("apply2_/gd_spe.zig").apply2_(o, x, y, op_),
                // .numeric => return @import("apply2_/gd_snu.zig").apply2_(o, x, y, op_),
                else => @compileError("Not implemented yet"),
            },
            .builder_sparse => unreachable,
            .diagonal => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gddigd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/gddi_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gddisd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/gddi_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gddihd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/gddi_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdditd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/gddi_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gddidi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gddipe.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/gddinu.zig").apply2_(o, x, y, op_),
                else => @compileError("Not implemented yet"),
            },
            .permutation => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdpegd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/gdpe_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdpesd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/gdpe_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdpehd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/gdpe_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdpetd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/gdpe_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gdpedi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gdpepe.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/gdpenu.zig").apply2_(o, x, y, op_),
                else => @compileError("Not implemented yet"),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .general_dense => return @import("apply2_/gdnugd.zig").apply2_(o, x, y, op_),
                // .general_sparse => return @import("apply2_/gdnu_s.zig").apply2_(o, x, y, op_),
                .symmetric_dense => return @import("apply2_/gdnusd.zig").apply2_(o, x, y, op_),
                // .symmetric_sparse => return @import("apply2_/gdnu_s.zig").apply2_(o, x, y, op_),
                .hermitian_dense => return @import("apply2_/gdnuhd.zig").apply2_(o, x, y, op_),
                // .hermitian_sparse => return @import("apply2_/gdnu_s.zig").apply2_(o, x, y, op_),
                .triangular_dense => return @import("apply2_/gdnutd.zig").apply2_(o, x, y, op_),
                // .triangular_sparse => return @import("apply2_/gdnu_s.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/gdnudi.zig").apply2_(o, x, y, op_),
                .permutation => return @import("apply2_/gdnupe.zig").apply2_(o, x, y, op_),
                .numeric => unreachable,
                else => @compileError("Not implemented yet"),
            },
        },
        .general_sparse => @compileError("Not implemented yet"),
        .symmetric_dense => switch (comptime meta.matrixType(X)) {
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => return @import("apply2_/sdsdsd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/sdsdss.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/sdsddi.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/sdsdnu.zig").apply2_(o, x, y, op_),
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => return @import("apply2_/sdsssd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/sdssss.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/sdssdi.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/sdssnu.zig").apply2_(o, x, y, op_),
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .diagonal => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => return @import("apply2_/sddisd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/sddiss.zig").apply2_(o, x, y, op_),
                .diagonal => return @import("apply2_/sddidi.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .numeric => return @import("apply2_/sddinu.zig").apply2_(o, x, y, op_),
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => return @import("apply2_/sdnusd.zig").apply2_(o, x, y, op_),
                .symmetric_sparse => return @import("apply2_/sdnuss.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/sdnudi.zig").apply2_(o, x, y, op_),
                .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
        },
        .symmetric_sparse => @compileError("Not implemented yet"),
        .hermitian_dense => switch (comptime meta.matrixType(X)) {
            .symmetric_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdsdsd.zig").apply2_(o, x, y, op_);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdsdss.zig").apply2_(o, x, y, op_);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdsdhd.zig").apply2_(o, x, y, op_);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdsdhs.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdsddi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdsdnu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .symmetric_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdsssd.zig").apply2_(o, x, y, op_);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdssss.zig").apply2_(o, x, y, op_);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdsshd.zig").apply2_(o, x, y, op_);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdsshs.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdssdi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdssnu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .hermitian_dense => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhdsd.zig").apply2_(o, x, y, op_);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhdss.zig").apply2_(o, x, y, op_);
                },
                .hermitian_dense => return @import("apply2_/hdhdhd.zig").apply2_(o, x, y, op_),
                .hermitian_sparse => return @import("apply2_/hdhdhs.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhddi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhdnu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .hermitian_sparse => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhssd.zig").apply2_(o, x, y, op_);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhsss.zig").apply2_(o, x, y, op_);
                },
                .hermitian_dense => return @import("apply2_/hdhshd.zig").apply2_(o, x, y, op_),
                .hermitian_sparse => return @import("apply2_/hdhshs.zig").apply2_(o, x, y, op_),
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhsdi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.isComplex(Y))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdhsnu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .builder_sparse => unreachable,
            .diagonal => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sddisd.zig").apply2_(o, x, y, op_);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sddiss.zig").apply2_(o, x, y, op_);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hddihd.zig").apply2_(o, x, y, op_);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(meta.Numeric(X)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hddihs.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sddidi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sddinu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .symmetric_dense => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdnusd.zig").apply2_(o, x, y, op_);
                },
                .symmetric_sparse => {
                    comptime if (meta.isComplex(X) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdnuss.zig").apply2_(o, x, y, op_);
                },
                .hermitian_dense => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdnuhd.zig").apply2_(o, x, y, op_);
                },
                .hermitian_sparse => {
                    comptime if (meta.isComplex(X))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/hdnuhs.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.isComplex(meta.Numeric(X)) or meta.isComplex(meta.Numeric(Y)))
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/sdnudi.zig").apply2_(o, x, y, op_);
                },
                .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
        },
        .hermitian_sparse => @compileError("Not implemented yet"),
        .triangular_dense => switch (comptime meta.matrixType(X)) {
            .triangular_dense => switch (comptime meta.matrixType(Y)) {
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtdtd.zig").apply2_(o, x, y, op_);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtdts.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtddi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtdnu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .triangular_sparse => switch (comptime meta.matrixType(Y)) {
                .triangular_dense => {
                    comptime if (meta.uploOf(X) != meta.uploOf(Y) or meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtstd.zig").apply2_(o, x, y, op_);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(X) != meta.uploOf(Y) or meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtsts.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtsdi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.uploOf(O) != meta.uploOf(X) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdtsnu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .builder_sparse => unreachable,
            .diagonal => switch (comptime meta.matrixType(Y)) {
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdditd.zig").apply2_(o, x, y, op_);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tddits.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tddidi.zig").apply2_(o, x, y, op_);
                },
                .numeric => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tddinu.zig").apply2_(o, x, y, op_);
                },
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .triangular_dense => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdnutd.zig").apply2_(o, x, y, op_);
                },
                .triangular_sparse => {
                    comptime if (meta.uploOf(O) != meta.uploOf(Y) or meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdnuts.zig").apply2_(o, x, y, op_);
                },
                .builder_sparse => unreachable,
                .diagonal => {
                    comptime if (meta.diagOf(O) == .unit)
                        @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n");

                    return @import("apply2_/tdnudi.zig").apply2_(o, x, y, op_);
                },
                .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
        },
        .triangular_sparse => @compileError("Not implemented yet"),
        .builder_sparse => @compileError("Not implemented yet"),
        .diagonal => switch (comptime meta.matrixType(X)) {
            .builder_sparse => unreachable,
            .diagonal => switch (comptime meta.matrixType(Y)) {
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/dididi.zig").apply2_(o, x, y, op_),
                .numeric => return @import("apply2_/didinu.zig").apply2_(o, x, y, op_),
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            .numeric => switch (comptime meta.matrixType(Y)) {
                .builder_sparse => unreachable,
                .diagonal => return @import("apply2_/dinudi.zig").apply2_(o, x, y, op_),
                .numeric => unreachable,
                else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
            },
            else => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
        },
        .permutation => @compileError("zsl.matrix.apply2_: the result of the operation is incompatible with o's type, got\n\to: *" ++ @typeName(O) ++ "\n\tx: " ++ @typeName(X) ++ "\n\ty: " ++ @typeName(Y) ++ "\n\top_: " ++ @typeName(Op) ++ "\n"),
        .numeric => unreachable,
    }
}
