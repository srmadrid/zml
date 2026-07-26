const linalg = @import("../linalg.zig");
const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

pub fn Trace(comptime X: type) type {
    comptime if (!meta.isMatrix(X))
        @compileError("zsl.linalg.Trace: X must be a matrix type, got\n\tX = " ++
            @typeName(X) ++ "\n");

    comptime if (meta.isBuilderMatrix(X))
        @compileError("zsl.linalg.Trace: X must not be a builder matrix type, got\n\tX = " ++
            @typeName(X) ++ "\n");

    comptime if (meta.isStaticMatrix(X)) {
        if (X.rows != X.cols)
            @compileError("zsl.linalg.Trace: static matrix type X must be square, got\n\tX = " ++
                @typeName(X) ++ "\n");
    };

    return numeric.Add(meta.Numeric(X), meta.Numeric(X));
}

/// Computes the trace of a matrix, tr(x) = ∑ᵢ xᵢᵢ.
///
/// ## Signature
/// ```zig
/// linalg.trace(x: X) !linalg.Trace(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The matrix.
///
/// ## Returns
/// `linalg.Trace(@TypeOf(x))`: The trace tr(x).
///
/// ## Errors
/// * `linalg.Error.DimensionMismatch`: If the input is not square.
pub fn trace(x: anytype) !linalg.Trace(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    if (comptime !meta.isStaticMatrix(X)) {
        if (x.rows != x.cols)
            return linalg.Error.DimensionMismatch;
    }

    return traceUnchecked(x);
}

/// Computes the trace of a matrix, tr(x) = ∑ᵢ xᵢᵢ, without performing dimension
/// checks.
///
/// ## Signature
/// ```zig
/// linalg.traceUnchecked(x: X) linalg.Trace(X)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The matrix.
///
/// ## Returns
/// `linalg.Trace(@TypeOf(x))`: The trace tr(x).
pub fn traceUnchecked(x: anytype) linalg.Trace(@TypeOf(x)) {
    const X: type = @TypeOf(x);

    return switch (comptime meta.matrixType(X)) {
        .general_static => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matgensta.zig").traceUnchecked(x),
        .general_dense => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matgenden.zig").traceUnchecked(x),
        .general_sparse => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matgenspa.zig").traceUnchecked(x),
        .symmetric_static => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matsymsta.zig").traceUnchecked(x),
        .symmetric_dense => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matsymsen.zig").traceUnchecked(x),
        .symmetric_sparse => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matsymspa.zig").traceUnchecked(x),
        .hermitian_static => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/mathersta.zig").traceUnchecked(x),
        .hermitian_dense => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matherden.zig").traceUnchecked(x),
        .hermitian_sparse => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matherspa.zig").traceUnchecked(x),
        .triangular_static => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/mattrista.zig").traceUnchecked(x),
        .triangular_dense => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/mattriden.zig").traceUnchecked(x),
        .triangular_sparse => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/mattrispa.zig").traceUnchecked(x),
        .diagonal_static => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matdiasta.zig").traceUnchecked(x),
        .diagonal_sparse => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matdiaspa.zig").traceUnchecked(x),
        .permutation_static => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matpersta.zig").traceUnchecked(x),
        .permutation_sparse => @import("trace/slow.zig").traceUnchecked(x), // @import("trace/matperspa.zig").traceUnchecked(x),
        else => unreachable,
    };
}
