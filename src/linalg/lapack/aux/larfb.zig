const int = @import("../../../int.zig");
const linalg = @import("../../../linalg.zig");
const matrix = @import("../../../matrix.zig");
const meta = @import("../../../meta.zig");
const numeric = @import("../../../numeric.zig");
const utils = @import("../utils.zig");

/// Applies a block reflector `H`, or its complex conjugate `Hᴴ`, to an `m × n`
/// matrix `C`, from either the left or the right.
///
/// ## Signature
/// ```zig
/// linalg.lapack.larfb(layout: matrix.Layout, side: linalg.blas.Size, trans: linalg.blas.Transpose, direct: linalg.lapack.Direction, storev: linalg.lapack.Storage, m: usize, n: usize, k: usize, v: [*]const V, ldv: usize, t: [*]const T, ldt: usize, c: [*]C, ldc: usize, work: [*]W, ldwork: usize) !void
/// ```
///
/// ## Arguments
/// * `layout` (`matrix.Layout`): Specifies whether two-dimensional array
///   storage is col-major or row-major.
/// * `side` (`linalg.blas.Side`): Specifies whether the block reflector is
///   applied from the left or the right:
///   * `left`: `C = H * C`.
///   * `right`: `C = C * H`.
/// * `transa` (`linalg.blas.Transpose`): Specifies whether the reflector or its
///   conjugate transpose is applied. (WHAT TO DO HERE?)
/// * `direct` (`linalg.lapack.Direction`): Specifies how `H` is formed from a
///   product of elementary reflectors:
///   * `forward`: `H = H₁ H₂ ⋯ Hₖ`.
///   * `backward`: `H = Hₖ ⋯ H₂ H₁`.
/// * `storev` (`linal.lapack.Storage`): Specifies how the vectors which define
///   the elementary reflectors are stored.
pub fn larfb(
    layout: matrix.Layout,
    side: linalg.blas.Side,
    trans: linalg.blas.Transpose,
    direct: linalg.lapack.Direction,
    storev: linalg.lapack.Storage,
    m: usize,
    n: usize,
    k: usize,
    v: anytype,
    ldv: usize,
    t: anytype,
    ldt: usize,
    c: anytype,
    ldc: usize,
    work: anytype,
    ldwork: usize,
) !void {
    if (m <= 0 or n <= 0)
        return;

    if (storev == .columnwise) {
        if (direct == .forward) {
            // Let v = ( v1 )    (first k rows)
            //         ( v2 )
            // where v1 is unit lower triangular.
            if (side == .left) {
                // Form h * c or h**h * c where c = ( c1 )
                //                                  ( c2 )

                // w := c**h * v = (c1**h * v1 + c2**h * v2) (stored in work)

                // w := c1**h
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        n,
                        c + utils.index(order, j, 0, ldc),
                        utils.row_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );

                    try lapack.lacgv(
                        n,
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                    );
                }

                // w := w * v1
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .no_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (m > k) {
                    // w := w + c2**h *v2
                    try blas.gemm(
                        order,
                        .conj_trans,
                        .no_trans,
                        n,
                        k,
                        m - k,
                        1,
                        c + utils.index(order, k, 0, ldc),
                        ldc,
                        v + utils.index(order, k, 0, ldv),
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t**h or w * t
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    trans.reverse(),
                    .non_unit,
                    n,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - v * w**h
                if (m > k) {
                    // c2 := c2 - v2 * w**h
                    try blas.gemm(
                        order,
                        .no_trans,
                        .conj_trans,
                        m - k,
                        n,
                        k,
                        -1,
                        v + utils.index(order, k, 0, ldv),
                        ldv,
                        work,
                        ldwork,
                        1,
                        c + utils.index(order, k, 0, ldc),
                        ldc,
                        ctx,
                    );
                }

                // w := w * v1**h
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .conj_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c1 := c1 - w**h
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < n) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, j, i, ldc)],
                            c[utils.index(order, j, i, ldc)],
                            try ops.conj(work[utils.index(order, i, j, ldwork)], ctx),
                            ctx,
                        );
                    }
                }
            } else {
                // Form c * h or c * h**h where c = ( c1 c2 )

                // w := c * v = (c1*v1 + c2*v2) (stored in work)

                // w := c1
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        m,
                        c + utils.index(order, 0, j, ldc),
                        utils.col_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );
                }

                // w := w * v1
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .no_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (n > k) {
                    // w := w + c2 * v2
                    try blas.gemm(
                        order,
                        .no_trans,
                        .no_trans,
                        m,
                        k,
                        n - k,
                        1,
                        c + utils.index(order, 0, k, ldc),
                        ldc,
                        v + utils.index(order, k, 0, ldv),
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t or w * t**h
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    trans,
                    .non_unit,
                    m,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - w * v**h
                if (n > k) {
                    // c2 := c2 - w * v2**h
                    try blas.gemm(
                        order,
                        .no_trans,
                        .conj_trans,
                        m,
                        n - k,
                        k,
                        -1,
                        work,
                        ldwork,
                        v + utils.index(order, k, 0, ldv),
                        ldv,
                        1,
                        c + utils.index(order, 0, k, ldc),
                        ldc,
                        ctx,
                    );
                }

                // w := w * v1**h
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .conj_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c1 := c1 - w
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, i, j, ldc)],
                            c[utils.index(order, i, j, ldc)],
                            work[utils.index(order, i, j, ldwork)],
                            ctx,
                        );
                    }
                }
            }
        } else {
            // Let v = ( v1 )
            //         ( v2 )    (last k rows)
            // where v2 is unit upper triangular.
            if (side == .left) {
                // Form h * c or h**h * c where c = ( c1 )
                //                                  ( c2 )

                // w := c**h * v = (c1**h * v1 + c2**h * v2) (stored in work)

                // w := c2**h
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        n,
                        c + utils.index(order, m - k + j, 0, ldc),
                        utils.row_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );

                    try lapack.lacgv(
                        n,
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                    );
                }

                // w := w * v2
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .no_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v + utils.index(order, m - k, 0, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (m > k) {
                    // w := w + c1**h * v1
                    try blas.gemm(
                        order,
                        .conj_trans,
                        .no_trans,
                        n,
                        k,
                        m - k,
                        1,
                        c,
                        ldc,
                        v,
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t**h or w * t
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    trans.reverse(),
                    .non_unit,
                    n,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - v * w**h
                if (m > k) {
                    // c1 := c1 - v1 * w**h
                    try blas.gemm(
                        order,
                        .no_trans,
                        .conj_trans,
                        m - k,
                        n,
                        k,
                        -1,
                        v,
                        ldv,
                        work,
                        ldwork,
                        1,
                        c,
                        ldc,
                        ctx,
                    );
                }

                // w := w * v2**h
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .conj_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v + utils.index(order, m - k, 0, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c2 := c2 - w**h
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < n) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, m - k + j, i, ldc)],
                            c[utils.index(order, m - k + j, i, ldc)],
                            try ops.conj(work[utils.index(order, i, j, ldwork)], ctx),
                            ctx,
                        );
                    }
                }
            } else {
                // Form c * h or c * h**h where c = ( c1 c2 )

                // w := c * v = (c1*v1 + c2*v2) (stored in work)

                // w := c2
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        m,
                        c + utils.index(order, 0, n - k + j, ldc),
                        utils.col_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );
                }

                // w := w * v2
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .no_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v + utils.index(order, 0, n - k, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (n > k) {
                    // w := w + c1 * v1
                    try blas.gemm(
                        order,
                        .no_trans,
                        .no_trans,
                        m,
                        k,
                        n - k,
                        1,
                        c,
                        ldc,
                        v,
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t or w * t**h
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    trans,
                    .non_unit,
                    m,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - w * v**h
                if (n > k) {
                    // c1 := c1 - w * v1**h
                    try blas.gemm(
                        order,
                        .no_trans,
                        .conj_trans,
                        m,
                        n - k,
                        k,
                        -1,
                        work,
                        ldwork,
                        v,
                        ldv,
                        1,
                        c,
                        ldc,
                        ctx,
                    );
                }

                // w := w * v2**h
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .conj_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v + utils.index(order, n - k, 0, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c2 := c2 - w
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, i, n - k + j, ldc)],
                            c[utils.index(order, i, n - k + j, ldc)],
                            work[utils.index(order, i, j, ldwork)],
                            ctx,
                        );
                    }
                }
            }
        }
    } else {
        if (direct == .forward) {
            // Let v = ( v1 v2 )    (v1: first k columns)
            // where v1 is unit upper triangular.
            if (side == .left) {
                // Form h * c or h**h * c where c = ( c1 )
                //                                  ( c2 )

                // w := c**h * v**h = (c1**h * v1**h + c2**h * v2**h) (stored in work)

                // w := c1**h
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        n,
                        c + utils.index(order, j, 0, ldc),
                        utils.row_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );

                    try lapack.lacgv(
                        n,
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                    );
                }

                // w := w * v1**h
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .conj_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (m > k) {
                    // w := w + c2**h * v2**h
                    try blas.gemm(
                        order,
                        .conj_trans,
                        .conj_trans,
                        n,
                        k,
                        m - k,
                        1,
                        c + utils.index(order, k, 0, ldc),
                        ldc,
                        v + utils.index(order, 0, k, ldv),
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t**h or w * t
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    trans.reverse(),
                    .non_unit,
                    n,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - v**h * w**h
                if (m > k) {
                    // c2 := c2 - v2**h * w**h
                    try blas.gemm(
                        order,
                        .conj_trans,
                        .conj_trans,
                        m - k,
                        n,
                        k,
                        -1,
                        v + utils.index(order, 0, k, ldv),
                        ldv,
                        work,
                        ldwork,
                        1,
                        c + utils.index(order, k, 0, ldc),
                        ldc,
                        ctx,
                    );
                }

                // w := w * v1
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .no_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c1 := c1 - w**h
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < n) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, j, i, ldc)],
                            c[utils.index(order, j, i, ldc)],
                            try ops.conj(work[utils.index(order, i, j, ldwork)], ctx),
                            ctx,
                        );
                    }
                }
            } else {
                // Form c * h or c * h**h where c = ( c1 c2 )

                // w := c * v**h = (c1*v1**h + c2*v2**h) (stored in work)

                // w := c1
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        m,
                        c + utils.index(order, 0, j, ldc),
                        utils.col_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );
                }

                // w := w * v1**h
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .conj_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (n > k) {
                    // w := w + c2 * v2**h
                    try blas.gemm(
                        order,
                        .no_trans,
                        .conj_trans,
                        m,
                        k,
                        n - k,
                        1,
                        c + utils.index(order, 0, k, ldc),
                        ldc,
                        v + utils.index(order, 0, k, ldv),
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t or w * t**h
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    trans,
                    .non_unit,
                    m,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - w * v
                if (n > k) {
                    // c2 := c2 - w * v2
                    try blas.gemm(
                        order,
                        .no_trans,
                        .no_trans,
                        m,
                        n - k,
                        k,
                        -1,
                        work,
                        ldwork,
                        v + utils.index(order, 0, k, ldv),
                        ldv,
                        1,
                        c + utils.index(order, 0, k, ldc),
                        ldc,
                        ctx,
                    );
                }

                // w := w * v1
                try blas.trmm(
                    order,
                    .right,
                    .upper,
                    .no_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v,
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c1 := c1 - w
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, i, j, ldc)],
                            c[utils.index(order, i, j, ldc)],
                            work[utils.index(order, i, j, ldwork)],
                            ctx,
                        );
                    }
                }
            }
        } else {
            // Let v = ( v1 v2 )    (v2: last k columns)
            // where v2 is unit lower triangular.
            if (side == .left) {
                // Form h * c or h**h * c where c = ( c1 )
                //                                  ( c2 )

                // w := c**h * v**h = (c1**h * v1**h + c2**h * v2**h) (stored in work)

                // w := c2**h
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        n,
                        c + utils.index(order, m - k + j, 0, ldc),
                        utils.row_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );

                    try lapack.lacgv(
                        n,
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                    );
                }

                // w := w * v2**h
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .conj_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v + utils.index(order, 0, m - k, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (m > k) {
                    // w := w + c1**h * v1**h
                    try blas.gemm(
                        order,
                        .conj_trans,
                        .conj_trans,
                        n,
                        k,
                        m - k,
                        1,
                        c,
                        ldc,
                        v,
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t**h or w * t
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    trans.reverse(),
                    .non_unit,
                    n,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - v**h * w**h
                if (m > k) {
                    // c1 := c1 - v1**h * w**h
                    try blas.gemm(
                        order,
                        .conj_trans,
                        .conj_trans,
                        m - k,
                        n,
                        k,
                        -1,
                        v,
                        ldv,
                        work,
                        ldwork,
                        1,
                        c,
                        ldc,
                        ctx,
                    );
                }

                // w := w * v2
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .no_trans,
                    .unit,
                    n,
                    k,
                    1,
                    v + utils.index(order, 0, m - k, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c2 := c2 - w**h
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < n) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, m - k + j, i, ldc)],
                            c[utils.index(order, m - k + j, i, ldc)],
                            try ops.conj(work[utils.index(order, i, j, ldwork)], ctx),
                            ctx,
                        );
                    }
                }
            } else {
                // Form c * h or c * h**h where c = ( c1 c2 )

                // w := c * v**h = (c1*v1**h + c2*v2**h) (stored in work)

                // w := c2
                var j: i32 = 0;
                while (j < k) : (j += 1) {
                    try blas.copy(
                        m,
                        c + utils.index(order, 0, n - k + j, ldc),
                        utils.col_ld(order, ldc),
                        work + utils.index(order, 0, j, ldwork),
                        utils.col_ld(order, ldwork),
                        ctx,
                    );
                }

                // w := w * v2**h
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .conj_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v + utils.index(order, 0, n - k, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                if (n > k) {
                    // w := w + c1 * v1**h
                    try blas.gemm(
                        order,
                        .no_trans,
                        .conj_trans,
                        m,
                        k,
                        n - k,
                        1,
                        c,
                        ldc,
                        v,
                        ldv,
                        1,
                        work,
                        ldwork,
                        ctx,
                    );
                }

                // w := w * t or w * t**h
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    trans,
                    .non_unit,
                    m,
                    k,
                    1,
                    t,
                    ldt,
                    work,
                    ldwork,
                    ctx,
                );

                // c := c - w * v
                if (n > k) {
                    // c1 := c1 - w * v1
                    try blas.gemm(
                        order,
                        .no_trans,
                        .no_trans,
                        m,
                        n - k,
                        k,
                        -1,
                        work,
                        ldwork,
                        v,
                        ldv,
                        1,
                        c,
                        ldc,
                        ctx,
                    );
                }

                // w := w * v2
                try blas.trmm(
                    order,
                    .right,
                    .lower,
                    .no_trans,
                    .unit,
                    m,
                    k,
                    1,
                    v + utils.index(order, 0, n - k, ldv),
                    ldv,
                    work,
                    ldwork,
                    ctx,
                );

                // c1 := c1 - w
                j = 0;
                while (j < k) : (j += 1) {
                    var i: i32 = 0;
                    while (i < m) : (i += 1) {
                        try ops.sub_(
                            &c[utils.index(order, i, n - k + j, ldc)],
                            c[utils.index(order, i, n - k + j, ldc)],
                            work[utils.index(order, i, j, ldwork)],
                            ctx,
                        );
                    }
                }
            }
        }
    }

    return;
}
