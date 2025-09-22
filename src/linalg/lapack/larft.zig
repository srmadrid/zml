const std = @import("std");

const types = @import("../../types.zig");
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Order = types.Order;
const Direction = lapack.Direction;
const Storage = lapack.Storage;

const utils = @import("../utils.zig");

pub fn larft(
    order: Order,
    direct: Direction,
    storev: Storage,
    n: i32,
    k: i32,
    v: anytype,
    ldv: i32,
    tau: anytype,
    t: anytype,
    ldt: i32,
    ctx: anytype,
) !void {
    // The general scheme used is inspired by the approach inside dgeqrt3
    // which was (at the time of writing this code):
    // based on the algorithm of elmroth and gustavson,
    // ibm j. res. develop. vol 44 no. 4 july 2000.

    if (n == 0 or k == 0)
        return;

    if (n == 1 or k == 1) {
        try ops.set(
            &t[utils.index(order, 0, 0, ldt)],
            tau[0],
            ctx,
        );

        return;
    }

    const l: i32 = int.div(k, 2);

    // Determine what kind of q we need to compute
    // qr happens when we have forward direction in column storage
    const qr: bool = direct == .forward and storev == .columnwise;

    // lq happens when we have forward direction in row storage
    const lq: bool = direct == .forward and storev == .rowwise;

    // ql happens when we have backward direction in column storage
    const ql: bool = direct == .backward and storev == .columnwise;

    // The last case is rq. due to how we structured this, if the
    // above 3 are false, then rq must be true, so we never store
    // this

    if (qr) {
        // Break v apart into 6 components
        //
        // v = |---------------|
        //     |v_{1,1} 0      |
        //     |v_{2,1} v_{2,2}|
        //     |v_{3,1} v_{3,2}|
        //     |---------------|
        //
        // v_{1,1}\in\c^{l,l}      unit lower triangular
        // v_{2,1}\in\c^{k-l,l}    rectangular
        // v_{3,1}\in\c^{n-k,l}    rectangular
        //
        // v_{2,2}\in\c^{k-l,k-l}  unit lower triangular
        // v_{3,2}\in\c^{n-k,k-l}  rectangular
        //
        // we will construct the t matrix
        // t = |---------------|
        //     |t_{1,1} t_{1,2}|
        //     |0       t_{2,2}|
        //     |---------------|
        //
        // t is the triangular factor obtained from block reflectors.
        // to motivate the structure, assume we have already computed t_{1,1}
        // and t_{2,2}. then collect the associated reflectors in v_1 and v_2
        //
        // t_{1,1}\in\c^{l, l}     upper triangular
        // t_{2,2}\in\c^{k-l, k-l} upper triangular
        // t_{1,2}\in\c^{l, k-l}   rectangular
        //
        // where l = floor(k/2)
        //
        // then, consider the product:
        //
        // (i - v_1*t_{1,1}*v_1')*(i - v_2*t_{2,2}*v_2')
        // = i - v_1*t_{1,1}*v_1' - v_2*t_{2,2}*v_2' + v_1*t_{1,1}*v_1'*v_2*t_{2,2}*v_2'
        //
        // define t{1,2} = -t_{1,1}*v_1'*v_2*t_{2,2}
        //
        // then, we can define the matrix v as
        // v = |-------|
        //     |v_1 v_2|
        //     |-------|
        //
        // so, our product is equivalent to the matrix product
        // i - v*t*v'
        // this means, we can compute t_{1,1} and t_{2,2}, then use this information
        // to compute t_{1,2}

        // Compute t_{1,1} recursively
        try larft(
            order,
            direct,
            storev,
            n,
            l,
            v,
            ldv,
            tau,
            t,
            ldt,
            ctx,
        );

        // compute t_{2,2} recursively
        try larft(
            order,
            direct,
            storev,
            n - l,
            k - l,
            v + utils.index(order, l, l, ldv),
            ldv,
            tau + types.scast(u32, l),
            t + utils.index(order, l, l, ldt),
            ldt,
            ctx,
        );

        // Compute t_{1,2}
        // t_{1,2} = v_{2,1}'
        var j: i32 = 0;
        while (j < l) : (j += 1) {
            var i: i32 = 0;
            while (i < k - l) : (i += 1) {
                try ops.set(
                    &t[utils.index(order, j, l + i, ldt)],
                    try ops.conjugate(v[utils.index(order, l + i, j, ldv)], ctx),
                    ctx,
                );
            }
        }

        // t_{1,2} = t_{1,2}*v_{2,2}
        try blas.trmm(
            order,
            .right,
            .lower,
            .no_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(order, l, l, ldv),
            ldv,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // t_{1,2} = v_{3,1}'*v_{3,2} + t_{1,2}
        // note: we assume k <= n, and gemm will do nothing if n=k
        try blas.gemm(
            order,
            .conj_trans,
            .no_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(order, k, 0, ldv),
            ldv,
            v + utils.index(order, k, l, ldv),
            ldv,
            1,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // At this point, we have that t_{1,2} = v_1'*v_2
        // all that is left is to pre and post multiply by -t_{1,1} and t_{2,2}
        // respectively.

        // t_{1,2} = -t_{1,1}*t_{1,2}
        try blas.trmm(
            order,
            .left,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t,
            ldt,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // t_{1,2} = t_{1,2}*t_{2,2}
        try blas.trmm(
            order,
            .right,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t + utils.index(order, l, l, ldt),
            ldt,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );
    } else if (lq) {
        // Break v apart into 6 components
        //
        // v = |----------------------|
        //     |v_{1,1} v_{1,2} v{1,3}|
        //     |0       v_{2,2} v{2,3}|
        //     |----------------------|
        //
        // v_{1,1}\in\c^{l,l}      unit upper triangular
        // v_{1,2}\in\c^{l,k-l}    rectangular
        // v_{1,3}\in\c^{l,n-k}    rectangular
        //
        // v_{2,2}\in\c^{k-l,k-l}  unit upper triangular
        // v_{2,3}\in\c^{k-l,n-k}  rectangular
        //
        // where l = floor(k/2)
        //
        // we will construct the t matrix
        // t = |---------------|
        //     |t_{1,1} t_{1,2}|
        //     |0       t_{2,2}|
        //     |---------------|
        //
        // t is the triangular factor obtained from block reflectors.
        // to motivate the structure, assume we have already computed t_{1,1}
        // and t_{2,2}. then collect the associated reflectors in v_1 and v_2
        //
        // t_{1,1}\in\c^{l, l}     upper triangular
        // t_{2,2}\in\c^{k-l, k-l} upper triangular
        // t_{1,2}\in\c^{l, k-l}   rectangular
        //
        // then, consider the product:
        //
        // (i - v_1'*t_{1,1}*v_1)*(i - v_2'*t_{2,2}*v_2)
        // = i - v_1'*t_{1,1}*v_1 - v_2'*t_{2,2}*v_2 + v_1'*t_{1,1}*v_1*v_2'*t_{2,2}*v_2
        //
        // define t_{1,2} = -t_{1,1}*v_1*v_2'*t_{2,2}
        //
        // then, we can define the matrix v as
        // v = |---|
        //     |v_1|
        //     |v_2|
        //     |---|
        //
        // so, our product is equivalent to the matrix product
        // i - v'*t*v
        // this means, we can compute t_{1,1} and t_{2,2}, then use this information
        // to compute t_{1,2}

        // Compute t_{1,1} recursively
        try larft(
            order,
            direct,
            storev,
            n,
            l,
            v,
            ldv,
            tau,
            t,
            ldt,
            ctx,
        );

        // Compute t_{2,2} recursively
        try larft(
            order,
            direct,
            storev,
            n - l,
            k - l,
            v + utils.index(order, l, l, ldv),
            ldv,
            tau + types.scast(u32, l),
            t + utils.index(order, l, l, ldt),
            ldt,
            ctx,
        );

        // Compute t_{1,2}
        // t_{1,2} = v_{1,2}
        try lapack.lacpy(
            order,
            .{ .full = {} },
            l,
            k - l,
            v + utils.index(order, 0, l, ldv),
            ldv,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // t_{1,2} = t_{1,2}*v_{2,2}'
        try blas.trmm(
            order,
            .right,
            .upper,
            .conj_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(order, l, l, ldv),
            ldv,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // t_{1,2} = v_{1,3}*v_{2,3}' + t_{1,2}
        // note: we assume k <= n, and gemm will do nothing if n=k
        try blas.gemm(
            order,
            .no_trans,
            .conj_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(order, 0, k, ldv),
            ldv,
            v + utils.index(order, l, k, ldv),
            ldv,
            1,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // At this point, we have that t_{1,2} = v_1*v_2'
        // all that is left is to pre and post multiply by -t_{1,1} and t_{2,2}
        // respectively.

        // t_{1,2} = -t_{1,1}*t_{1,2}
        try blas.trmm(
            order,
            .left,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t,
            ldt,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );

        // t_{1,2} = t_{1,2}*t_{2,2}
        try blas.trmm(
            order,
            .right,
            .upper,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t + utils.index(order, l, l, ldt),
            ldt,
            t + utils.index(order, 0, l, ldt),
            ldt,
            ctx,
        );
    } else if (ql) {
        // Break v apart into 6 components
        //
        // v = |---------------|
        //     |v_{1,1} v_{1,2}|
        //     |v_{2,1} v_{2,2}|
        //     |0       v_{3,2}|
        //     |---------------|
        //
        // v_{1,1}\in\c^{n-k,k-l}  rectangular
        // v_{2,1}\in\c^{k-l,k-l}  unit upper triangular
        //
        // v_{1,2}\in\c^{n-k,l}    rectangular
        // v_{2,2}\in\c^{k-l,l}    rectangular
        // v_{3,2}\in\c^{l,l}      unit upper triangular
        //
        // we will construct the t matrix
        // t = |---------------|
        //     |t_{1,1} 0      |
        //     |t_{2,1} t_{2,2}|
        //     |---------------|
        //
        // t is the triangular factor obtained from block reflectors.
        // to motivate the structure, assume we have already computed t_{1,1}
        // and t_{2,2}. then collect the associated reflectors in v_1 and v_2
        //
        // t_{1,1}\in\c^{k-l, k-l} non-unit lower triangular
        // t_{2,2}\in\c^{l, l}     non-unit lower triangular
        // t_{2,1}\in\c^{k-l, l}   rectangular
        //
        // where l = floor(k/2)
        //
        // then, consider the product:
        //
        // (i - v_2*t_{2,2}*v_2')*(i - v_1*t_{1,1}*v_1')
        // = i - v_2*t_{2,2}*v_2' - v_1*t_{1,1}*v_1' + v_2*t_{2,2}*v_2'*v_1*t_{1,1}*v_1'
        //
        // define t_{2,1} = -t_{2,2}*v_2'*v_1*t_{1,1}
        //
        // then, we can define the matrix v as
        // v = |-------|
        //     |v_1 v_2|
        //     |-------|
        //
        // so, our product is equivalent to the matrix product
        // i - v*t*v'
        // this means, we can compute t_{1,1} and t_{2,2}, then use this information
        // to compute t_{2,1}

        // Compute t_{1,1} recursively
        try larft(
            order,
            direct,
            storev,
            n - l,
            k - l,
            v,
            ldv,
            tau,
            t,
            ldt,
            ctx,
        );

        // Compute t_{2,2} recursively
        try larft(
            order,
            direct,
            storev,
            n,
            l,
            v + utils.index(order, 0, k - l, ldv),
            ldv,
            tau + types.scast(u32, k - l),
            t + utils.index(order, k - l, k - l, ldt),
            ldt,
            ctx,
        );

        // Compute t_{2,1}
        // t_{2,1} = v_{2,2}'
        var j: i32 = 0;
        while (j < k - l) : (j += 1) {
            var i: i32 = 0;
            while (i < l) : (i += 1) {
                try ops.set(
                    &t[utils.index(order, k - l + i, j, ldt)],
                    try ops.conjugate(v[utils.index(order, n - k + j, k - l + i, ldv)], ctx),
                    ctx,
                );
            }
        }

        // t_{2,1} = t_{2,1}*v_{2,1}
        try blas.trmm(
            order,
            .right,
            .upper,
            .no_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(order, n - k, 0, ldv),
            ldv,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // t_{2,1} = v_{2,2}'*v_{2,1} + t_{2,1}
        // note: we assume k <= n, and gemm will do nothing if n=k
        try blas.gemm(
            order,
            .conj_trans,
            .no_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(order, 0, k - l, ldv),
            ldv,
            v,
            ldv,
            1,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // At this point, we have that t_{2,1} = v_2'*v_1
        // all that is left is to pre and post multiply by -t_{2,2} and t_{1,1}
        // respectively.

        // t_{2,1} = -t_{2,2}*t_{2,1}
        try blas.trmm(
            order,
            .left,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t + utils.index(order, k - l, k - l, ldt),
            ldt,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // t_{2,1} = t_{2,1}*t_{1,1}
        try blas.trmm(
            order,
            .right,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t,
            ldt,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );
    } else {
        // Else means rq case
        //
        // break v apart into 6 components
        //
        // v = |-----------------------|
        //     |v_{1,1} v_{1,2} 0      |
        //     |v_{2,1} v_{2,2} v_{2,3}|
        //     |-----------------------|
        //
        // v_{1,1}\in\c^{k-l,n-k}  rectangular
        // v_{1,2}\in\c^{k-l,k-l}  unit lower triangular
        //
        // v_{2,1}\in\c^{l,n-k}    rectangular
        // v_{2,2}\in\c^{l,k-l}    rectangular
        // v_{2,3}\in\c^{l,l}      unit lower triangular
        //
        // we will construct the t matrix
        // t = |---------------|
        //     |t_{1,1} 0      |
        //     |t_{2,1} t_{2,2}|
        //     |---------------|
        //
        // t is the triangular factor obtained from block reflectors.
        // to motivate the structure, assume we have already computed t_{1,1}
        // and t_{2,2}. then collect the associated reflectors in v_1 and v_2
        //
        // t_{1,1}\in\c^{k-l, k-l} non-unit lower triangular
        // t_{2,2}\in\c^{l, l}     non-unit lower triangular
        // t_{2,1}\in\c^{k-l, l}   rectangular
        //
        // where l = floor(k/2)
        //
        // then, consider the product:
        //
        // (i - v_2'*t_{2,2}*v_2)*(i - v_1'*t_{1,1}*v_1)
        // = i - v_2'*t_{2,2}*v_2 - v_1'*t_{1,1}*v_1 + v_2'*t_{2,2}*v_2*v_1'*t_{1,1}*v_1
        //
        // define t_{2,1} = -t_{2,2}*v_2*v_1'*t_{1,1}
        //
        // then, we can define the matrix v as
        // v = |---|
        //     |v_1|
        //     |v_2|
        //     |---|
        //
        // so, our product is equivalent to the matrix product
        // i - v'*t*v
        // this means, we can compute t_{1,1} and t_{2,2}, then use this information
        // to compute t_{2,1}

        // Compute t_{1,1} recursively
        try larft(
            order,
            direct,
            storev,
            n - l,
            k - l,
            v,
            ldv,
            tau,
            t,
            ldt,
            ctx,
        );

        // Compute t_{2,2} recursively
        try larft(
            order,
            direct,
            storev,
            n,
            l,
            v + utils.index(order, k - l, 0, ldv),
            ldv,
            tau + types.scast(u32, k - l),
            t + utils.index(order, k - l, k - l, ldt),
            ldt,
            ctx,
        );

        // Compute t_{2,1}
        // t_{2,1} = v_{2,2}
        try lapack.lacpy(
            order,
            .{ .full = {} },
            l,
            k - l,
            v + utils.index(order, k - l, n - k, ldv),
            ldv,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // t_{2,1} = t_{2,1}*v_{1,2}'
        try blas.trmm(
            order,
            .right,
            .lower,
            .conj_trans,
            .unit,
            l,
            k - l,
            1,
            v + utils.index(order, 0, n - k, ldv),
            ldv,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // t_{2,1} = v_{2,1}*v_{1,1}' + t_{2,1}
        // note: we assume k <= n, and gemm will do nothing if n=k
        try blas.gemm(
            order,
            .no_trans,
            .conj_trans,
            l,
            k - l,
            n - k,
            1,
            v + utils.index(order, k - l, 0, ldv),
            ldv,
            v,
            ldv,
            1,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // At this point, we have that t_{2,1} = v_2*v_1'
        // all that is left is to pre and post multiply by -t_{2,2} and t_{1,1}
        // respectively.

        // t_{2,1} = -t_{2,2}*t_{2,1}
        try blas.trmm(
            order,
            .left,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            -1,
            t + utils.index(order, k - l, k - l, ldt),
            ldt,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );

        // t_{2,1} = t_{2,1}*t_{1,1}
        try blas.trmm(
            order,
            .right,
            .lower,
            .no_trans,
            .non_unit,
            l,
            k - l,
            1,
            t,
            ldt,
            t + utils.index(order, k - l, 0, ldt),
            ldt,
            ctx,
        );
    }
}
