const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const blas = @import("../blas.zig");

pub fn rotmg(
    d1: anytype,
    d2: anytype,
    x1: anytype,
    y1: anytype,
    param: anytype,
    ctx: anytype,
) !void {
    const D1: type = types.Child(@TypeOf(d1));
    const D2: type = types.Child(@TypeOf(d2));
    const X1: type = types.Child(@TypeOf(x1));
    const Y1: type = @TypeOf(y1);
    const P: type = types.Child(@TypeOf(param));
    const C: type = types.Coerce(D1, types.Coerce(D2, types.Coerce(X1, types.Coerce(Y1, P))));

    if (comptime types.isArbitraryPrecision(D1) or
        types.isArbitraryPrecision(D2) or
        types.isArbitraryPrecision(X1) or
        types.isArbitraryPrecision(Y1) or
        types.isArbitraryPrecision(P))
    {
        @compileError("zml.linalg.blas.rotmg not implemented for arbitrary precision types yet");
    } else {
        const gam: types.EnsureFloat(C) = 4096;
        const gamsq: types.EnsureFloat(C) = 16777216;
        const rgamsq: types.EnsureFloat(C) = 5.9604645e-8;

        var flag: types.EnsureFloat(C) = 0;
        var h11: types.EnsureFloat(C) = 0;
        var h12: types.EnsureFloat(C) = 0;
        var h21: types.EnsureFloat(C) = 0;
        var h22: types.EnsureFloat(C) = 0;
        var p1: types.EnsureFloat(C) = 0;
        var p2: types.EnsureFloat(C) = 0;
        var q1: types.EnsureFloat(C) = 0;
        var q2: types.EnsureFloat(C) = 0;
        var temp: types.EnsureFloat(C) = 0;
        var u: types.EnsureFloat(C) = 0;

        if (d1.* < 0) {
            flag = -1;
            h11 = 0;
            h12 = 0;
            h21 = 0;
            h22 = 0;

            d1.* = 0;
            d2.* = 0;
            x1.* = 0;
        } else {
            ops.mul_( // p2 = d2 * y1
                &p2,
                d2.*,
                y1,
                ctx,
            ) catch unreachable;
            if (p2 == 0) {
                flag = -2;

                ops.set( // param[0] = flag
                    &param[0],
                    flag,
                    ctx,
                ) catch unreachable;

                return;
            }

            ops.mul_( // p1 = d1 * x1
                &p1,
                d1.*,
                x1.*,
                ctx,
            ) catch unreachable;

            ops.mul_( // q2 = p2 * y1
                &q2,
                p2,
                y1,
                ctx,
            ) catch unreachable;

            ops.mul_( // q1 = p1 * x1
                &q1,
                p1,
                x1.*,
                ctx,
            ) catch unreachable;

            if (ops.abs(q1, ctx) catch unreachable > ops.abs(q2, ctx) catch unreachable) {
                ops.div_( // h21 = -y1 / x1
                    &h21,
                    -y1,
                    x1.*,
                    ctx,
                ) catch unreachable;

                ops.div_( // h12 = p1 / p2
                    &h12,
                    p1,
                    p2,
                    ctx,
                ) catch unreachable;

                ops.sub_( // u = 1 - h12 * h21
                    &u,
                    1,
                    ops.mul(
                        h12,
                        h21,
                        ctx,
                    ) catch unreachable,
                    ctx,
                ) catch unreachable;

                if (u > 0) {
                    flag = 0;

                    ops.div_( // d1 /= u
                        d1,
                        d1.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // d2 /= u
                        d2,
                        d2.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.mul_( // x1 *= u
                        x1,
                        x1.*,
                        u,
                        ctx,
                    ) catch unreachable;
                } else {
                    flag = -1;
                    h11 = 0;
                    h12 = 0;
                    h21 = 0;
                    h22 = 0;

                    d1.* = 0;
                    d2.* = 0;
                    x1.* = 0;
                }
            } else {
                if (q2 < 0) {
                    flag = -1;
                    h11 = 0;
                    h12 = 0;
                    h21 = 0;
                    h22 = 0;

                    d1.* = 0;
                    d2.* = 0;
                    x1.* = 0;
                } else {
                    flag = 1;

                    ops.div_( // h11 = p1 / p2
                        &h11,
                        p1,
                        p2,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // h22 = x1 / y1
                        &h22,
                        x1.*,
                        y1,
                        ctx,
                    ) catch unreachable;

                    ops.add_( // u = 1 + h11 * h22
                        &u,
                        1,
                        ops.mul(
                            h11,
                            h22,
                            ctx,
                        ) catch unreachable,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // temp = d2 / u
                        &temp,
                        d2.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.div_( // d2 = d1 / u
                        d2,
                        d1.*,
                        u,
                        ctx,
                    ) catch unreachable;

                    ops.set( // d1 = temp
                        d1,
                        temp,
                        ctx,
                    ) catch unreachable;

                    ops.mul_( // x1 = y1 * u
                        x1,
                        y1,
                        u,
                        ctx,
                    ) catch unreachable;
                }
            }

            if (d1.* != 0) {
                while ((d1.* <= rgamsq) or (d1.* >= gamsq)) {
                    if (flag == 0) {
                        h11 = 1;
                        h22 = 1;
                        flag = -1;
                    } else {
                        h21 = -1;
                        h12 = 1;
                        flag = -1;
                    }
                    if (d1.* <= rgamsq) {
                        ops.mul_( // d1 *= gam * gam
                            d1,
                            d1.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // x1 /= gam
                            x1,
                            x1.*,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h11 /= gam
                            &h11,
                            h11,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h12 /= gam
                            &h12,
                            h12,
                            gam,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.div_( // d1 /= gam * gam
                            d1,
                            d1.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // x1 *= gam
                            x1,
                            x1.*,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h11 *= gam
                            &h11,
                            h11,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h12 *= gam
                            &h12,
                            h12,
                            gam,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }

            if (d2.* != 0) {
                while ((ops.abs(d2.*, ctx) catch unreachable <= rgamsq) or
                    (ops.abs(d2.*, ctx) catch unreachable >= gamsq))
                {
                    if (flag == 0) {
                        h11 = 1;
                        h22 = 1;
                        flag = -1;
                    } else {
                        h21 = -1;
                        h12 = 1;
                        flag = -1;
                    }
                    if (ops.abs(d2.*, ctx) catch unreachable <= rgamsq) {
                        ops.mul_( // d2 *= gam * gam
                            d2,
                            d2.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h21 /= gam
                            &h21,
                            h21,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.div_( // h22 /= gam
                            &h22,
                            h22,
                            gam,
                            ctx,
                        ) catch unreachable;
                    } else {
                        ops.div_( // d2 /= gam * gam
                            d2,
                            d2.*,
                            ops.mul(
                                gam,
                                gam,
                                ctx,
                            ) catch unreachable,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h21 *= gam
                            &h21,
                            h21,
                            gam,
                            ctx,
                        ) catch unreachable;

                        ops.mul_( // h22 *= gam
                            &h22,
                            h22,
                            gam,
                            ctx,
                        ) catch unreachable;
                    }
                }
            }
        }

        if (flag < 0) {
            ops.set( // param[1] = h11
                &param[1],
                h11,
                ctx,
            ) catch unreachable;

            ops.set( // param[2] = h21
                &param[2],
                h21,
                ctx,
            ) catch unreachable;

            ops.set( // param[3] = h12
                &param[3],
                h12,
                ctx,
            ) catch unreachable;

            ops.set( // param[4] = h22
                &param[4],
                h22,
                ctx,
            ) catch unreachable;
        } else if (flag == 0) {
            ops.set( // param[2] = h21
                &param[2],
                h21,
                ctx,
            ) catch unreachable;

            ops.set( // param[3] = h12
                &param[3],
                h12,
                ctx,
            ) catch unreachable;
        } else {
            ops.set( // param[1] = h11
                &param[1],
                h11,
                ctx,
            ) catch unreachable;

            ops.set( // param[4] = h22
                &param[4],
                h22,
                ctx,
            ) catch unreachable;
        }

        ops.set( // param[0] = flag
            &param[0],
            flag,
            ctx,
        ) catch unreachable;
    }
}
