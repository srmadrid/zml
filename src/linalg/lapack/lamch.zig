const std = @import("std");

const types = @import("../../types.zig");
const scast = types.scast;
const Scalar = types.Scalar;
const ops = @import("../../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const linalg = @import("../../linalg.zig");
const blas = @import("../blas.zig");
const lapack = @import("../lapack.zig");
const Mach = lapack.Mach;

pub fn lamch(
    comptime T: type,
    cmach: Mach,
) Scalar(T) {
    if (comptime !types.isArbitraryPrecision(T)) {
        const rnd: T = constants.one(T, .{}) catch unreachable;

        if (comptime types.numericType(T) == .float or
            types.numericType(T) == .cfloat)
        {
            var eps: Scalar(T) = undefined;
            if (ops.eq(rnd, 1, .{}) catch unreachable) {
                eps = std.math.floatEps(Scalar(T)) * 0.5;
            } else {
                eps = std.math.floatEps(Scalar(T));
            }

            switch (cmach) {
                .eps => return eps,
                .sfmin => {
                    var sfmin: Scalar(T) = std.math.floatMin(Scalar(T));
                    const samll: Scalar(T) = 1 / std.math.floatMax(Scalar(T));
                    if (samll >= sfmin) {
                        sfmin = samll * (1 + eps);
                    }

                    return sfmin;
                },
                .base => return 2,
                .prec => return eps * 2,
                .t => return std.math.floatMantissaBits(Scalar(T)),
                .rnd => return rnd,
                .emin => return std.math.floatExponentMin(Scalar(T)),
                .rmin => return std.math.floatMin(Scalar(T)),
                .emax => return std.math.floatExponentMax(Scalar(T)),
                .rmax => return std.math.floatMax(Scalar(T)),
            }
        } else {
            switch (cmach) {
                .eps => return constants.one(T, .{}) catch unreachable,
                .sfmin => return if (comptime types.numericType(T) == .int and @typeInfo(T).int.signedness == .signed) {
                    int.minVal(T);
                } else {
                    constants.one(T, .{}) catch unreachable;
                },
                .base => return constants.two(T, .{}) catch unreachable,
                .prec => return constants.one(T, .{}) catch unreachable,
                .t => return scast(T, @bitSizeOf(T)),
                .rnd => return constants.one(T, .{}) catch unreachable,
                .emin => return scast(T, int.minVal(T)),
                .rmin => return if (comptime types.numericType(T) == .int and @typeInfo(T).int.signedness == .signed) {
                    int.minVal(T);
                } else {
                    constants.one(T, .{}) catch unreachable;
                },
                .emax => return scast(T, int.maxVal(T)),
                .rmax => return scast(T, int.maxVal(T)),
            }
        }
    } else {
        // Arbitrary precision types not supported yet
        @compileError("zml.linalg.lapack.lamch not implemented for arbitrary precision types yet");
    }
}
