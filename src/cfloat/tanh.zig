const types = @import("../types.zig");
const ops = @import("../ops.zig");

pub fn tanh(z: anytype) @TypeOf(z) {
    const Z = @TypeOf(z);

    comptime if (!types.isNumeric(Z) or types.numericType(Z) != .cfloat)
        @compileError("zml.cfloat.tanh: z must be a cfloat, got \n\tz: " ++ @typeName(Z) ++ "\n");

    const d: @TypeOf(z.re) = ops.add(
        ops.cosh(
            ops.mul(
                2,
                z.re,
                .{},
            ) catch unreachable,
            .{},
        ) catch unreachable,
        ops.cos(
            ops.mul(
                2,
                z.im,
                .{},
            ) catch unreachable,
            .{},
        ) catch unreachable,
        .{},
    ) catch unreachable;
    return .{
        .re = ops.div(
            ops.cosh(
                ops.mul(
                    2,
                    z.re,
                    .{},
                ) catch unreachable,
                .{},
            ) catch unreachable,
            d,
            .{},
        ) catch unreachable,
        .im = ops.div(
            ops.sinh(
                ops.mul(
                    2,
                    z.im,
                    .{},
                ) catch unreachable,
                .{},
            ) catch unreachable,
            d,
            .{},
        ) catch unreachable,
    };
}
