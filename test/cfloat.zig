test {
    const test_ops = false;
    const tmp = false;

    if (test_ops) {
        _ = @import("cfloat/add.zig");
        _ = @import("cfloat/sub.zig");
        _ = @import("cfloat/mul.zig");
        _ = @import("cfloat/div.zig");
    }

    if (tmp) {
        _ = @import("cfloat/arg.zig");
        _ = @import("cfloat/abs.zig");
        _ = @import("cfloat/sqrt.zig");
        _ = @import("cfloat/exp.zig");
        _ = @import("cfloat/log.zig");
        _ = @import("cfloat/log10.zig");
        _ = @import("cfloat/sin.zig");
        _ = @import("cfloat/cos.zig");
        _ = @import("cfloat/tan.zig");
        _ = @import("cfloat/asin.zig");
        _ = @import("cfloat/acos.zig");
        _ = @import("cfloat/atan.zig");
        _ = @import("cfloat/sinh.zig");
        _ = @import("cfloat/cosh.zig");
        _ = @import("cfloat/tanh.zig");
        _ = @import("cfloat/asinh.zig");
        _ = @import("cfloat/acosh.zig");
        _ = @import("cfloat/atanh.zig");
    }

    _ = @import("cfloat/pow.zig");
}
