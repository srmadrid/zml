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
    }

    _ = @import("cfloat/log10.zig");
}
