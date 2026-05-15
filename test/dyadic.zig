test {
    // Override test flags
    const test_all = true;

    // Individual test flags
    const test_casts = true;
    const test_arithmetic = false;

    if (test_all or test_casts) {
        _ = @import("dyadic/initBool.zig");
        _ = @import("dyadic/initInt.zig");
        _ = @import("dyadic/initFloat.zig");

        _ = @import("dyadic/toFloat.zig");
    }

    if (test_all or test_arithmetic) {
        _ = @import("dyadic/add.zig");
        _ = @import("dyadic/sub.zig");
        _ = @import("dyadic/mul.zig");
        _ = @import("dyadic/fma.zig");
        _ = @import("dyadic/div.zig");
    }
}
