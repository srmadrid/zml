test {
    const test_types = true;
    const test_float = true;
    const test_cfloat = true;
    const test_linalg = true;

    if (test_types) {
        _ = @import("types.zig");
    }

    if (test_float) {
        _ = @import("float.zig");
    }

    if (test_cfloat) {
        _ = @import("cfloat.zig");
    }

    if (test_linalg) {
        _ = @import("linalg.zig");
    }
}
