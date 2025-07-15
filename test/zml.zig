test {
    const test_types = false;
    const test_float = false;
    const test_cfloat = false;
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
