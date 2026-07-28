const matrix = @import("../../matrix.zig");

pub fn index(layout: matrix.Layout, i: usize, j: usize, ld: usize) usize {
    return if (layout == .col_major) i + j * ld else i * ld + j;
}

pub fn col_ld(layout: matrix.Layout, ld: usize) usize {
    return if (layout == .col_major) 1 else ld;
}

pub fn row_ld(layout: matrix.Layout, ld: usize) usize {
    return if (layout == .col_major) ld else 1;
}
