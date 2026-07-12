const std = @import("std");

pub fn hasNZ(v: anytype, row: usize, col: usize) bool {
    if (row > 0 and col > 0) return false;

    const target = if (row > 0) row else col;

    var lo: usize = 0;
    var hi: usize = v.nnz;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (v.idx[mid] == target) return true;
        if (v.idx[mid] < target) lo = mid + 1 else hi = mid;
    }
    return false;
}

pub fn format(formatter: anytype, comptime num_fmt: []const u8, len: usize, writer: *std.Io.Writer) !void {
    const gap_size: usize = 3;
    const gap_str = " " ** gap_size;

    var cell_width: usize = 0;
    for (0..len) |i| {
        const val = formatter.vector.get(i) catch unreachable;
        const val_len = std.fmt.count(num_fmt, .{val});
        if (val_len > cell_width)
            cell_width = val_len;
    }

    if (cell_width == 0)
        cell_width = 1;

    const inner_width = (len * cell_width) + ((len - 1) * gap_size);

    // Top border
    try writer.writeAll("          ┌");
    for (0..inner_width + 2) |_| try writer.writeAll(" ");
    try writer.writeAll("┐\n");

    // Data
    try writer.writeAll("          │ ");
    for (0..len) |i| {
        if (i > 0) try writer.writeAll(gap_str);
        const val = formatter.vector.get(i) catch unreachable;
        const val_len = std.fmt.count(num_fmt, .{val});

        for (0..(cell_width - val_len)) |_| try writer.writeAll(" ");
        try writer.print(num_fmt, .{val});
    }
    try writer.writeAll(" │\n");

    // Bot Border
    try writer.writeAll("          └");
    for (0..inner_width + 2) |_| try writer.writeAll(" ");
    try writer.writeAll("┘");
}
