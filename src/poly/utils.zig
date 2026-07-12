const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

pub fn writeSuperscript(writer: *std.Io.Writer, value: usize) !void {
    const superscripts = [10][]const u8{ "⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹" };

    var buf: [32]u8 = undefined;
    const str = std.fmt.bufPrint(&buf, "{d}", .{value}) catch unreachable;

    for (str) |c| {
        const digit_index = c - '0';
        try writer.writeAll(superscripts[digit_index]);
    }
}

pub fn format(formatter: anytype, comptime num_fmt: []const u8, degree: usize, writer: *std.Io.Writer) !void {
    const gap_size: usize = 2;
    const gap_str = " " ** gap_size ++ "+" ++ " " ** gap_size;

    try writer.writeAll("          ");

    var printed_any_terms = false;

    var i: usize = degree + 1;
    while (i > 0) {
        i -= 1;

        const val = formatter.poly.get(i) catch unreachable;
        if (numeric.eq(val, numeric.zero(@TypeOf(val))))
            continue;

        if (printed_any_terms)
            try writer.writeAll(gap_str);

        try writer.print(num_fmt, .{val});

        if (i == 1) {
            try writer.writeAll(" x");
        } else if (i > 1) {
            try writer.writeAll(" x");
            try writeSuperscript(writer, i);
        }

        printed_any_terms = true;
    }

    if (!printed_any_terms)
        try writer.print(num_fmt, .{numeric.zero(meta.Numeric(meta.Child(@TypeOf(formatter.poly))))});
}
