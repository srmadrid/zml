const std = @import("std");
const zml = @import("zml");
//const ci = @import("c.zig");

pub fn main() !void {
    // const a: std.mem.Allocator = std.heap.page_allocator;
    var gpa = std.heap.DebugAllocator(.{}){};
    defer _ = gpa.deinit();
    const a = gpa.allocator();
    //_ = a;

    // try symbolicTesting(a);

    // try generalTesting(a);

    try addTesting(a);

    // try iterTesting(a);

    // try iterPerfTesting(a);

    // try multiIterTesting(a);

    // try perfTesting(a);

    // try typeTesting(a);

    // try transposeTesting(a);

    // try blasPerfTesting(a);

    // coreTesting();
}

fn ask_user(default: usize) !usize {
    const stdin = std.io.getStdIn().reader();
    const stdout = std.io.getStdOut().writer();

    var buf: [10]u8 = undefined;

    try stdout.print("Enter the number of iterations: ", .{});

    if (try (stdin.readUntilDelimiterOrEof(buf[0..], '\n'))) |user_input| {
        return std.fmt.parseInt(usize, user_input, 10);
    } else {
        return default;
    }
}

fn blasPerfTesting(a: std.mem.Allocator) !void {
    var A: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 10000, 10000 }, .{ .order = .ColumnMajor });
    defer A.deinit();

    // Ask for the number of iterations.
    const n = try ask_user(10);
    //const n = 10;
    std.debug.print("Number of iterations: {}\n", .{n});

    // Profiling
    var start_time: i128 = undefined;
    var end_time: i128 = undefined;
    var tot_time: i128 = 0;

    // Zig implementation 1
    std.debug.print("Zig implementation 1\n", .{});
    for (0..n) |_| {
        std.debug.print(".", .{});
        start_time = std.time.nanoTimestamp();
        std.mem.doNotOptimizeAway(zml.linalg.blas.asum(f64, @intCast(A.size), A.data.ptr, 1));
        //std.mem.doNotOptimizeAway(zml.BLAS.asum(f64, @divTrunc(@as(isize, @intCast(A.size)), @as(isize, 2)), A.data.ptr, 2));
        end_time = std.time.nanoTimestamp();

        tot_time += end_time - start_time;
    }
    std.debug.print("\n", .{});

    // Convert nanoseconds to seconds as a floating-point number.
    const duration_ns = tot_time;
    tot_time = 0;
    const duration_s: f128 = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Zig implementation 1 took: {d:.9} seconds\n", .{duration_s});

    // Zig implementation 2
    std.debug.print("Zig implementation 2\n", .{});
    for (0..n) |_| {
        std.debug.print(".", .{});
        start_time = std.time.nanoTimestamp();
        //std.mem.doNotOptimizeAway(zml.BLAS.asum2(f64, @intCast(A.size), A.data.ptr, 1));
        //std.mem.doNotOptimizeAway(zml.BLAS.asum2(f64, @divTrunc(@as(isize, @intCast(A.size)), @as(isize, 2)), A.data.ptr, 2));
        end_time = std.time.nanoTimestamp();

        tot_time += end_time - start_time;
    }
    std.debug.print("\n", .{});

    // Convert nanoseconds to seconds as a floating-point number.
    const duration_ns2 = tot_time;
    tot_time = 0;
    const duration_s2: f128 = @as(f128, @floatFromInt(duration_ns2)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Zig implementation 2 took: {d:.9} seconds\n", .{duration_s2});

    // Cblas implementation
    std.debug.print("Cblas implementation\n", .{});
    for (0..n) |_| {
        std.debug.print(".", .{});
        start_time = std.time.nanoTimestamp();
        //std.mem.doNotOptimizeAway(ci.cblas_dasum(@intCast(A.size), A.data.ptr, 1));
        //std.mem.doNotOptimizeAway(ci.cblas_dasum(@divTrunc(@as(c_int, @intCast(A.size)), @as(c_int, 2)), A.data.ptr, 2));
        end_time = std.time.nanoTimestamp();

        tot_time += end_time - start_time;
    }
    std.debug.print("\n", .{});

    // Convert nanoseconds to seconds as a floating-point number.
    const duration_ns3 = tot_time;
    tot_time = 0;
    const duration_s3: f128 = @as(f128, @floatFromInt(duration_ns3)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Cblas implementation took: {d:.9} seconds\n", .{duration_s3});
}

fn transposeTesting(a: std.mem.Allocator) !void {
    if (false) {
        var A: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{ .order = .RowMajor });
        defer A.deinit();
        const At = try A.transpose(null);

        std.debug.print("A dimentions = {}\n", .{A.shape.len});

        std.debug.print("A.shape = [  ", .{});
        for (A.shape[0..A.ndim]) |dim| {
            std.debug.print("{}  ", .{dim});
        }
        std.debug.print("]\n", .{});

        std.debug.print("A.strides = [  ", .{});
        for (A.strides[0..A.ndim]) |stride| {
            std.debug.print("{}  ", .{stride});
        }
        std.debug.print("]\n", .{});

        std.debug.print("A.size = {}\n", .{A.size});

        std.debug.print("At dimentions = {}\n", .{At.shape.len});

        std.debug.print("At.shape = [  ", .{});
        for (At.shape[0..At.ndim]) |dim| {
            std.debug.print("{}  ", .{dim});
        }
        std.debug.print("]\n", .{});

        std.debug.print("At.strides = [  ", .{});
        for (At.strides[0..At.ndim]) |stride| {
            std.debug.print("{}  ", .{stride});
        }
        std.debug.print("]\n", .{});

        std.debug.print("At.size = {}\n", .{At.size});

        for (0..A.size) |i| {
            A.data[i] = @floatFromInt(i + 1);
        }
        std.debug.print("A =\n", .{});
        for (0..A.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..A.shape[1]) |j| {
                std.debug.print("{!d:.2}  ", .{A.get(&[_]usize{ i, j })});
            }
            std.debug.print("\n", .{});
        }

        std.debug.print("At =\n", .{});
        for (0..At.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..At.shape[1]) |j| {
                std.debug.print("{!d:.2}  ", .{At.get(&[_]usize{ i, j })});
            }
            std.debug.print("\n", .{});
        }

        var B: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{ .order = .ColumnMajor });
        defer B.deinit();

        try B.add(A, try A.transpose(null));
        // or try B.add(A, At);

        std.debug.print("B =\n", .{});
        for (0..B.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..B.shape[1]) |j| {
                std.debug.print("{!d:.2}  ", .{B.get(&[_]usize{ i, j })});
            }
            std.debug.print("\n", .{});
        }

        const C = At.flatten();

        std.debug.print("C.shape = [  ", .{});
        for (C.shape[0..C.ndim]) |dim| {
            std.debug.print("{}  ", .{dim});
        }
        std.debug.print("]\n", .{});

        std.debug.print("C.size = {}\n", .{C.size});

        std.debug.print("C =\n", .{});
        for (0..C.size) |i| {
            std.debug.print("{!d:.2}  ", .{C.data[i]});
        }
        std.debug.print("\n", .{});

        var D: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{});
        defer D.deinit();
        D.setAll(1);
        var E: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 6, 6 }, .{});
        defer E.deinit();
        E.setAll(0);
    }

    var x = try zml.Array(f64).init(a, &.{ 5, 9, 4, 8 }, .{});

    std.debug.print("x.shape = [\n", .{});
    for (0..x.ndim) |i| {
        std.debug.print("x.shape[{}] = {}\n", .{ i, x.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x.strides = [\n", .{});
    for (0..x.ndim) |i| {
        std.debug.print("x.strides[{}] = {}\n", .{ i, x.metadata.dense.strides[i] });
    }
    std.debug.print("]\n\n", .{});

    for (0..x.shape[0]) |i| {
        for (0..x.shape[1]) |j| {
            for (0..x.shape[2]) |k| {
                for (0..x.shape[3]) |l| {
                    try x.set(&[_]usize{ i, j, k, l }, @floatFromInt(i * 1000 + j * 100 + k * 10 + l));
                }
            }
        }
    }

    const x1 = try x.slice(&.{ try .init(1, 5, 1), try .init(2, 5, 1), .single(2) });

    std.debug.print("x1 = x[1:5, 2:5, 2]\n", .{});
    std.debug.print("x1.shape = [\n", .{});
    for (0..x1.ndim) |i| {
        std.debug.print("x1.shape[{}] = {}\n", .{ i, x1.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x1.strides = [\n", .{});
    for (0..x1.ndim) |i| {
        std.debug.print("x1.strides[{}] = {}\n", .{ i, x1.metadata.strided.strides[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x1.offset = {}\n\n", .{x1.metadata.strided.offset});

    const x2 = try x.slice(&.{ .all, try .init(2, 9, 1) });

    std.debug.print("x2 = x[:, 2:10]\n", .{});
    std.debug.print("x2.shape = [\n", .{});
    for (0..x2.ndim) |i| {
        std.debug.print("x2.shape[{}] = {}\n", .{ i, x2.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x2.strides = [\n", .{});
    for (0..x2.ndim) |i| {
        std.debug.print("x2.strides[{}] = {}\n", .{ i, x2.metadata.strided.strides[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x2.offset = {}\n\n", .{x2.metadata.strided.offset});

    const x3 = try x.slice(&.{ try .init(1, 5, 1), try .init(8, 4, -3), .all_reverse, .all });

    std.debug.print("x3 = x[1:5, 8:4:-3, :, :]\n", .{});
    std.debug.print("x3.shape = [\n", .{});
    for (0..x3.ndim) |i| {
        std.debug.print("x3.shape[{}] = {}\n", .{ i, x3.shape[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x3.strides = [\n", .{});
    for (0..x3.ndim) |i| {
        std.debug.print("x3.strides[{}] = {}\n", .{ i, x3.metadata.strided.strides[i] });
    }
    std.debug.print("]\n", .{});
    std.debug.print("x3.offset = {}\n\n", .{x3.metadata.strided.offset});

    std.debug.print("x3 =\n", .{});
    for (0..x3.shape[2]) |k| {
        std.debug.print("Slice {}:\n", .{k});
        for (0..x3.shape[0]) |i| {
            std.debug.print("\t", .{});
            for (0..x3.shape[1]) |j| {
                std.debug.print("{d}  ", .{(try x3.get(&.{ i, j, k, 0 })).*});
            }
            std.debug.print("\n", .{});
        }
        std.debug.print("\n", .{});
    }
}

fn typeTesting(a: std.mem.Allocator) !void {
    const Complex = std.math.Complex;

    //const z = Complex(f64).init(10, 10);
    //const w = Complex(f64).init(20, 20);
    //var u: Complex(f64) = undefined;
    //const _add = @import("core/core.zig").supported._add;
    //_add(&u, z, w);
    //std.debug.print("u: {}\n", .{u});

    var B: zml.Array(Complex(f64)) = try zml.Array(Complex(f64)).init(a, &.{ 1, 1 }, .{ .order = .ColumnMajor });
    defer B.deinit();
    for (0..B.size) |i| {
        B.data[i].re = @floatFromInt(i + 1);
        B.data[i].im = @floatFromInt((i + 1) * 2);
    }
    std.debug.print("B =\n", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B.shape[1]) |j| {
            const z = try B.get(&.{ i, j });
            std.debug.print("{d:.2} + {d:.2}i  ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }

    var C: zml.Array(Complex(f64)) = try zml.Array(Complex(f64)).init(a, &.{ 5, 8 }, .{ .order = .RowMajor });
    defer C.deinit();
    for (0..C.size) |i| {
        C.data[i].re = @floatFromInt(i + 1);
        C.data[i].im = @floatFromInt((i + 1) * 2);
    }
    std.debug.print("\nC =\n", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C.shape[1]) |j| {
            const z = try C.get(&.{ i, j });
            std.debug.print("{d:.2} + {d:.2}i  ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }

    var D: zml.Array(Complex(f64)) = try zml.Array(Complex(f64)).init(a, &.{ 5, 8 }, .{ .order = .ColumnMajor });
    defer D.deinit();
    try D.add(B, C);
    // try zml.Array(u64).add(&D, B, C);
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[1]) |j| {
            const z = try D.get(&.{ i, j });
            std.debug.print("{d:.2} + {d:.2}i  ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }
}

fn symbolicTesting(a: std.mem.Allocator) !void {
    // Does not work, just for design.
    //const x = try zml.Variable.init(a, "x", &zml.Set.RealNumbers);
    //const S = try zml.Set.init.builder(a, "S", x, &[_]zml.Expression{zml.Expression.init.fromString(a, "sin(x) = 0", &[_]*zml.Symbol{&x})});
    //_ = S;
    const expr = "S = \\{x\\in\\mathbb{R}\\mid x > 0, \\arcsin(x) = 2\\pi k, \\forall k\\in\\mathbb{N}\\}";
    const arr = try zml.Expression.tokenize(a, expr);
    std.debug.print("{s}\n\nTokenized:\n", .{expr});
    for (arr.items) |token| {
        std.debug.print("<string = {s}, type = {}>\n", .{ token.string, token.type });
    }
}

fn generalTesting(a: std.mem.Allocator) !void {
    std.debug.print("Size of flags: {}\n", .{@sizeOf(zml.array.Flags)});

    var A: zml.Array(f64) = try .init(a, &.{ 20, 15, 8, 18 }, .{});
    defer A.deinit(a);

    std.debug.print("A dimentions = {}\n", .{A.ndim});

    std.debug.print("A.shape = [  ", .{});
    for (A.shape[0..A.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});

    std.debug.print("A.strides = [  ", .{});
    for (A.metadata.dense.strides[0..A.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});

    std.debug.print("A.size = {}\n", .{A.size});

    var B: zml.Array(f64) = try .linspace(a, 0, -1, 7, .{});
    defer B.deinit(a);

    std.debug.print("B = [", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("{d}, ", .{B.data[i]});
    }
    std.debug.print("]\n", .{});

    var C: zml.Array(f64) = try zml.abs(B, .{ .allocator = a });
    defer C.deinit(a);

    std.debug.print("C = [", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("{d}, ", .{C.data[i]});
    }
    std.debug.print("]\n", .{});

    var D: zml.Array(zml.cf64) = try zml.Array(zml.cf64).init(a, &.{5}, .{});
    defer D.deinit(a);
    for (0..D.size) |i| {
        D.data[i] = zml.cf64.init(@floatFromInt(i + 1), @floatFromInt((i + 1) * 2));
    }

    std.debug.print("D = [", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("{d} + {d}i, ", .{ D.data[i].re, D.data[i].im });
    }
    std.debug.print("]\n", .{});

    const D_broadcasted: zml.Array(zml.cf64) = try D.broadcast(&.{ 3, 5 });

    std.debug.print("D_broadcasted.shape = [  ", .{});
    for (D_broadcasted.shape[0..D_broadcasted.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("D_broadcasted = [\n", .{});
    for (0..D_broadcasted.shape[0]) |i| {
        for (0..D_broadcasted.shape[1]) |j| {
            const z = try D_broadcasted.get(&.{ i, j });
            std.debug.print("\t{d} + {d}i, ", .{ z.re, z.im });
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("]\n", .{});

    var E: zml.Array(f64) = try zml.abs(D, .{ .allocator = a });
    defer E.deinit(a);

    std.debug.print("E = [", .{});
    for (0..E.shape[0]) |i| {
        std.debug.print("{d}, ", .{E.data[i]});
    }
    std.debug.print("]\n", .{});

    var F = try E.slice(&.{try .init(null, null, -2)});
    std.debug.print("F = E[::-2]\n", .{});
    std.debug.print("F = [", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("{d}, ", .{(try F.get(&.{i})).*});
    }
    std.debug.print("]\n", .{});

    try zml.ceil_(&F, .{});

    std.debug.print("E after ceil = [", .{});
    for (0..E.shape[0]) |i| {
        std.debug.print("{d}, ", .{E.data[i]});
    }
    std.debug.print("]\n", .{});

    std.debug.print("F after ceil = [", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("{d}, ", .{(try F.get(&.{i})).*});
    }
    std.debug.print("]\n", .{});

    var G: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 5, 5 }, .{});
    defer G.deinit(a);

    for (0..G.shape[0]) |i| {
        for (0..G.shape[1]) |j| {
            try G.set(&.{ i, j }, -zml.scast(f64, (i + 1) * 10 + (j + 1)) / 100);
        }
    }

    std.debug.print("G =\n", .{});
    for (0..G.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try G.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var H: zml.Array(f64) = try G.slice(&.{ try .init(2, 5, 1), try .init(null, null, -2) });
    defer H.deinit(a);

    std.debug.print("H = G[2:5, ::-2]\n", .{});
    std.debug.print("H =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    try zml.abs_(&H, .{});

    std.debug.print("G after abs =\n", .{});
    for (0..G.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try G.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    std.debug.print("H after abs =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const shapes: []const []const usize = &.{
        &.{ 15, 3, 1 },
        &.{ 3, 1 },
        &.{ 1, 1, 5 },
    };
    _ = shapes;

    const G_T = try G.transpose(.{});
    std.debug.print("G_T = G.transpose()\n", .{});
    std.debug.print("G.shape = [  ", .{});
    for (G.shape[0..G.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G.strides = [  ", .{});
    for (G.metadata.dense.strides[0..G.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G_T.shape = [  ", .{});
    for (G_T.shape[0..G_T.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G_T.strides = [  ", .{});
    for (G_T.metadata.strided.strides[0..G_T.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});
    std.debug.print("G_T =\n", .{});
    for (0..G_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G_T.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try G_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const H_T = try H.transpose(.{});
    std.debug.print("H_T = H.transpose()\n", .{});
    std.debug.print("H_T =\n", .{});
    for (0..H_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H_T.shape[1]) |j| {
            std.debug.print("{d}  ", .{(try H_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }
}

fn addTesting(a: std.mem.Allocator) !void {
    var B: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 5, 1, 8 }, .{ .order = .rowMajor });
    defer B.deinit(a);
    for (0..B.size) |i| {
        B.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("B =\n", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B.shape[2]) |j| {
            std.debug.print("{!d:.2}  ", .{(try B.get(&.{ i, 0, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var C: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 1, 8 }, .{ .order = .rowMajor });
    defer C.deinit(a);
    for (0..C.size) |i| {
        C.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("\nC =\n", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try C.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var D: zml.Array(f64) = try zml.add(B, C, .{ .allocator = a });
    defer D.deinit(a);
    std.debug.print("\nD = B + C\n", .{});
    std.debug.print("D.shape = [  ", .{});
    for (D.shape[0..D.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[2]) |j| {
            std.debug.print("{!d:.2}  ", .{(try D.get(&.{ i, 0, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const B_T = try (try B.reshape(&.{ 5, 8 })).transpose(null);
    std.debug.print("B_T = B.transpose(null)\n", .{});
    std.debug.print("B_T =\n", .{});
    for (0..B_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B_T.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try B_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    const C_T = try C.transpose(null);
    std.debug.print("C_T = C.transpose(null)\n", .{});
    std.debug.print("C_T =\n", .{});
    for (0..C_T.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C_T.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try C_T.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var E: zml.Array(f64) = try zml.add(B_T, C_T, .{ .allocator = a });
    defer E.deinit(a);
    std.debug.print("\nE = B_T + C_T\n", .{});
    std.debug.print("\nE =\n", .{});
    for (0..E.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..E.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try E.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var F: zml.Array(f64) = try zml.add(5, E, .{ .allocator = a });
    defer F.deinit(a);
    std.debug.print("\nF = 5 + E\n", .{});
    std.debug.print("\nF =\n", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..F.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try F.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var G: zml.Array(f64) = try zml.array.apply2(a, B, C, zml.float.pow, .{});
    defer G.deinit(a);
    std.debug.print("\nG = B ** C\n", .{});
    std.debug.print("G =\n", .{});
    for (0..G.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..G.shape[2]) |j| {
            std.debug.print("{!d:.2}  ", .{(try G.get(&.{ i, 0, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var H: zml.Array(f64) = try F.slice(&.{ try .init(2, 5, 1), .all_reverse });
    std.debug.print("\nH = F[2:5, ::-1]\n", .{});
    std.debug.print("H =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    var I: zml.Array(f64) = try .init(a, &.{ H.shape[0], H.shape[1] }, .{});
    defer I.deinit(a);
    for (0..I.size) |i| {
        I.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("\nI =\n", .{});
    for (0..I.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..I.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try I.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    //try zml.div_(&H, I, .{});
    try zml.gamma_to(&H, I, .{});
    std.debug.print("\nH after add(I) =\n", .{});
    for (0..H.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..H.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try H.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }

    std.debug.print("\nF after add(I) =\n", .{});
    for (0..F.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..F.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{(try F.get(&.{ i, j })).*});
        }
        std.debug.print("\n", .{});
    }
}

fn iterTesting(a: std.mem.Allocator) !void {
    var iterArrR: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .rowMajor });
    defer iterArrR.deinit(a);
    var iterArrC: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .columnMajor });
    defer iterArrC.deinit(a);

    var iterR: zml.array.Iterator(f64) = zml.array.Iterator(f64).init(&iterArrR);
    var iterC: zml.array.Iterator(f64) = zml.array.Iterator(f64).init(&iterArrC);
    std.debug.print("Position(R)    Position(C)     R   C\n", .{});
    std.debug.print("----------------------\n", .{});
    std.debug.print("[  ", .{});
    for (0..iterR.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterR.position[i]});
    }
    std.debug.print("]  [  ", .{});
    for (0..iterC.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterC.position[i]});
    }
    std.debug.print("],  {},  {}\n", .{ iterR.index, iterC.index });
    while (iterR.nextAO(0, .leftToRight) != null and iterC.nextAO(0, .leftToRight) != null) {
        std.debug.print("[  ", .{});
        for (0..iterR.arr.ndim) |i| {
            std.debug.print("{}  ", .{iterR.position[i]});
        }
        std.debug.print("]  [  ", .{});
        for (0..iterC.arr.ndim) |i| {
            std.debug.print("{}  ", .{iterC.position[i]});
        }
        std.debug.print("],  {},  {}\n", .{ iterR.index, iterC.index });
    }
    _ = iterC.nextAO(0, .leftToRight); // Once R has reached the end it is null and the if never reaches C, so C has to run again
    std.debug.print("Final state:\n", .{});
    std.debug.print("[  ", .{});
    for (0..iterR.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterR.position[i]});
    }
    std.debug.print("]  [  ", .{});
    for (0..iterC.arr.ndim) |i| {
        std.debug.print("{}  ", .{iterC.position[i]});
    }
    std.debug.print("],  {},  {}\n\n", .{ iterR.index, iterC.index });
}

fn iterPerfTesting(a: std.mem.Allocator) !void {
    std.debug.print("Column major array, long next, release fast\n", .{});

    var arrBig: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 100, 100, 100, 100 }, .{ .order = .columnMajor });
    defer arrBig.deinit(a);
    var iterBig: zml.array.Iterator(f64) = zml.array.Iterator(f64).init(&arrBig);

    //
    const n: usize = 10;
    var start_time = std.time.nanoTimestamp();

    var count: u128 = 0;
    for (0..n) |_| {
        while (iterBig.nextAO(iterBig.array.ndim - 1, .rightToLeft) != null) {
            count += 1;
        }
    }

    var end_time = std.time.nanoTimestamp();
    var duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    var duration_s: f128 = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Row major iteration took: {d:.9} seconds\n", .{duration_s});

    start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        while (iterBig.nextAO(0, .leftToRight) != null) {
            count += 1;
        }
    }

    end_time = std.time.nanoTimestamp();
    duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Column major iteration took: {d:.9} seconds\n", .{duration_s});

    start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        while (iterBig.next() != null) {
            count += 1;
        }
    }

    end_time = std.time.nanoTimestamp();
    duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    duration_s = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("Default iteration took: {d:.9} seconds\n", .{duration_s});
}

fn multiIterTesting(a: std.mem.Allocator) !void {
    var arr1: zml.Array(f64) = try zml.Array(f64).init(a, &[_]usize{ 4, 2, 1 }, .{ .order = .RowMajor });
    defer arr1.deinit();
    var arr2: zml.Array(f64) = try zml.Array(f64).init(a, &[_]usize{ 2, 3 }, .{ .order = .ColumnMajor });
    defer arr2.deinit();

    // Other arrays for broadcasting testing.
    //var arr3: zml.Array(f64) = try zml.Array(f64).init(a, &[_]usize{ 2, 1, 1, 3 }, zml.array.Flags{ .RowMajorContiguous = true, .ColumnMajorContiguous = false });
    var arr3: zml.Array(f64) = try zml.Array(f64).init(a, &[_]usize{ 2, 4, 2, 3 }, .{ .order = .RowMajor });
    defer arr3.deinit();
    var scalar: zml.Array(f64) = try zml.Array(f64).init(a, &[_]usize{}, .{});
    defer scalar.deinit();

    var iter: zml.array.MultiIterator(f64) = try zml.array.MultiIterator(f64).init(&[_]zml.Array(f64){ arr1, arr2, arr3, scalar }, .{ .order = .ColumnMajor });

    std.debug.print("iter.shape = [  ", .{});
    for (0..iter.ndim) |i| {
        std.debug.print("{}  ", .{iter.shape[i]});
    }
    std.debug.print("]\n", .{});
    std.debug.print("arr1.shape = [  ", .{});
    for (0..arr1.ndim) |i| {
        std.debug.print("{}  ", .{arr1.shape[i]});
    }
    std.debug.print("]\n", .{});
    std.debug.print("arr2.shape = [  ", .{});
    for (0..arr2.ndim) |i| {
        std.debug.print("{}  ", .{arr2.shape[i]});
    }
    std.debug.print("]\n", .{});
    std.debug.print("arr3.shape = [  ", .{});
    for (0..arr3.ndim) |i| {
        std.debug.print("{}  ", .{arr3.shape[i]});
    }
    std.debug.print("]\n", .{});

    std.debug.print("F                      1                   2                3                      Scalar\n", .{});
    std.debug.print("-------------------------------------------------------------------------------------------\n", .{});
    std.debug.print("[  ", .{});
    for (0..iter.ndim) |i| {
        std.debug.print("{}  ", .{iter.position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.index});
    for (0..iter.iterators[0].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[0].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[0].index});
    for (0..iter.iterators[1].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[1].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[1].index});
    for (0..iter.iterators[2].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[2].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[2].index});
    for (0..iter.iterators[3].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[3].position[i]});
    }
    std.debug.print("]: {:0>2}\n", .{iter.iterators[3].index});
    while (iter.nextOrder(.RowMajor) != null) {
        std.debug.print("[  ", .{});
        for (0..iter.ndim) |i| {
            std.debug.print("{}  ", .{iter.position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.index});
        for (0..iter.iterators[0].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[0].position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.iterators[0].index});
        for (0..iter.iterators[1].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[1].position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.iterators[1].index});
        for (0..iter.iterators[2].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[2].position[i]});
        }
        std.debug.print("]: {:0>2}   [  ", .{iter.iterators[2].index});
        for (0..iter.iterators[3].ndim) |i| {
            std.debug.print("{}  ", .{iter.iterators[3].position[i]});
        }
        std.debug.print("]: {:0>2}\n", .{iter.iterators[3].index});
    }
    std.debug.print("Final state:\n", .{});
    std.debug.print("[  ", .{});
    for (0..iter.ndim) |i| {
        std.debug.print("{}  ", .{iter.position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.index});
    for (0..iter.iterators[0].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[0].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[0].index});
    for (0..iter.iterators[1].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[1].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[1].index});
    for (0..iter.iterators[2].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[2].position[i]});
    }
    std.debug.print("]: {:0>2}   [  ", .{iter.iterators[2].index});
    for (0..iter.iterators[3].ndim) |i| {
        std.debug.print("{}  ", .{iter.iterators[3].position[i]});
    }
    std.debug.print("]: {:0>2}\n", .{iter.iterators[3].index});
}

fn formatValueWithCustomPadding(value: i32, padding: u8) []const u8 {
    const buffer = std.mem.Allocator.buffer(@sizeOf(value), padding);
    const fmt = "{02d}";
    const len = std.fmt.bufPrint(buffer[0..], fmt, .{value}) catch return "";
    return buffer[0..len];
}

fn perfTesting(a: std.mem.Allocator) !void {
    // Perf testing
    var big1: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 10000, 10000 }, .{ .order = .ColumnMajor });
    defer big1.deinit();
    var big2: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 1, 1 }, .{ .order = .ColumnMajor });
    defer big2.deinit();
    var big3: zml.Array(f64) = try zml.Array(f64).init(a, &.{ 10000, 10000 }, .{ .order = .ColumnMajor });
    defer big3.deinit();

    std.debug.print("big3 dimentions = {}\n", .{big3.ndim});

    std.debug.print("big3.shape = [  ", .{});
    for (big3.shape[0..big3.ndim]) |dim| {
        std.debug.print("{}  ", .{dim});
    }
    std.debug.print("]\n", .{});

    std.debug.print("big3.strides = [  ", .{});
    for (big3.strides[0..big3.ndim]) |stride| {
        std.debug.print("{}  ", .{stride});
    }
    std.debug.print("]\n", .{});

    std.debug.print("big3.size = {}\n", .{big3.size});

    // Profiling
    const n: usize = 10;
    const start_time = std.time.nanoTimestamp();

    for (0..n) |_| {
        std.debug.print(".", .{});
        try @call(.auto, zml.Array(f64).add, .{ &big3, big1, big2 });
        //try big3.add(big1, big2);
    }
    std.debug.print("\n", .{});

    const end_time = std.time.nanoTimestamp();
    const duration_ns = end_time - start_time;

    // Convert nanoseconds to seconds as a floating-point number.
    const duration_s: f128 = @as(f128, @floatFromInt(duration_ns)) / (1_000_000_000.0 * @as(f128, @floatFromInt(n)));

    // Print the duration in seconds with high precision (e.g., 9 decimal places).
    std.debug.print("add took: {d:.9} seconds\n", .{duration_s});
}
