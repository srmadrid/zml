const std = @import("std");
const zml = @import("zml.zig");

pub fn main() !void {
    // const a: std.mem.Allocator = std.heap.page_allocator;
    var gpa = std.heap.GeneralPurposeAllocator(.{}){};
    const a = gpa.allocator();

    // try symbolicTesting(a);

    // try generalTesting(a);

    // try addTesting(a);

    // try iterTesting(a);

    // try iterPerfTesting(a);

    // try multiIterTesting(a);

    // try perfTesting(a);

    // try typeTesting(a);

    try transposeTesting(a);
}

fn transposeTesting(a: std.mem.Allocator) !void {
    var A: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 6, 6 }, .{ .order = .RowMajor });
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

    var B: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 6, 6 }, .{ .order = .ColumnMajor });
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
}

fn typeTesting(a: std.mem.Allocator) !void {
    const Complex = std.math.Complex;

    //const z = Complex(f64).init(10, 10);
    //const w = Complex(f64).init(20, 20);
    //var u: Complex(f64) = undefined;
    //const _add = @import("core/core.zig").supported._add;
    //_add(&u, z, w);
    //std.debug.print("u: {}\n", .{u});

    var B: zml.NDArray(Complex(f64)) = try zml.NDArray(Complex(f64)).init(a, &.{ 1, 1 }, .{ .order = .ColumnMajor });
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

    var C: zml.NDArray(Complex(f64)) = try zml.NDArray(Complex(f64)).init(a, &.{ 5, 8 }, .{ .order = .RowMajor });
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

    var D: zml.NDArray(Complex(f64)) = try zml.NDArray(Complex(f64)).init(a, &.{ 5, 8 }, .{ .order = .ColumnMajor });
    defer D.deinit();
    try D.add(B, C);
    // try zml.NDArray(u64).add(&D, B, C);
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
    std.debug.print("Size of flags: {}\n", .{@sizeOf(zml.ndarray.Flags)});

    var A: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 20, 15, 8, 18 }, .{ .order = .ColumnMajor });
    defer A.deinit();

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

    //std.debug.print("Location of (2, 13, 3, 16) = {}\n", .{A._index(&[_]usize{ 2, 13, 3, 16 })});

    //var pos = [_]usize{ 0, 0, 0, 0 };
    //A._position(39562, &pos);
    //std.debug.print("Location of 39562 = [  ", .{});
    //for (pos) |p| {
    //    std.debug.print("{}  ", .{p});
    //}
    //std.debug.print("]\n\n", .{});
}

fn addTesting(a: std.mem.Allocator) !void {
    var B: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 1, 1 }, .{ .order = .ColumnMajor });
    defer B.deinit();
    for (0..B.size) |i| {
        B.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("B =\n", .{});
    for (0..B.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..B.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{B.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n", .{});
    }

    var C: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 5, 8 }, .{ .order = .RowMajor });
    defer C.deinit();
    for (0..C.size) |i| {
        C.data[i] = @floatFromInt(i + 1);
    }
    std.debug.print("\nC =\n", .{});
    for (0..C.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..C.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{C.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n", .{});
    }

    var D: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 5, 8 }, .{ .order = .ColumnMajor });
    defer D.deinit();
    try D.add(B, C);
    // try zml.NDArray(u64).add(&D, B, C);
    std.debug.print("\nD =\n", .{});
    for (0..D.shape[0]) |i| {
        std.debug.print("\t", .{});
        for (0..D.shape[1]) |j| {
            std.debug.print("{!d:.2}  ", .{D.get(&[_]usize{ i, j })});
        }
        std.debug.print("\n", .{});
    }
}

fn iterTesting(a: std.mem.Allocator) !void {
    var iterArrR: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .RowMajor });
    defer iterArrR.deinit();
    var iterArrC: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 3, 2, 4 }, .{ .order = .ColumnMajor });
    defer iterArrC.deinit();
    var iterR: zml.ndarray.Iterator(f64) = zml.ndarray.Iterator(f64).init(iterArrR);
    var iterC: zml.ndarray.Iterator(f64) = zml.ndarray.Iterator(f64).init(iterArrC);
    std.debug.print("Position(R)    Position(C)     R   C\n", .{});
    std.debug.print("----------------------\n", .{});
    std.debug.print("[  ", .{});
    for (0..iterR.ndim) |i| {
        std.debug.print("{}  ", .{iterR.position[i]});
    }
    std.debug.print("]  [  ", .{});
    for (0..iterC.ndim) |i| {
        std.debug.print("{}  ", .{iterC.position[i]});
    }
    std.debug.print("],  {},  {}\n", .{ iterR.index, iterC.index });
    while (iterR.nextOrder(.ColumnMajor) != null and iterC.nextOrder(.ColumnMajor) != null) {
        std.debug.print("[  ", .{});
        for (0..iterR.ndim) |i| {
            std.debug.print("{}  ", .{iterR.position[i]});
        }
        std.debug.print("]  [  ", .{});
        for (0..iterC.ndim) |i| {
            std.debug.print("{}  ", .{iterC.position[i]});
        }
        std.debug.print("],  {},  {}\n", .{ iterR.index, iterC.index });
    }
    _ = iterC.nextOrder(.ColumnMajor);
    std.debug.print("Final state:\n", .{});
    std.debug.print("[  ", .{});
    for (0..iterR.ndim) |i| {
        std.debug.print("{}  ", .{iterR.position[i]});
    }
    std.debug.print("]  [  ", .{});
    for (0..iterC.ndim) |i| {
        std.debug.print("{}  ", .{iterC.position[i]});
    }
    std.debug.print("],  {},  {}\n\n", .{ iterR.index, iterC.index });
}

fn iterPerfTesting(a: std.mem.Allocator) !void {
    std.debug.print("Column major array, long next, release fast\n", .{});

    var arrBig: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 100, 100, 100, 100 }, .{ .order = .ColumnMajor });
    defer arrBig.deinit();
    var iterBig: zml.ndarray.Iterator(f64) = zml.ndarray.Iterator(f64).init(arrBig);

    //
    const n: usize = 10;
    var start_time = std.time.nanoTimestamp();

    var count: u128 = 0;
    for (0..n) |_| {
        while (iterBig.nextOrder(.RowMajor) != null) {
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
        while (iterBig.nextOrder(.ColumnMajor) != null) {
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
    var arr1: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &[_]usize{ 4, 2, 1 }, .{ .order = .RowMajor });
    defer arr1.deinit();
    var arr2: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &[_]usize{ 2, 3 }, .{ .order = .ColumnMajor });
    defer arr2.deinit();

    // Other arrays for broadcasting testing.
    //var arr3: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &[_]usize{ 2, 1, 1, 3 }, zml.ndarray.Flags{ .RowMajorContiguous = true, .ColumnMajorContiguous = false });
    var arr3: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &[_]usize{ 2, 4, 2, 3 }, .{ .order = .RowMajor });
    defer arr3.deinit();
    var scalar: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &[_]usize{}, .{});
    defer scalar.deinit();

    var iter: zml.ndarray.MultiIterator(f64) = try zml.ndarray.MultiIterator(f64).init(&[_]zml.NDArray(f64){ arr1, arr2, arr3, scalar }, .{ .order = .ColumnMajor });

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
    var big1: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 10000, 10000 }, .{ .order = .ColumnMajor });
    defer big1.deinit();
    var big2: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 1, 1 }, .{ .order = .ColumnMajor });
    defer big2.deinit();
    var big3: zml.NDArray(f64) = try zml.NDArray(f64).init(a, &.{ 10000, 10000 }, .{ .order = .ColumnMajor });
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
        try @call(.auto, zml.NDArray(f64).add, .{ &big3, big1, big2 });
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
