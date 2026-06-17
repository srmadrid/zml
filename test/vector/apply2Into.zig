const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("../zsl.zig");

const combinations = [_][3]type{
    // vecsta_vecsta_vecsta
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // vecsta_vecsta_vecsta: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecsta_vecsta_vecden
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // vecsta_vecsta_vecden: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecsta_vecsta_vecspa
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },
    // vecsta_vecsta_vecspa: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecsta_vecsta_num
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.cf64 },
    // vecsta_vecsta_num: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.cf64 },

    // vecsta_vecden_vecsta
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // vecsta_vecden_vecsta: aliasing
    .{ tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecsta_vecden_vecden
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecsta_vecden_vecspa
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecsta_vecden_num
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.cf64 },

    // vecsta_vecspa_vecsta
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // vecsta_vecspa_vecsta: aliasing
    .{ tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecsta_vecspa_vecden
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecsta_vecspa_vecspa
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecsta_vecspa_num
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // vecsta_num_vecsta
    .{ tzsl.vector.Static(zsl.cf32), zsl.cf64, tzsl.vector.Static(zsl.cf64) },
    // vecsta_num_vecsta: aliasing
    .{ tzsl.vector.Static(zsl.cf64), zsl.cf64, tzsl.vector.Static(zsl.cf64) },

    // vecsta_num_vecden
    .{ tzsl.vector.Static(zsl.cf32), zsl.cf64, zsl.vector.Dense(zsl.cf64) },

    // vecsta_num_vecspa
    .{ tzsl.vector.Static(zsl.cf32), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },

    // vecden_vecsta_vecsta
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecden_vecsta_vecden
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // vecden_vecsta_vecden: aliasing
    .{ zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecden_vecsta_vecspa
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecden_vecsta_num
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.cf64 },

    // vecden_vecden_vecsta
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // vecden_vecden_vecsta: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecden_vecden_vecden
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // vecden_vecden_vecden: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecden_vecden_vecspa
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },
    // vecden_vecden_vecspa: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecden_vecden_num
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.cf64 },
    // vecden_vecden_num: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), zsl.cf64 },

    // vecden_vecspa_vecsta
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecden_vecspa_vecden
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // vecden_vecspa_vecden: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecden_vecspa_vecspa
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecden_vecspa_num
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // vecden_num_vecsta
    .{ zsl.vector.Dense(zsl.cf32), zsl.cf64, tzsl.vector.Static(zsl.cf64) },

    // vecden_num_vecden
    .{ zsl.vector.Dense(zsl.cf32), zsl.cf64, zsl.vector.Dense(zsl.cf64) },
    // vecden_num_vecden: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.cf64, zsl.vector.Dense(zsl.cf64) },

    // vecden_num_vecspa
    .{ zsl.vector.Dense(zsl.cf32), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },

    // vecspa_vecspa_vecspa
    .{ zsl.vector.Sparse(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },
    // vecspa_vecspa_vecspa: aliasing
    .{ zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecspa_vecspa_num
    .{ zsl.vector.Sparse(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },
    // vecspa_vecspa_num: aliasing
    .{ zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // vecspa_num_vecspa
    .{ zsl.vector.Sparse(zsl.cf32), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },
    // vecspa_num_vecspa: aliasing
    .{ zsl.vector.Sparse(zsl.cf64), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },
};

const len_limits: [3]usize = .{
    7,
    16,
    33,
};

const inc_combinations = [_][3]isize{
    .{ 1, 1, 2 },
    .{ 1, 2, 1 },
    .{ 2, 1, 1 },
    .{ 1, -3, 2 },
    .{ -3, 1, 2 },
    .{ -3, 2, 1 },
    .{ 2, 1, 1 },
    .{ 5, -3, 2 },
    .{ 5, 2, -3 },
    .{ 2, -3, 5 },
    .{ 2, 5, -3 },
};

test "zsl.vector.apply2Into" {
    @setEvalBranchQuota(3 * combinations.len);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.testing.io).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |fake_combination| {
        inline for (len_limits) |len| {
            comptime var combination = fake_combination;
            inline for (0..combination.len) |i| {
                if (comptime !zsl.meta.isNumeric(combination[i]) and @hasDecl(combination[i], "instantiate"))
                    combination[i] = zsl.vector.Static(len, combination[i].Numeric);
            }

            const can_alias_B = combination[0] == combination[1];
            const can_alias_C = combination[0] == combination[2];

            const ops_to_test =
                if (comptime zsl.meta.isNumeric(combination[1]))
                    .{"mul"}
                else if (comptime zsl.meta.isNumeric(combination[2]))
                    .{ "mul", "div" }
                else
                    .{ "add", "sub" };

            try executeTestBlock(allocator, rand, combination, ops_to_test, len, 1, 1, 1, false, false);
            if (can_alias_B) try executeTestBlock(allocator, rand, combination, ops_to_test, len, 1, 1, 1, true, false);
            if (can_alias_C) try executeTestBlock(allocator, rand, combination, ops_to_test, len, 1, 1, 1, false, true);

            if (zsl.meta.isDenseVector(combination[0]) or zsl.meta.isDenseVector(combination[1]) or zsl.meta.isDenseVector(combination[2])) {
                for (inc_combinations) |incs| {
                    try executeTestBlock(allocator, rand, combination, ops_to_test, len, incs[0], incs[1], incs[2], false, false);
                    if (can_alias_B) try executeTestBlock(allocator, rand, combination, ops_to_test, len, incs[0], incs[1], incs[2], true, false);
                    if (can_alias_C) try executeTestBlock(allocator, rand, combination, ops_to_test, len, incs[0], incs[1], incs[2], false, true);
                }
            }
        }
    }
}

fn executeTestBlock(
    allocator: std.mem.Allocator,
    rand: std.Random,
    comptime combination: [3]type,
    comptime ops: anytype,
    comptime len: usize,
    incA: isize,
    incB: isize,
    incC: isize,
    comptime alias_B: bool,
    comptime alias_C: bool,
) !void {
    inline for (ops) |op| {
        var A = try tzsl.vector.randomVector(
            combination[0],
            allocator,
            rand,
            len,
            incA,
            len, // Ensure if A is sparse it always has space (A.nnz > B.nnz + C.nnz)
        );
        defer tzsl.deinit(allocator, &A);

        var B = if (comptime alias_B)
            A
        else if (comptime zsl.meta.isNumeric(combination[1]))
            tzsl.randomNumber(combination[1], rand)
        else
            try tzsl.vector.randomVector(
                combination[1],
                allocator,
                rand,
                len,
                incB,
                len / 4,
            );
        defer if (comptime !alias_B) tzsl.deinit(allocator, &B);

        var C = if (comptime alias_C)
            A
        else if (comptime zsl.meta.isNumeric(combination[2]))
            tzsl.randomNumber(combination[2], rand)
        else
            try tzsl.vector.randomVector(
                combination[2],
                allocator,
                rand,
                len,
                incC,
                len / 4,
            );
        defer if (comptime !alias_C) tzsl.deinit(allocator, &C);

        var D = if (comptime std.mem.eql(u8, op, "add"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.addInto)
        else if (comptime std.mem.eql(u8, op, "sub"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.subInto)
        else if (comptime std.mem.eql(u8, op, "mul"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.mulInto)
        else
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.divInto);
        defer D.deinit(allocator);

        const err = if (comptime std.mem.eql(u8, op, "add"))
            zsl.vector.addInto(&A, B, C)
        else if (comptime std.mem.eql(u8, op, "sub"))
            zsl.vector.subInto(&A, B, C)
        else if (comptime std.mem.eql(u8, op, "mul"))
            zsl.vector.mulInto(&A, B, C)
        else
            zsl.vector.divInto(&A, B, C);

        err catch |e| {
            std.debug.print(
                "Failed on A with dimension mistmatch: {s} = B: {s} {s} C: {s}, case len = {}, incA = {}, incB = {}, incC = {}\n",
                .{ @typeName(combination[0]), @typeName(combination[1]), op, @typeName(combination[2]), len, incA, incB, incC },
            );

            return e;
        };

        tzsl.vector.areEql(A, D) catch |e| {
            const aliasing = if (comptime alias_B) "B" else if (comptime alias_C) "C" else "no";
            std.debug.print(
                "Failed on A: {s} = B: {s} {s} C: {s}, case len = {}, incA = {}, incB = {}, incC = {}, aliasing = {s}\n",
                .{ @typeName(combination[0]), @typeName(combination[1]), op, @typeName(combination[2]), len, incA, incB, incC, aliasing },
            );

            tzsl.vector.printVector("A", A);
            if (comptime zsl.meta.isVector(@TypeOf(B))) tzsl.vector.printVector("B", B) else std.debug.print("B: {}\n", .{B});
            if (comptime zsl.meta.isVector(@TypeOf(C))) tzsl.vector.printVector("C", C) else std.debug.print("C: {}\n", .{C});
            tzsl.vector.printVector("D", D);

            return e;
        };
    }
}
