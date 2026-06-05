const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("../zsl.zig");

const combinations = [_][3]type{
    // ststst
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // ststst: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // ststde
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // ststde: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // ststsp
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },
    // ststsp: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // ststnu
    .{ tzsl.vector.Static(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.cf64 },
    // ststnu: aliasing
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.cf64 },

    // stdest
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // stdest: aliasing
    .{ tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // stdede
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // stdesp
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // stdenu
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.cf64 },

    // stspst
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // stspst: aliasing
    .{ tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // stspde
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // stspsp
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // stspnu
    .{ tzsl.vector.Static(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // stnust
    .{ tzsl.vector.Static(zsl.cf32), zsl.cf64, tzsl.vector.Static(zsl.cf64) },
    // stnust: aliasing
    .{ tzsl.vector.Static(zsl.cf64), zsl.cf64, tzsl.vector.Static(zsl.cf64) },

    // stnude
    .{ tzsl.vector.Static(zsl.cf32), zsl.cf64, zsl.vector.Dense(zsl.cf64) },

    // stnusp
    .{ tzsl.vector.Static(zsl.cf32), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },

    // destst
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // destde
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // destde: aliasing
    .{ zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // destsp
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // destnu
    .{ zsl.vector.Dense(zsl.cf32), tzsl.vector.Static(zsl.cf64), zsl.cf64 },

    // dedest
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },
    // dedest: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // dedede
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // dedede: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // dedesp
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },
    // dedesp: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // dedenu
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Dense(zsl.cf64), zsl.cf64 },
    // dedenu: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64), zsl.cf64 },

    // despst
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // despde
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },
    // despde: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // despsp
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // despnu
    .{ zsl.vector.Dense(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // denust
    .{ zsl.vector.Dense(zsl.cf32), zsl.cf64, tzsl.vector.Static(zsl.cf64) },

    // denude
    .{ zsl.vector.Dense(zsl.cf32), zsl.cf64, zsl.vector.Dense(zsl.cf64) },
    // denude: aliasing
    .{ zsl.vector.Dense(zsl.cf64), zsl.cf64, zsl.vector.Dense(zsl.cf64) },

    // denusp
    .{ zsl.vector.Dense(zsl.cf32), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },

    // spspsp
    .{ zsl.vector.Sparse(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },
    // spspsp: aliasing
    .{ zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // spspnu
    .{ zsl.vector.Sparse(zsl.cf32), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },
    // spspnu: aliasing
    .{ zsl.vector.Sparse(zsl.cf64), zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // spnusp
    .{ zsl.vector.Sparse(zsl.cf32), zsl.cf64, zsl.vector.Sparse(zsl.cf64) },
    // spnusp: aliasing
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

test "zsl.vector.apply2_" {
    @setEvalBranchQuota(3000);

    const allocator = std.testing.allocator;

    var prng = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(std.testing.io).toSeconds()));
    const rand = prng.random();

    inline for (combinations) |fake_combination| {
        inline for (len_limits) |len| {
            comptime var combination = fake_combination;
            inline for (0..combination.len) |i| {
                if (@hasDecl(combination[i], "instantiate"))
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
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.add_)
        else if (comptime std.mem.eql(u8, op, "sub"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.sub_)
        else if (comptime std.mem.eql(u8, op, "mul"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.mul_)
        else
            try tzsl.vector.correctApply2(zsl.meta.Numeric(combination[0]), allocator, len, B, C, zsl.numeric.div_);
        defer D.deinit(allocator);

        const err = if (comptime std.mem.eql(u8, op, "add"))
            zsl.vector.add_(&A, B, C)
        else if (comptime std.mem.eql(u8, op, "sub"))
            zsl.vector.sub_(&A, B, C)
        else if (comptime std.mem.eql(u8, op, "mul"))
            zsl.vector.mul_(&A, B, C)
        else
            zsl.vector.div_(&A, B, C);

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
