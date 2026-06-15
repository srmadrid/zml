const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("../zsl.zig");

const combinations = [_][2]type{
    // vecsta_vecsta
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecsta_vecden
    .{ tzsl.vector.Static(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecsta_vecspa
    .{ tzsl.vector.Static(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecsta_num
    .{ tzsl.vector.Static(zsl.cf64), zsl.cf64 },

    // vecden_vecsta
    .{ zsl.vector.Dense(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecden_vecden
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecden_vecspa
    .{ zsl.vector.Dense(zsl.cf64), zsl.vector.Sparse(zsl.cf64) },

    // vecden_num
    .{ zsl.vector.Dense(zsl.cf64), zsl.cf64 },

    // vecspa_vecsta
    .{ zsl.vector.Sparse(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // vecspa_vecden
    .{ zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecspa_vecspa
    .{ zsl.vector.Sparse(zsl.cf64), zsl.vector.Dense(zsl.cf64) },

    // vecspa_num
    .{ zsl.vector.Sparse(zsl.cf64), zsl.cf64 },

    // num_vecsta
    .{ zsl.cf64, tzsl.vector.Static(zsl.cf64) },

    // num_vecden
    .{ zsl.cf64, zsl.vector.Dense(zsl.cf64) },

    // num_vecspa
    .{ zsl.cf64, zsl.vector.Sparse(zsl.cf64) },
};

const len_limits: [3]usize = .{
    7,
    16,
    33,
};

const inc_combinations = [_][2]isize{
    .{ 1, 2 },
    .{ 2, 1 },
    .{ -3, 2 },
    .{ 2, -3 },
    .{ -3, 5 },
    .{ 5, -3 },
};

test "zsl.vector.apply2Alloc" {
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

            const ops_to_test =
                if (comptime zsl.meta.isNumeric(combination[0]))
                    .{"mul"}
                else if (comptime zsl.meta.isNumeric(combination[1]))
                    .{ "mul", "div" }
                else
                    .{ "add", "sub" };

            try executeTestBlock(allocator, rand, combination, ops_to_test, len, 1, 1);

            if (zsl.meta.isDenseVector(combination[0]) or zsl.meta.isDenseVector(combination[1])) {
                for (inc_combinations) |incs| {
                    try executeTestBlock(allocator, rand, combination, ops_to_test, len, incs[0], incs[1]);
                }
            }
        }
    }
}

fn executeTestBlock(
    allocator: std.mem.Allocator,
    rand: std.Random,
    comptime combination: [2]type,
    comptime ops: anytype,
    comptime len: usize,
    incB: isize,
    incC: isize,
) !void {
    inline for (ops) |op| {
        var B = if (comptime zsl.meta.isNumeric(combination[0])) tzsl.randomNumber(combination[0], rand) else try tzsl.vector.randomVector(
            combination[0],
            allocator,
            rand,
            len,
            incB,
            len / 4,
        );
        defer tzsl.deinit(allocator, &B);

        var C = if (comptime zsl.meta.isNumeric(combination[1])) tzsl.randomNumber(combination[1], rand) else try tzsl.vector.randomVector(
            combination[1],
            allocator,
            rand,
            len,
            incC,
            len / 4,
        );
        defer tzsl.deinit(allocator, &C);

        var A = if (comptime std.mem.eql(u8, op, "add"))
            try zsl.vector.addAlloc(allocator, B, C)
        else if (comptime std.mem.eql(u8, op, "sub"))
            try zsl.vector.subAlloc(allocator, B, C)
        else if (comptime std.mem.eql(u8, op, "mul"))
            try zsl.vector.mulAlloc(allocator, B, C)
        else
            try zsl.vector.divAlloc(allocator, B, C);
        defer tzsl.deinit(allocator, &A);

        var D = if (comptime std.mem.eql(u8, op, "add"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(@TypeOf(A)), allocator, len, B, C, zsl.numeric.addInto)
        else if (comptime std.mem.eql(u8, op, "sub"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(@TypeOf(A)), allocator, len, B, C, zsl.numeric.subInto)
        else if (comptime std.mem.eql(u8, op, "mul"))
            try tzsl.vector.correctApply2(zsl.meta.Numeric(@TypeOf(A)), allocator, len, B, C, zsl.numeric.mulInto)
        else
            try tzsl.vector.correctApply2(zsl.meta.Numeric(@TypeOf(A)), allocator, len, B, C, zsl.numeric.divInto);
        defer D.deinit(allocator);

        tzsl.vector.areEql(A, D) catch |e| {
            std.debug.print(
                "Failed on A: {s} = B: {s} {s} C: {s}, case len = {}\n",
                .{ @typeName(@TypeOf(A)), @typeName(combination[0]), op, @typeName(combination[1]), len },
            );
            return e;
        };
    }
}
