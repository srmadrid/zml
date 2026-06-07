const std = @import("std");

const zsl = @import("zsl");

const tzsl = @import("../zsl.zig");

const combinations = [_][2]type{
    // __stst
    .{ tzsl.vector.Static(zsl.cf64), tzsl.vector.Static(zsl.cf64) },

    // __stnu
    .{ tzsl.vector.Static(zsl.cf64), zsl.cf64 },

    // __nust
    .{ zsl.cf64, tzsl.vector.Static(zsl.cf64) },
};

const len_limits: [3]usize = .{
    7,
    16,
    33,
};

test "zsl.vector.apply2" {
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
        const B = if (comptime zsl.meta.isNumeric(combination[0])) tzsl.randomNumber(combination[0], rand) else try tzsl.vector.randomVector(
            combination[0],
            undefined,
            rand,
            len,
            incB,
            len / 4,
        );

        const C = if (comptime zsl.meta.isNumeric(combination[1])) tzsl.randomNumber(combination[1], rand) else try tzsl.vector.randomVector(
            combination[1],
            undefined,
            rand,
            len,
            incC,
            len / 4,
        );

        const A = if (comptime std.mem.eql(u8, op, "add"))
            zsl.vector.add(B, C)
        else if (comptime std.mem.eql(u8, op, "sub"))
            zsl.vector.sub(B, C)
        else if (comptime std.mem.eql(u8, op, "mul"))
            zsl.vector.mul(B, C)
        else
            zsl.vector.div(B, C);

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
