const std = @import("std");

const zml = @import("zml");

const float = @import("float.zig");

/// Prints a report of the test results for a function.
pub fn printReport(
    name: []const u8,
    data_cf32: anytype,
    results_cf32: anytype,
    data_cf64: anytype,
    results_cf64: anytype,
    data_cf80: anytype,
    results_cf80: anytype,
    data_cf128: anytype,
    results_cf128: anytype,
) void {
    std.debug.print("\n==========================================================\n", .{});
    std.debug.print(" Function: {s}\n", .{name});
    std.debug.print("==========================================================\n", .{});
    std.debug.print("| Type  | Max ULP | Mean ULP | 99%ile | % Exact | Status |\n", .{});
    std.debug.print("|-------|---------|----------|--------|---------|--------|\n", .{});

    // We process each precision block sequentially
    processPrecision(zml.cf32, data_cf32, results_cf32);
    processPrecision(zml.cf64, data_cf64, results_cf64);
    processPrecision(zml.cf80, data_cf80, results_cf80);
    processPrecision(zml.cf128, data_cf128, results_cf128);

    std.debug.print("----------------------------------------------------------\n", .{});
}

fn processPrecision(
    comptime T: type,
    data: anytype,
    results: anytype,
) void {
    const data_fields: comptime_int = @typeInfo(@TypeOf(data[0])).@"struct".fields.len;
    var ulps: [data.len * (if (zml.types.isComplex(@TypeOf(results[0]))) 2 else 1)]f128 = undefined;

    var sum_ulp: f128 = 0.0;
    var max_ulp: f128 = 0;
    var exact_count: u32 = 0;
    var subnormal_count: u32 = 0;

    if (comptime zml.types.isComplex(@TypeOf(results[0]))) {
        for (0..data.len) |i| {
            var ulp: f128 = undefined;
            if (float.isSubnormal(data[i][1].re) or // input 1 real part
                float.isSubnormal(data[i][1].im) or // input 1 imag part
                (data_fields == 3 and float.isSubnormal(data[i][2].re)) or // input 2 real part if exists
                (data_fields == 3 and float.isSubnormal(data[i][2].im)) or // input 2 imag part if exists
                float.isSubnormal(data[i][0].re) or // expected result real part
                float.isSubnormal(data[i][0].im) or // expected result imag part
                float.isSubnormal(results[i].re)) // actual result real part
            {
                ulp = 0.0;
                subnormal_count += 1;
            } else {
                ulp = float.ulpDistance(results[i].re, data[i][0].re);
            }

            if (ulp == 0) exact_count += 1;
            if (ulp > max_ulp) max_ulp = ulp;

            sum_ulp += ulp;
            ulps[2 * i] = ulp;

            if (float.isSubnormal(data[i][1].re) or // input 1 real part
                float.isSubnormal(data[i][1].im) or // input 1 imag part
                (data_fields == 3 and float.isSubnormal(data[i][2].re)) or // input 2 real part if exists
                (data_fields == 3 and float.isSubnormal(data[i][2].im)) or // input 2 imag part if exists
                float.isSubnormal(data[i][0].im) or // expected result imag part
                float.isSubnormal(results[i].im)) // actual result imag part
            {
                ulp = 0.0;
                subnormal_count += 1;
            } else {
                ulp = float.ulpDistance(results[i].im, data[i][0].im);
            }

            if (ulp == 0) exact_count += 1;
            if (ulp > max_ulp) max_ulp = ulp;

            sum_ulp += ulp;
            ulps[2 * i + 1] = ulp;
        }
    } else {
        for (0..data.len) |i| {
            var ulp: f128 = undefined;
            if (float.isSubnormal(data[i][1].re) or float.isSubnormal(data[i][1].im) or
                float.isSubnormal(data[i][0]) or float.isSubnormal(results[i]))
            {
                ulp = 0.0;
                subnormal_count += 1;
            } else {
                ulp = float.ulpDistance(results[i], data[i][0]);
            }

            if (ulp == 0) exact_count += 1;
            if (ulp > max_ulp) max_ulp = ulp;

            sum_ulp += ulp;
            ulps[i] = ulp;
        }
    }

    std.mem.sort(f128, &ulps, {}, std.sort.asc(f128));

    const mean_ulp: f128 = sum_ulp / zml.scast(f128, zml.int.sub(ulps.len, subnormal_count));
    const p99_index: u32 = zml.int.min(
        zml.scast(u32, std.math.ceil(zml.scast(f128, ulps.len) * 0.99)) - 1,
        ulps.len - 1,
    );
    const p99_ulp: f128 = ulps[p99_index];
    const exact_percentage: f128 = (zml.scast(f128, exact_count) / zml.scast(f128, ulps.len)) * 100.0;
    const status: []const u8 =
        if (max_ulp <= 2.0)
            "\x1b[32mPASS\x1b[0m"
        else if (max_ulp <= 4.0)
            "\x1b[33mWARN\x1b[0m"
        else
            "\x1b[31mFAIL\x1b[0m";

    std.debug.print(
        "| {s} |",
        .{
            switch (T) {
                zml.cf32 => "cf32 ",
                zml.cf64 => "cf64 ",
                zml.cf80 => "cf80 ",
                zml.cf128 => "cf128",
                else => unreachable,
            },
        },
    );

    if (max_ulp < 10)
        std.debug.print(" ", .{});

    std.debug.print(
        "      {d} |     {d:.2} |      {d} |",
        .{
            zml.scast(u128, max_ulp),
            mean_ulp,
            p99_ulp,
        },
    );

    if (exact_percentage < 100.0)
        std.debug.print(" ", .{});

    std.debug.print(" {d:.2}% |   {s} |\n", .{ exact_percentage, status });
}

test {
    // _ = @import("cfloat/add.zig");
    // _ = @import("cfloat/sub.zig");
    // _ = @import("cfloat/mul.zig");
    // _ = @import("cfloat/div.zig");

    _ = @import("cfloat/arg.zig");
    _ = @import("cfloat/abs.zig");
    _ = @import("cfloat/exp.zig");
    _ = @import("cfloat/log.zig");
    // _ = @import("cfloat/pow.zig");
    _ = @import("cfloat/sqrt.zig");
    _ = @import("cfloat/sin.zig");
    _ = @import("cfloat/cos.zig");
    _ = @import("cfloat/tan.zig");
    _ = @import("cfloat/asin.zig");
    _ = @import("cfloat/acos.zig");
    _ = @import("cfloat/atan.zig");
    // _ = @import("cfloat/sinh.zig");
    // _ = @import("cfloat/cosh.zig");
    // _ = @import("cfloat/tanh.zig");
    // _ = @import("cfloat/asinh.zig");
    // _ = @import("cfloat/acosh.zig");
    // _ = @import("cfloat/atanh.zig");
}
