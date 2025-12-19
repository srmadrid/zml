const std = @import("std");

const zml = @import("zml");

pub fn isSubnormal(x: anytype) bool {
    return !std.math.isNormal(x) and
        x != 0.0 and x != -0.0 and // zero is not subnormal
        !std.math.isInf(x) and // inf is not subnormal
        !std.math.isNan(x); // nan is not subnormal
}

pub fn ulpDistance(a: anytype, b: anytype) f128 {
    const ux: type = switch (@TypeOf(a)) {
        f16 => u16,
        f32 => u32,
        f64 => u64,
        f80 => u80,
        f128 => u128,
        else => unreachable,
    };

    var shift: u32 = switch (@TypeOf(a)) {
        f16 => 15,
        f32 => 31,
        f64 => 63,
        f80 => 79,
        f128 => 127,
        else => unreachable,
    };

    var ua: ux = @bitCast(a);
    var ub: ux = @bitCast(b);

    if (@TypeOf(a) == f80) { // Remove explicit integer bit
        const frac_mask: u80 = 0x7fffffffffffffff;

        const ua_high = (ua >> 64) << 63;
        const ub_high = (ub >> 64) << 63;

        ua = ua_high | (ua & frac_mask);
        ub = ub_high | (ub & frac_mask);

        shift = 78;
    }

    if (ua >> @intCast(shift) != 0)
        ua = ~ua
    else
        ua ^= @as(ux, 1) << @intCast(shift);

    if (ub >> @intCast(shift) != 0)
        ub = ~ub
    else
        ub ^= @as(ux, 1) << @intCast(shift);

    if (ua >= ub)
        return @floatFromInt(ua - ub)
    else
        return @floatFromInt(ub - ua);
}

/// Prints a report of the test results for a float function.
pub fn printReport(
    name: []const u8,
    data_f32: anytype,
    results_f32: anytype,
    data_f64: anytype,
    results_f64: anytype,
    data_f80: anytype,
    results_f80: anytype,
    data_f128: anytype,
    results_f128: anytype,
) void {
    std.debug.print("\n=========================================================\n", .{});
    std.debug.print(" Function: {s}\n", .{name});
    std.debug.print("=========================================================\n", .{});
    std.debug.print("| Type | Max ULP | Mean ULP | 99%ile | % Exact | Status |\n", .{});
    std.debug.print("|------|---------|----------|--------|---------|--------|\n", .{});

    // We process each precision block sequentially
    processPrecision(f32, data_f32, results_f32);
    processPrecision(f64, data_f64, results_f64);
    processPrecision(f80, data_f80, results_f80);
    processPrecision(f128, data_f128, results_f128);

    std.debug.print("---------------------------------------------------------\n", .{});
}

fn processPrecision(
    comptime T: type,
    data: anytype,
    results: anytype,
) void {
    const data_fields: comptime_int = @typeInfo(@TypeOf(data[0])).@"struct".fields.len;
    var ulps: [data.len]f128 = undefined;

    var sum_ulp: f128 = 0.0;
    var max_ulp: f128 = 0;
    var exact_count: u32 = 0;
    var subnormal_count: u32 = 0;

    for (0..data.len) |i| {
        var ulp: f128 = undefined;
        if (isSubnormal(data[i][1]) or // input 1
            (data_fields == 3 and isSubnormal(data[i][2])) or // input 2 if exists
            isSubnormal(data[i][0]) or // expected result
            isSubnormal(results[i])) // actual result
        {
            ulp = 0.0;
            subnormal_count += 1;
        } else {
            ulp = ulpDistance(results[i], data[i][0]);
        }

        if (ulp == 0) exact_count += 1;
        if (ulp > max_ulp) max_ulp = ulp;

        sum_ulp += ulp;
        ulps[i] = ulp;
    }

    std.mem.sort(f128, &ulps, {}, std.sort.asc(f128));

    const mean_ulp: f128 = sum_ulp / zml.scast(f128, zml.int.sub(data.len, subnormal_count));
    const p99_index: u32 = zml.int.min(
        zml.scast(u32, std.math.ceil(zml.scast(f128, data.len) * 0.99)) - 1,
        data.len - 1,
    );
    const p99_ulp: f128 = ulps[p99_index];
    const exact_percentage: f128 = (zml.scast(f128, exact_count) / zml.scast(f128, data.len)) * 100.0;
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
            if (T == f128) @typeName(T) else @typeName(T) ++ " ",
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
    const test_basic = false;
    const test_exponential = true;
    const test_power = true;
    const test_trigonometric = true;
    const test_hyperbolic = true;
    const test_error_gamma = true;
    const test_bessel = false;
    const test_nearest_integer = false;
    const test_floating_point = false;

    if (test_basic) {
        //_ = @import("float/fmod.zig");
        //_ = @import("float/remainder.zig");
        //_ = @import("float/remquo.zig");
        //_ = @import("float/fdim.zig");
    }

    if (test_exponential) {
        _ = @import("float/exp.zig");
        _ = @import("float/exp2.zig");
        _ = @import("float/expm1.zig");
        _ = @import("float/log.zig");
        _ = @import("float/log10.zig");
        _ = @import("float/log2.zig");
        _ = @import("float/log1p.zig");
    }

    if (test_power) {
        _ = @import("float/pow.zig");
        _ = @import("float/sqrt.zig");
        _ = @import("float/cbrt.zig");
        _ = @import("float/hypot.zig");
    }

    if (test_trigonometric) {
        _ = @import("float/sin.zig");
        _ = @import("float/cos.zig");
        _ = @import("float/tan.zig");
        _ = @import("float/asin.zig");
        _ = @import("float/acos.zig");
        _ = @import("float/atan.zig");
        _ = @import("float/atan2.zig");
        _ = @import("float/sincos.zig");
    }

    if (test_hyperbolic) {
        _ = @import("float/sinh.zig");
        _ = @import("float/cosh.zig");
        _ = @import("float/tanh.zig");
        _ = @import("float/asinh.zig");
        _ = @import("float/acosh.zig");
        _ = @import("float/atanh.zig");
    }

    if (test_error_gamma) {
        _ = @import("float/erf.zig");
        _ = @import("float/erfc.zig");
        _ = @import("float/gamma.zig");
        _ = @import("float/lgamma.zig");
    }

    if (test_bessel) {
        //_ = @import("float/j0.zig");
        //_ = @import("float/j1.zig");
        //_ = @import("float/jn.zig");
        //_ = @import("float/y0.zig");
        //_ = @import("float/y1.zig");
        //_ = @import("float/yn.zig");
    }

    if (test_nearest_integer) {
        //_ = @import("float/nearbyint.zig");
        //_ = @import("float/lrint.zig");
        //_ = @import("float/llrint.zig");
    }

    if (test_floating_point) {
        //_ = @import("float/ilogb.zig");
        //_ = @import("float/logb.zig");
        //_ = @import("float/nextafter.zig");
        //_ = @import("float/nexttoward.zig");
    }
}
