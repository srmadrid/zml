const std = @import("std");
const zsl = @import("zsl");

pub fn main(init: std.process.Init) !void {
    @setEvalBranchQuota(10_000);

    // try blas_lv1_threshold_calibration(init);
    // try blas_lv2_threshold_calibration(init);

    // const arena = init.arena.allocator();
    const gpa = init.gpa;

    // const io = init.io;

    // var xoshiro = std.Random.DefaultPrng.init(@bitCast(std.Io.Clock.real.now(io).toMicroseconds()));
    // const prng = xoshiro.random();
    // const normal = zsl.stats.Normal(f64).init(0.0, 1.0);

    var tape: zsl.autodiff.Tape(f64) = try .init(gpa, 1000);
    defer tape.deinit(gpa);

    const f = struct {
        pub fn f(x: anytype, y: anytype) @TypeOf(x) {
            return zsl.numeric.mul(x, zsl.numeric.sqrt(zsl.numeric.add(x, y)));
        }
    }.f;

    const x: zsl.autodiff.Var(f64) = .init(&tape, 3.0);
    const y: zsl.autodiff.Var(f64) = .init(&tape, 2.0);

    const r = f(x, y);

    std.debug.print("f(3.0, 2.0) = {}\n", .{r.val()});

    r.backward();

    std.debug.print("df/dx = {}\n", .{x.grad()});
    std.debug.print("df/dy = {}\n", .{y.grad()});

    // const m = 4;
    // const n = 4;

    // var a: zsl.matrix.general.Dense(f64, .col_major) = try .initFn(gpa, m, n, zsl.stats.Normal(f64).sample, .{ normal, prng });
    // defer a.deinit(gpa);

    // var a_clone = try a.clone(gpa);
    // defer a_clone.deinit(gpa);

    // var s: zsl.matrix.diagonal.Sparse(f64) = try .init(gpa, m, n);
    // defer s.deinit(gpa);
    // var u: zsl.matrix.general.Dense(f64, .col_major) = try .init(gpa, m, m);
    // defer u.deinit(gpa);
    // var vt: zsl.matrix.general.Dense(f64, .col_major) = try .init(gpa, n, n);
    // defer vt.deinit(gpa);

    // const superb_len = zsl.int.max(1, zsl.int.min(m, n));
    // const superb = try gpa.alloc(f64, superb_len);
    // defer gpa.free(superb);

    // const info = zsl.linalg.lapacke.dgesvd(
    //     zsl.linalg.cblas.layout.col_major,
    //     zsl.linalg.lapacke.job.all,
    //     zsl.linalg.lapacke.job.all,
    //     m,
    //     n,
    //     a_clone.data,
    //     zsl.numeric.cast(isize, a_clone.ld),
    //     s.data,
    //     u.data,
    //     zsl.numeric.cast(isize, u.ld),
    //     vt.data,
    //     zsl.numeric.cast(isize, vt.ld),
    //     superb.ptr,
    // );

    // if (info != 0)
    //     return error.SVD;

    // var us = try zsl.linalg.matmulAlloc(gpa, u, s);
    // defer us.deinit(gpa);
    // var a_reconstructed = try zsl.linalg.matmulAlloc(gpa, us, vt);
    // defer a_reconstructed.deinit(gpa);

    // printMatrix("A = U S V^T", a);
    // printMatrix("U", u);
    // printMatrix("S", s);
    // printMatrix("V^T", vt);
    // printMatrix("A = U S V^T (reconstructed)", a_reconstructed);

    // var diff = try zsl.matrix.subAlloc(gpa, a_reconstructed, a);
    // defer diff.deinit(gpa);

    // const diff_norm = try zsl.linalg.normAlloc(gpa, diff, .frobenius);

    // std.debug.print("‖A - U S V^T‖ = {e}\n", .{diff_norm});
}

pub fn blas_lv1_threshold_calibration(init: std.process.Init) !void {
    @setEvalBranchQuota(10000);
    const io = init.io;
    var gpa: std.heap.DebugAllocator(.{}) = .{ .backing_allocator = init.gpa };
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const max_n: usize = 268_435_457;
    // const alpha: f64 = 2.0;
    const a = try allocator.alloc(f64, zsl.numeric.cast(usize, max_n));
    const b = try allocator.alloc(f64, zsl.numeric.cast(usize, max_n));
    defer allocator.free(a);
    defer allocator.free(b);
    for (a, b, 0..) |*a_val, *b_val, i| {
        a_val.* = @as(f64, @floatFromInt(i % 100)) / 100.0;
        b_val.* = @as(f64, @floatFromInt(i % 33)) / 33.0;
    }

    const thresholds = [_]usize{
        16_384,
        32_768,
        65_536,
        131_072,
        262_144,
        524_288,
        1_048_576,
        2_097_152,
        4_194_304,
        8_388_608,
        16_777_216,
    };

    const options_max_threads = 64; // options.max_threads

    const hw_threads: usize = std.Thread.getCpuCount() catch 1;
    const effective_cap: usize = @min(hw_threads, options_max_threads);

    // Count how many N values we'll iterate over, so we can allocate per-cell
    // ratio storage for the summary pass.
    var n_rows: usize = 0;
    {
        var x: usize = 128;
        while (x <= max_n) : (x *= 2) n_rows += 1;
    }

    // ratios[k][row] = parallel_time / serial_time for threshold k at row `row`.
    // thread_counts[k][row] = how many threads actually spawned there.
    const ratios = try allocator.alloc([]f64, thresholds.len);
    defer allocator.free(ratios);
    const thread_counts = try allocator.alloc([]usize, thresholds.len);
    defer allocator.free(thread_counts);
    const row_ns = try allocator.alloc(usize, n_rows);
    defer allocator.free(row_ns);
    for (ratios, thread_counts) |*r, *t| {
        r.* = try allocator.alloc(f64, n_rows);
        t.* = try allocator.alloc(usize, n_rows);
    }
    defer for (ratios, thread_counts) |r, t| {
        allocator.free(r);
        allocator.free(t);
    };

    std.debug.print("\n=== IAMIN Threshold Calibration ===\n", .{});
    std.debug.print("Hardware threads: {d}, options.max_threads: {d}, effective cap: {d}\n", .{ hw_threads, options_max_threads, effective_cap });
    std.debug.print("Automatic-mode rule: threads = max(1, min(effective_cap, n / T))\n\n", .{});
    std.debug.print("Cell format: 'R.RRx(Nt)' -> R = parallel_time/serial_time, N = threads spawned\n", .{});
    std.debug.print("  R < 1.00 -> parallel wins            R > 1.00 -> parallel loses\n", .{});
    std.debug.print("  (1t)     -> n/T < 2, stayed serial   higher T -> fewer threads on large N\n\n", .{});

    // Header
    std.debug.print("{s:>11} | {s:>12}", .{ "N", "Serial (ns)" });
    for (thresholds) |th| {
        var buf: [16]u8 = undefined;
        std.debug.print(" | {s:>12}", .{fmtT(&buf, th)});
    }
    std.debug.print(" | {s:>22}\n", .{"Best for this N"});

    const sep_len: usize = 11 + 3 + 12 + thresholds.len * 15 + 3 + 22;
    var i: usize = 0;
    while (i < sep_len) : (i += 1) std.debug.print("-", .{});
    std.debug.print("\n", .{});

    var row: usize = 0;
    var current_n: usize = 128;
    while (current_n <= max_n) : (current_n *= 2) {
        row_ns[row] = current_n;

        const iters: usize =
            if (current_n < 1_000_000) 100 else if (current_n < 16_000_000) 50 else if (current_n < 256_000_000) 10 else 5;

        // Serial baseline
        std.mem.doNotOptimizeAway(
            try zsl.linalg.blas.iamin(current_n, a.ptr, 1, .{ .num_threads = 1 }),
        );
        const s0 = std.Io.Clock.real.now(io);
        for (0..iters) |_| {
            std.mem.doNotOptimizeAway(
                try zsl.linalg.blas.iamin(current_n, a.ptr, 1, .{ .num_threads = 1 }),
            );
        }
        const s1 = std.Io.Clock.real.now(io);
        const serial_ns = (zsl.numeric.cast(f128, s1.toNanoseconds()) - zsl.numeric.cast(f128, s0.toNanoseconds())) / @as(f128, iters);

        std.debug.print("{d:>11} | {d:>12.0}", .{ current_n, @as(f64, @floatCast(serial_ns)) });

        var best_ratio: f128 = 1.0;
        var best_idx: ?usize = null;
        var best_threads: usize = 1;

        for (thresholds, 0..) |T, k| {
            std.mem.doNotOptimizeAway(
                try zsl.linalg.blas.iamin(current_n, a.ptr, 1, .{ .num_threads = 0, .parallel_threshold = T }),
            );
            const p0 = std.Io.Clock.real.now(io);
            for (0..iters) |_| {
                std.mem.doNotOptimizeAway(
                    try zsl.linalg.blas.iamin(current_n, a.ptr, 1, .{ .num_threads = 0, .parallel_threshold = T }),
                );
            }
            const p1 = std.Io.Clock.real.now(io);
            const ns = (zsl.numeric.cast(f128, p1.toNanoseconds()) - zsl.numeric.cast(f128, p0.toNanoseconds())) / @as(f128, iters);
            const ratio: f128 = ns / serial_ns;

            const n_usize = zsl.numeric.cast(usize, current_n);
            const threads = @max(@as(usize, 1), @min(n_usize / T, effective_cap));

            ratios[k][row] = @as(f64, @floatCast(ratio));
            thread_counts[k][row] = threads;

            std.debug.print(" | {d:>6.2}x({d:>2}t)", .{
                @as(f64, @floatCast(ratio)), threads,
            });

            if (ratio < best_ratio) {
                best_ratio = ratio;
                best_idx = k;
                best_threads = threads;
            }
        }

        if (best_idx) |k| {
            var tbuf: [16]u8 = undefined;
            var sbuf: [48]u8 = undefined;
            const tl = fmtT(&tbuf, thresholds[k]);
            const sum = std.fmt.bufPrint(&sbuf, "{s} @ {d}t ({d:.2}x)", .{ tl, best_threads, @as(f64, @floatCast(best_ratio)) }) catch "?";
            std.debug.print(" | {s:>22}\n", .{sum});
        } else {
            std.debug.print(" | {s:>22}\n", .{"serial (no parallel won)"});
        }

        row += 1;
    }

    // Summary metrics
    //
    // Noise band: ratios within [1 - NOISE, 1 + NOISE] are treated as "tied
    // with serial" — neither a win nor a regression. Calibrated from the
    // small-N rows where every threshold is degenerately serial yet still
    // shows ~5% wobble.
    const NOISE: f64 = 0.05;
    // Regression threshold: ratios above this are counted as real losses.
    const REGRESSION: f64 = 1.10;

    std.debug.print("\n=== Summary (noise band = ±{d:.0}%, regression = >{d:.0}%) ===\n\n", .{
        NOISE * 100, (REGRESSION - 1.0) * 100,
    });

    std.debug.print("{s:>10} | {s:>9} | {s:>9} | {s:>9} | {s:>12} | {s:>9} | {s:>14}\n", .{
        "Threshold", "Geomean*", "Best", "Worst", "Worst @ N", "Reg. N", "First win @ N",
    });
    var j: usize = 0;
    while (j < 10 + 3 + 9 + 3 + 9 + 3 + 9 + 3 + 12 + 3 + 9 + 3 + 14) : (j += 1) std.debug.print("-", .{});
    std.debug.print("\n", .{});

    for (thresholds, 0..) |T, k| {
        var buf: [16]u8 = undefined;
        const label = fmtT(&buf, T);

        var log_sum: f64 = 0;
        var spawn_count: usize = 0;
        var best: f64 = std.math.inf(f64);
        var worst: f64 = -std.math.inf(f64);
        var worst_n: usize = 0;
        var regression_count: usize = 0;
        var first_win_n: ?usize = null;

        for (0..n_rows) |r| {
            const ratio = ratios[k][r];
            const threads = thread_counts[k][r];

            // Only consider rows where this threshold actually caused a spawn.
            // Serial-degenerate rows carry no information about the threshold.
            if (threads <= 1) continue;

            spawn_count += 1;
            log_sum += @log(ratio);
            if (ratio < best) best = ratio;
            if (ratio > worst) {
                worst = ratio;
                worst_n = row_ns[r];
            }
            if (ratio > REGRESSION) regression_count += 1;
            if (first_win_n == null and ratio < 1.0 - NOISE) first_win_n = row_ns[r];
        }

        if (spawn_count == 0) {
            std.debug.print("{s:>10} | {s:>71}\n", .{ label, "never spawned" });
            continue;
        }

        const geo = @exp(log_sum / @as(f64, @floatFromInt(spawn_count)));

        var nbuf: [24]u8 = undefined;
        const worst_n_str = fmtN(&nbuf, worst_n);
        var nbuf2: [24]u8 = undefined;
        const first_win_str = if (first_win_n) |n_| fmtN(&nbuf2, n_) else "never";

        std.debug.print("{s:>10} | {d:>8.3}x | {d:>8.2}x | {d:>8.2}x | {s:>12} | {d:>4}/{d:<4} | {s:>14}\n", .{
            label, geo, best, worst, worst_n_str, regression_count, spawn_count, first_win_str,
        });
    }

    std.debug.print("\nColumn guide:\n", .{});
    std.debug.print("  Geomean*      — geo mean of ratios over rows where T actually spawned.\n", .{});
    std.debug.print("                  Biased: each T is averaged over a different set of rows,\n", .{});
    std.debug.print("                  so low-T numbers include small-N disasters that high-T\n", .{});
    std.debug.print("                  values skip. Use as a rough sort, not a verdict.\n", .{});
    std.debug.print("  Best          — best ratio observed (lowest is good).\n", .{});
    std.debug.print("  Worst         — worst ratio observed. This is the number a default has\n", .{});
    std.debug.print("                  to live with. A default that shows 1.50x here will hurt\n", .{});
    std.debug.print("                  real users whose N falls in its bad zone.\n", .{});
    std.debug.print("  Worst @ N     — N value where the worst ratio occurred.\n", .{});
    std.debug.print("  Reg. N        — (# rows with ratio > regression) / (# rows with spawn).\n", .{});
    std.debug.print("                  Fraction of spawn cases that were actual slowdowns.\n", .{});
    std.debug.print("  First win @ N — smallest N where T beat serial by more than the noise\n", .{});
    std.debug.print("                  band. Lower = more aggressive threshold pays off sooner.\n\n", .{});

    std.debug.print("How to choose a default:\n", .{});
    std.debug.print("  1. Rule out any T with Worst > ~1.15x — it's an unforced error.\n", .{});
    std.debug.print("  2. Among survivors, prefer lowest 'First win @ N' (earlier speedups).\n", .{});
    std.debug.print("  3. Break ties with Geomean*. Small differences at the bandwidth floor\n", .{});
    std.debug.print("     don't matter.\n", .{});
}

pub fn blas_lv2_threshold_calibration(init: std.process.Init) !void {
    @setEvalBranchQuota(10000);
    const io = init.io;
    var gpa: std.heap.DebugAllocator(.{}) = .{ .backing_allocator = init.gpa };
    defer _ = gpa.deinit();
    const allocator = gpa.allocator();

    const max_n: usize = 16_384;

    const a = try allocator.alloc(zsl.cf64, max_n * max_n);
    const x_vec = try allocator.alloc(zsl.cf64, max_n);
    const y_vec = try allocator.alloc(zsl.cf64, max_n);
    defer allocator.free(a);
    defer allocator.free(x_vec);
    defer allocator.free(y_vec);

    // Initialize with small values to prevent f64 infinity/NaN degradation
    // over hundreds of operations if beta is used.
    for (a, 0..) |*val, i| val.* = .{ .re = @as(f64, @floatFromInt(i % 100)) / 1000.0, .im = @as(f64, @floatFromInt(i % 66)) / 660.0 };
    for (x_vec, 0..) |*val, i| val.* = .{ .re = @as(f64, @floatFromInt(i % 33)) / 330.0, .im = @as(f64, @floatFromInt(i % 17)) / 170.0 };
    for (y_vec, 0..) |*val, i| val.* = .{ .re = @as(f64, @floatFromInt(i % 8)) / 80.0, .im = @as(f64, @floatFromInt(i % 4)) / 40.0 };

    const thresholds = [_]usize{
        16_384,
        32_768,
        65_536,
        131_072,
        262_144,
        524_288,
        1_048_576,
        2_097_152,
        4_194_304,
        8_388_608,
        16_777_216,
    };

    const options_max_threads = 64;

    const hw_threads: usize = std.Thread.getCpuCount() catch 1;
    const effective_cap: usize = @min(hw_threads, options_max_threads);

    // Count how many N values we'll iterate over.
    var n_rows: usize = 0;
    {
        var dim: usize = 64;
        while (dim <= max_n) : (dim *= 2) n_rows += 1;
    }

    const ratios = try allocator.alloc([]f64, thresholds.len);
    defer allocator.free(ratios);
    const thread_counts = try allocator.alloc([]usize, thresholds.len);
    defer allocator.free(thread_counts);
    const row_ns = try allocator.alloc(usize, n_rows);
    defer allocator.free(row_ns);
    for (ratios, thread_counts) |*r, *t| {
        r.* = try allocator.alloc(f64, n_rows);
        t.* = try allocator.alloc(usize, n_rows);
    }
    defer for (ratios, thread_counts) |r, t| {
        allocator.free(r);
        allocator.free(t);
    };

    std.debug.print("\n=== TRMV (Level 2) Threshold Calibration ===\n", .{});
    std.debug.print("Hardware threads: {d}, options.max_threads: {d}, effective cap: {d}\n", .{ hw_threads, options_max_threads, effective_cap });
    std.debug.print("Automatic-mode rule: threads = max(1, min(effective_cap, (M*N) / T))\n\n", .{});
    std.debug.print("Cell format: 'R.RRx(Nt)' -> R = parallel_time/serial_time, N = threads spawned\n", .{});
    std.debug.print("  R < 1.00 -> parallel wins            R > 1.00 -> parallel loses\n", .{});
    std.debug.print("  (1t)     -> (M*N)/T < 2, stayed serial\n\n", .{});

    // Header
    std.debug.print("{s:>11} | {s:>12}", .{ "Dim (NxN)", "Serial (ns)" });
    for (thresholds) |th| {
        var buf: [16]u8 = undefined;
        std.debug.print(" | {s:>12}", .{fmtT(&buf, th)});
    }
    std.debug.print(" | {s:>22}\n", .{"Best for this Dim"});

    const sep_len: usize = 11 + 3 + 12 + thresholds.len * 15 + 3 + 22;
    var i: usize = 0;
    while (i < sep_len) : (i += 1) std.debug.print("-", .{});
    std.debug.print("\n", .{});

    var row: usize = 0;
    var current_n: usize = 64;
    while (current_n <= max_n) : (current_n *= 2) {
        row_ns[row] = current_n;
        const current_elements = current_n * current_n;

        // Scaling iterations down drastically as N grows quadratically
        const iters: usize =
            if (current_elements < 1_000_000) 250 else if (current_elements < 16_000_000) 100 else if (current_elements < 64_000_000) 25 else 10;

        // Serial baseline
        try zsl.linalg.blas.trmv(.col_major, .upper, .no_trans, .non_unit, current_n, a.ptr, max_n, x_vec.ptr, 1, .{ .num_threads = 1 });
        const s0 = std.Io.Clock.real.now(io);
        for (0..iters) |_| {
            try zsl.linalg.blas.trmv(.col_major, .upper, .no_trans, .non_unit, current_n, a.ptr, max_n, x_vec.ptr, 1, .{ .num_threads = 1 });
        }
        const s1 = std.Io.Clock.real.now(io);
        const serial_ns = (zsl.numeric.cast(f128, s1.toNanoseconds()) - zsl.numeric.cast(f128, s0.toNanoseconds())) / @as(f128, iters);

        std.debug.print("{d:>11} | {d:>12.0}", .{ current_n, @as(f64, @floatCast(serial_ns)) });

        var best_ratio: f128 = 1.0;
        var best_idx: ?usize = null;
        var best_threads: usize = 1;

        for (thresholds, 0..) |T, k| {
            try zsl.linalg.blas.trmv(.col_major, .upper, .no_trans, .non_unit, current_n, a.ptr, max_n, x_vec.ptr, 1, .{ .num_threads = 0, .parallel_threshold = T, .work = y_vec.ptr });
            const p0 = std.Io.Clock.real.now(io);
            for (0..iters) |_| {
                try zsl.linalg.blas.trmv(.col_major, .upper, .no_trans, .non_unit, current_n, a.ptr, max_n, x_vec.ptr, 1, .{ .num_threads = 0, .parallel_threshold = T, .work = y_vec.ptr });
            }
            const p1 = std.Io.Clock.real.now(io);
            const ns = (zsl.numeric.cast(f128, p1.toNanoseconds()) - zsl.numeric.cast(f128, p0.toNanoseconds())) / @as(f128, iters);
            const ratio: f128 = ns / serial_ns;

            const threads = @max(@as(usize, 1), @min(current_elements / T, effective_cap));

            ratios[k][row] = @as(f64, @floatCast(ratio));
            thread_counts[k][row] = threads;

            std.debug.print(" | {d:>6.2}x({d:>2}t)", .{
                @as(f64, @floatCast(ratio)), threads,
            });

            if (ratio < best_ratio) {
                best_ratio = ratio;
                best_idx = k;
                best_threads = threads;
            }
        }

        if (best_idx) |k| {
            var tbuf: [16]u8 = undefined;
            var sbuf: [48]u8 = undefined;
            const tl = fmtT(&tbuf, thresholds[k]);
            const sum = std.fmt.bufPrint(&sbuf, "{s} @ {d}t ({d:.2}x)", .{ tl, best_threads, @as(f64, @floatCast(best_ratio)) }) catch "?";
            std.debug.print(" | {s:>22}\n", .{sum});
        } else {
            std.debug.print(" | {s:>22}\n", .{"serial (no parallel won)"});
        }

        row += 1;
    }

    // Summary metrics
    const NOISE: f64 = 0.05;
    const REGRESSION: f64 = 1.10;

    std.debug.print("\n=== Summary (noise band = ±{d:.0}%, regression = >{d:.0}%) ===\n\n", .{
        NOISE * 100, (REGRESSION - 1.0) * 100,
    });

    std.debug.print("{s:>10} | {s:>9} | {s:>9} | {s:>9} | {s:>12} | {s:>9} | {s:>14}\n", .{
        "Threshold", "Geomean*", "Best", "Worst", "Worst @ Dim", "Reg. Dim", "First win @ N",
    });
    var j: usize = 0;
    while (j < 10 + 3 + 9 + 3 + 9 + 3 + 9 + 3 + 12 + 3 + 9 + 3 + 14) : (j += 1) std.debug.print("-", .{});
    std.debug.print("\n", .{});

    for (thresholds, 0..) |T, k| {
        var buf: [16]u8 = undefined;
        const label = fmtT(&buf, T);

        var log_sum: f64 = 0;
        var spawn_count: usize = 0;
        var best: f64 = std.math.inf(f64);
        var worst: f64 = -std.math.inf(f64);
        var worst_n: usize = 0;
        var regression_count: usize = 0;
        var first_win_n: ?usize = null;

        for (0..n_rows) |r| {
            const ratio = ratios[k][r];
            const threads = thread_counts[k][r];

            if (threads <= 1) continue;

            spawn_count += 1;
            log_sum += @log(ratio);
            if (ratio < best) best = ratio;
            if (ratio > worst) {
                worst = ratio;
                worst_n = row_ns[r];
            }
            if (ratio > REGRESSION) regression_count += 1;
            if (first_win_n == null and ratio < 1.0 - NOISE) first_win_n = row_ns[r];
        }

        if (spawn_count == 0) {
            std.debug.print("{s:>10} | {s:>71}\n", .{ label, "never spawned" });
            continue;
        }

        const geo = @exp(log_sum / @as(f64, @floatFromInt(spawn_count)));

        var nbuf: [24]u8 = undefined;
        const worst_n_str = fmtN(&nbuf, worst_n);
        var nbuf2: [24]u8 = undefined;
        const first_win_str = if (first_win_n) |n_| fmtN(&nbuf2, n_) else "never";

        std.debug.print("{s:>10} | {d:>8.3}x | {d:>8.2}x | {d:>8.2}x | {s:>12} | {d:>4}/{d:<4} | {s:>14}\n", .{
            label, geo, best, worst, worst_n_str, regression_count, spawn_count, first_win_str,
        });
    }
}

fn fmtT(buf: []u8, x: usize) []const u8 {
    if (x >= 1024 * 1024)
        return std.fmt.bufPrint(buf, "T={d}Mi", .{x / (1024 * 1024)}) catch "?";
    if (x >= 1024)
        return std.fmt.bufPrint(buf, "T={d}Ki", .{x / 1024}) catch "?";
    return std.fmt.bufPrint(buf, "T={d}", .{x}) catch "?";
}

fn fmtN(buf: []u8, x: usize) []const u8 {
    const ux: usize = @intCast(x);
    if (ux >= 1024 * 1024)
        return std.fmt.bufPrint(buf, "{d}Mi", .{ux / (1024 * 1024)}) catch "?";
    if (ux >= 1024)
        return std.fmt.bufPrint(buf, "{d}Ki", .{ux / 1024}) catch "?";
    return std.fmt.bufPrint(buf, "{d}", .{ux}) catch "?";
}

/// Formats the dyadic exactly into a base-10 string.
/// The caller takes ownership of the returned slice and must free it.
pub fn dyadicToString(allocator: std.mem.Allocator, x: anytype) ![]u8 {
    // Handle edge cases
    if (x.isNan())
        return allocator.dupe(u8, "NaN");

    if (x.isInf())
        return allocator.dupe(u8, if (x.positive) "Inf" else "-Inf");

    if (x.isZero())
        return allocator.dupe(u8, if (x.positive) "0.0" else "-0.0");

    const Mantissa = @TypeOf(x).Mantissa;
    const bits = @typeInfo(Mantissa).int.bits;

    // Calculate max significant figures at compile time
    const calculated_sig_figs = comptime @as(usize, @intFromFloat(@as(f64, @floatFromInt(bits)) * 0.3010299956639812 * 1.05));
    const max_sig_figs: usize = if (calculated_sig_figs == 0) 1 else calculated_sig_figs;

    const Helpers = struct {
        // Add two base-10 digit arrays
        fn addDigits(alloc: std.mem.Allocator, dest: *std.ArrayList(u8), src: []const u8) !void {
            var carry: u8 = 0;
            var i: usize = 0;
            while (i < src.len or carry > 0) : (i += 1) {
                if (i >= dest.items.len) {
                    try dest.append(alloc, 0);
                }
                const a = dest.items[i];
                const b = if (i < src.len) src[i] else 0;
                const sum = a + b + carry;
                dest.items[i] = sum % 10;
                carry = sum / 10;
            }
        }

        // Multiply a base-10 digit array by a small scalar
        fn mulByDigit(alloc: std.mem.Allocator, dest: *std.ArrayList(u8), multiplier: u32) !void {
            var carry: u32 = 0;
            for (dest.items) |*digit| {
                const prod = @as(u32, digit.*) * multiplier + carry;
                digit.* = @intCast(prod % 10);
                carry = prod / 10;
            }
            while (carry > 0) {
                try dest.append(alloc, @intCast(carry % 10));
                carry /= 10;
            }
        }
    };

    var num: std.ArrayList(u8) = .empty;
    defer num.deinit(allocator);
    try num.append(allocator, 0);

    var power_of_2: std.ArrayList(u8) = .empty;
    defer power_of_2.deinit(allocator);
    try power_of_2.append(allocator, 1);

    // Step 1: Load the base mantissa into the `num` bignum array
    var i: u16 = 0;
    while (i < bits) : (i += 1) {
        const shift: std.math.Log2Int(Mantissa) = @intCast(i);
        const bit = (x.mantissa >> shift) & 1;
        if (bit == 1) {
            try Helpers.addDigits(allocator, &num, power_of_2.items);
        }
        try Helpers.mulByDigit(allocator, &power_of_2, 2);
    }

    // Step 2: Apply the exponent
    const e_i32: i32 = x.exponent;
    const abs_e: u32 = @intCast(if (e_i32 < 0) -e_i32 else e_i32);

    if (e_i32 > 0) {
        // Number is Mantissa * 2^E
        var j: u32 = 0;
        while (j < abs_e) : (j += 1) {
            try Helpers.mulByDigit(allocator, &num, 2);
        }
    } else if (e_i32 < 0) {
        // Number is Mantissa * 5^|E| / 10^|E|
        // We compute the numerator exactly to avoid floating-point logic
        var j: u32 = 0;
        while (j < abs_e) : (j += 1) {
            try Helpers.mulByDigit(allocator, &num, 5);
        }
    }

    // Pad with trailing zeros to guarantee we safely cross the decimal point
    if (e_i32 < 0) {
        while (num.items.len <= abs_e) {
            try num.append(allocator, 0);
        }
    }

    // Step 3: Format the final string
    var result: std.ArrayList(u8) = .empty;
    if (!x.positive) {
        try result.append(allocator, '-');
    }

    var sig_figs: usize = 0;
    var started_counting = false;

    if (e_i32 >= 0) {
        // Integer representation
        var idx: usize = num.items.len;
        while (idx > 0) {
            idx -= 1;
            const digit = num.items[idx];

            if (digit != 0) started_counting = true;
            if (started_counting) sig_figs += 1;

            if (sig_figs <= max_sig_figs) {
                try result.append(allocator, digit + '0');
            } else {
                // Cap exceeded: Pad integer part with zeros for magnitude safety
                try result.append(allocator, '0');
            }
        }
        try result.appendSlice(allocator, ".0");
    } else {
        // Fractional representation, inject decimal point
        var idx: usize = num.items.len;
        while (idx > 0) {
            idx -= 1;
            const digit = num.items[idx];

            if (digit != 0 and !started_counting) started_counting = true;
            if (started_counting) sig_figs += 1;

            if (sig_figs <= max_sig_figs or !started_counting) {
                try result.append(allocator, digit + '0');
            } else {
                // Cap exceeded
                if (idx >= abs_e) {
                    // Haven't hit the decimal yet, pad to preserve scale
                    try result.append(allocator, '0');
                } else {
                    // Safely right of the decimal point, just truncate
                    break;
                }
            }

            if (idx == abs_e) {
                try result.append(allocator, '.');
            }
        }

        // Cleanup: If the cap caused the string to stop exactly on the decimal, append a '0'
        if (result.items.len > 0 and result.items[result.items.len - 1] == '.') {
            try result.append(allocator, '0');
        }
    }

    return result.toOwnedSlice(allocator);
}

// fn avg(values: []const f64) f64 {
//     var sum: f64 = 0;
//     for (values) |value| {
//         sum += value;
//     }
//     return zsl.float.div(sum, values.len);
// }

// fn avg_complex(values: []const zsl.cf64) f64 {
//     var sum: f64 = 0;
//     for (values) |value| {
//         sum += zsl.cfloat.abs(value);
//     }
//     return zsl.float.div(sum, values.len);
// }

// fn random_buffer(
//     allocator: std.mem.Allocator,
//     size: u32,
// ) ![]f64 {
//     var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
//     const rand = prng.random();

//     const buffer = try allocator.alloc(f64, size);
//     for (0..buffer.len) |i| {
//         buffer[i] = rand.float(f64);
//     }
//     return buffer;
// }

// fn random_buffer_fill(
//     buffer: []f64,
// ) void {
//     var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
//     //var prng = std.Random.DefaultPrng.init(2); // fixed seed for reproducibility
//     const rand = prng.random();

//     for (0..buffer.len) |i| {
//         buffer[i] = rand.float(f64);
//     }
// }

// fn random_buffer_fill_complex(
//     buffer: []zsl.cf64,
// ) void {
//     //var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
//     var prng = std.Random.DefaultPrng.init(2); // fixed seed for reproducibility
//     const rand = prng.random();

//     for (0..buffer.len) |i| {
//         buffer[i] = zsl.cf64.init(rand.float(f64), rand.float(f64));
//     }
// }

// /// Generate a random m×n matrix A with specified 2-norm condition number `kappa`.
// pub fn random_matrix(
//     allocator: std.mem.Allocator,
//     m: u32,
//     n: u32,
//     kappa: f64,
//     order: zsl.Layout,
// ) ![]f64 {
//     // - `allocator`: memory allocator
//     // - `m`, `n`: dimensions
//     // - `kappa`: desired condition number κ₂(A)
//     // - `order`: either .RowMajor or .ColMajor
//     //
//     // Method:
//     // 1. Construct singular values σᵢ geometrically spaced between 1 and κ^(1/min(m,n)).
//     // 2. Generate random m×min(m,n) matrix X, QR factor → Q₁ (m×r).
//     // 3. Generate random n×min(m,n) matrix Y, QR factor → Q₂ (n×r).
//     // 4. Form A = Q₁ * diag(σ) * Q₂ᵀ.
//     //     - Compute B = diag(σ) * Q₂ᵀ.
//     //     - Then A = Q₁ × B.
//     const r = zsl.int.min(m, n);

//     // 1) geometric singular values σ[i]
//     var sig = try allocator.alloc(f64, r);
//     defer allocator.free(sig);
//     for (0..r) |i| {
//         sig[i] = zsl.float.pow(kappa, zsl.float.div(zsl.scast(f64, i), r - 1));
//     }

//     // 2) random X: m×r → QR → Q₁
//     const X = try random_buffer(allocator, m * r);
//     defer allocator.free(X);
//     const tauX = try allocator.alloc(f64, r);
//     defer allocator.free(tauX);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     // X now holds Q₁

//     // 3) random Y: n×r → QR → Q₂
//     const Y = try random_buffer(allocator, n * r);
//     defer allocator.free(Y);
//     const tauY = try allocator.alloc(f64, r);
//     defer allocator.free(tauY);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     // Y now holds Q₂

//     // 4a) form B = diag(sig) * Q2ᵀ (B is r×n in same layout)
//     const B = try allocator.alloc(f64, r * n);
//     defer allocator.free(B);
//     for (0..r) |i| {
//         const row_start = if (order == .col_major)
//             B.ptr + i
//         else
//             B.ptr + i * n;

//         const stride = if (order == .col_major) r else 1;

//         for (0..n) |j| {
//             // Q2[j,i] is at Y.ptr + j*ldY + i
//             const q2ji = (Y.ptr + j * if (order == .col_major) n else r)[i];
//             row_start[j * stride] = sig[i] * q2ji;
//         }
//     }

//     // 4b) A = Q1 (m×r) × B (r×n) → A (m×n)
//     const A = try allocator.alloc(f64, m * n);
//     zsl.linalg.blas.dgemm(
//         order,
//         .no_trans,
//         .no_trans,
//         zsl.scast(i32, m),
//         zsl.scast(i32, n),
//         zsl.scast(i32, r),
//         1,
//         X.ptr,
//         zsl.scast(i32, if (order == .col_major) m else r),
//         B.ptr,
//         zsl.scast(i32, if (order == .col_major) r else n),
//         0,
//         A.ptr,
//         zsl.scast(i32, if (order == .col_major) m else n),
//     );

//     return A;
// }

// pub fn random_matrix_buffer(
//     allocator: std.mem.Allocator,
//     m: u32,
//     n: u32,
//     kappa: f64,
//     A: []f64,
//     order: zsl.Order,
// ) !void {
//     // - `allocator`: memory allocator
//     // - `m`, `n`: dimensions
//     // - `kappa`: desired condition number κ₂(A)
//     // - `order`: either .RowMajor or .ColMajor
//     //
//     // Method:
//     // 1. Construct singular values σᵢ geometrically spaced between 1 and κ^(1/min(m,n)).
//     // 2. Generate random m×min(m,n) matrix X, QR factor → Q₁ (m×r).
//     // 3. Generate random n×min(m,n) matrix Y, QR factor → Q₂ (n×r).
//     // 4. Form A = Q₁ * diag(σ) * Q₂ᵀ.
//     //     - Compute B = diag(σ) * Q₂ᵀ.
//     //     - Then A = Q₁ × B.
//     const r = zsl.int.min(m, n);

//     // 1) geometric singular values σ[i]
//     var sig = try allocator.alloc(f64, r);
//     defer allocator.free(sig);
//     for (0..r) |i| {
//         sig[i] = zsl.float.pow(kappa, zsl.float.div(zsl.scast(f64, i), r - 1));
//     }

//     // 2) random X: m×r → QR → Q₁
//     const X = try random_buffer(allocator, m * r);
//     defer allocator.free(X);
//     const tauX = try allocator.alloc(f64, r);
//     defer allocator.free(tauX);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, m),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         X.ptr,
//         zsl.scast(c_int, if (order == .col_major) m else r),
//         tauX.ptr,
//     );
//     // X now holds Q₁

//     // 3) random Y: n×r → QR → Q₂
//     const Y = try random_buffer(allocator, n * r);
//     defer allocator.free(Y);
//     const tauY = try allocator.alloc(f64, r);
//     defer allocator.free(tauY);
//     _ = ci.LAPACKE_dgeqrf(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         order.toCInt(),
//         zsl.scast(c_int, n),
//         zsl.scast(c_int, r),
//         zsl.scast(c_int, r),
//         Y.ptr,
//         zsl.scast(c_int, if (order == .col_major) n else r),
//         tauY.ptr,
//     );
//     // Y now holds Q₂

//     // 4a) form B = diag(sig) * Q2ᵀ (B is r×n in same layout)
//     const B = try allocator.alloc(f64, r * n);
//     defer allocator.free(B);
//     for (0..r) |i| {
//         const row_start = if (order == .col_major)
//             B.ptr + i
//         else
//             B.ptr + i * n;

//         const stride = if (order == .col_major) r else 1;

//         for (0..n) |j| {
//             // Q2[j,i] is at Y.ptr + j*ldY + i
//             const q2ji = (Y.ptr + j * if (order == .col_major) n else r)[i];
//             row_start[j * stride] = sig[i] * q2ji;
//         }
//     }

//     // 4b) A = Q1 (m×r) × B (r×n) → A (m×n)
//     zsl.linalg.blas.dgemm(
//         order,
//         .no_trans,
//         .no_trans,
//         zsl.scast(i32, m),
//         zsl.scast(i32, n),
//         zsl.scast(i32, r),
//         1,
//         X.ptr,
//         zsl.scast(i32, if (order == .col_major) m else r),
//         B.ptr,
//         zsl.scast(i32, if (order == .col_major) r else n),
//         0,
//         A.ptr,
//         zsl.scast(i32, if (order == .col_major) m else n),
//     );
// }

// fn max_difference(a: []const f64, b: []const f64) struct {
//     index: u32,
//     value: f64,
// } {
//     std.debug.assert(a.len == b.len);

//     var max_diff: f64 = 0;
//     var max_index: u32 = 0;

//     var i: u32 = 0;
//     while (i < a.len) : (i += 1) {
//         const diff = zsl.float.abs(a[i] - b[i]);
//         if (diff > max_diff) {
//             max_diff = diff;
//             max_index = i;
//         }
//     }

//     return .{ .index = max_index, .value = max_diff };
// }

// fn random_symmetric_matrix(
//     allocator: std.mem.Allocator,
//     size: u32,
//     factor: f64,
// ) ![]f64 {
//     var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
//     const rand = prng.random();

//     const matrix = try allocator.alloc(f64, size * size);

//     for (0..size) |i| {
//         for (i + 1..size) |j| {
//             const value = rand.float(f64) * factor;
//             matrix[i * size + j] = value;
//             matrix[j * size + i] = value;
//         }
//     }

//     return matrix;
// }

// /// Generates a random symmetric positive definite matrix with a specified condition number.
// fn random_symmetric_positive_definite_matrix(
//     allocator: std.mem.Allocator,
//     size: u32,
//     cond_target: f64,
// ) ![]f64 {
//     // 1) compute λ_i geometrically between 1 and cond_target
//     var lambdas = try allocator.alloc(f64, size);
//     defer allocator.free(lambdas);
//     const l_min = 1;
//     const l_max = cond_target;
//     for (0..size) |i| {
//         lambdas[i] = l_min * zsl.float.pow(l_max / l_min, zsl.float.div(zsl.scast(f64, i), size - 1));
//     }

//     // 2) random X, then QR -> get Q (n×n)
//     const X = try random_buffer(allocator, size * size);
//     defer allocator.free(X);
//     const tau = try allocator.alloc(f64, size);
//     defer allocator.free(tau);
//     // QR factorization in-place: X = Q*R
//     _ = ci.LAPACKE_dgeqrf(
//         zsl.Order.col_major.toCInt(),
//         zsl.scast(c_int, size),
//         zsl.scast(c_int, size),
//         X.ptr,
//         zsl.scast(c_int, size),
//         tau.ptr,
//     );
//     _ = ci.LAPACKE_dorgqr(
//         zsl.Order.col_major.toCInt(),
//         zsl.scast(c_int, size),
//         zsl.scast(c_int, size),
//         zsl.scast(c_int, size),
//         X.ptr,
//         zsl.scast(c_int, size),
//         tau.ptr,
//     );
//     // now X holds Q

//     // 3) form B = D * Q^T, where D = diag(λ_1, λ_2, ..., λ_n)
//     var B = try allocator.alloc(f64, size * size);
//     defer allocator.free(B);
//     // copy Q^T into B
//     for (0..size) |i| {
//         for (0..size) |j|
//             B[i * size + j] = X[j * size + i];
//     }
//     // scale row i of B by λ_i
//     for (0..size) |i| {
//         zsl.linalg.blas.dscal(zsl.scast(i32, size), lambdas[i], B.ptr + i * size, 1);
//     }

//     // 4) Assemble A = Q * D * Q^T as A = Q * B
//     const A = try allocator.alloc(f64, size * size);
//     zsl.linalg.blas.dgemm(
//         .col_major,
//         .no_trans,
//         .no_trans,
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         1,
//         X.ptr,
//         zsl.scast(i32, size),
//         B.ptr,
//         zsl.scast(i32, size),
//         0,
//         A.ptr,
//         zsl.scast(i32, size),
//     );

//     return A;
// }

// fn frobernius_norm_difference(a: []const f64, b: []const f64) f64 {
//     std.debug.assert(a.len == b.len);

//     var norm: f64 = 0;
//     for (0..a.len) |i| {
//         const diff = a[i] - b[i];
//         norm += diff * diff;
//     }
//     return zsl.float.sqrt(norm);
// }

// fn frobernius_norm_difference_matrix(a: anytype, b: anytype) !f64 {
//     const m = if (comptime zsl.meta.isSymmetricMatrix(@TypeOf(a)) or
//         zsl.meta.isHermitianMatrix(@TypeOf(a)) or
//         zsl.meta.isTridiagonalMatrix(@TypeOf(a)) or
//         zsl.meta.isPermutationMatrix(@TypeOf(a)))
//         a.size
//     else
//         a.rows;
//     const n = if (comptime zsl.meta.isSymmetricMatrix(@TypeOf(a)) or
//         zsl.meta.isHermitianMatrix(@TypeOf(a)) or
//         zsl.meta.isTridiagonalMatrix(@TypeOf(a)) or
//         zsl.meta.isPermutationMatrix(@TypeOf(a)))
//         a.size
//     else
//         a.cols;

//     var norm: f64 = 0;

//     var i: u32 = 0;
//     while (i < m) : (i += 1) {
//         var j: u32 = 0;
//         while (j < n) : (j += 1) {
//             if (comptime !zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(a))) and !zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(b)))) {
//                 const diff = try a.get(i, j) - try b.get(i, j);
//                 norm += diff * diff;
//             } else {
//                 const diff = try zsl.sub(try a.get(i, j), try b.get(i, j), .{});
//                 norm += try zsl.abs2(diff, .{});
//             }
//         }
//     }

//     return zsl.float.sqrt(norm);
// }

// fn is_symmetric(a: []const f64, size: u32) bool {
//     for (0..size) |i| {
//         for (i + 1..size) |j| {
//             if (!std.math.approxEqRel(f64, a[i * size + j], a[j * size + i], 1e-9)) {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

// fn random_complex_matrix(
//     allocator: std.mem.Allocator,
//     rows: u32,
//     cols: u32,
//     factor: f64,
// ) ![]zsl.cf64 {
//     var prng = std.Random.DefaultPrng.init(@bitCast(std.time.timestamp()));
//     const rand = prng.random();

//     const matrix = try allocator.alloc(zsl.cf64, rows * cols);
//     for (0..matrix.len) |i| {
//         matrix[i] = zsl.cf64.init(rand.float(f64) * factor, rand.float(f64) * factor);
//     }
//     return matrix;
// }

// fn random_complex_hermitian_positive_definite_matrix(
//     allocator: std.mem.Allocator,
//     size: u32,
//     factor: f64,
// ) ![]zsl.cf64 {
//     // allocate M
//     const M = try random_complex_matrix(allocator, size, size, factor);
//     // allocate A = M M^H
//     const A = try allocator.alloc(zsl.cf64, size * size);

//     // A = M × M^H
//     zsl.linalg.blas.zgemm(
//         .row_major,
//         .no_trans,
//         .conj_trans,
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         zsl.scast(i32, size),
//         zsl.cf64.init(1, 0),
//         M.ptr,
//         zsl.scast(i32, size),
//         M.ptr,
//         zsl.scast(i32, size),
//         zsl.cf64.init(0, 0),
//         A.ptr,
//         zsl.scast(i32, size),
//     );

//     allocator.free(M);
//     return A;
// }

// fn frobenius_norm_complex_difference(a: []const zsl.cf64, b: []const zsl.cf64) f64 {
//     std.debug.assert(a.len == b.len);

//     var norm: f64 = 0;
//     for (0..a.len) |i| {
//         const diff = zsl.cfloat.sub(a[i], b[i]);
//         norm += diff.re * diff.re + diff.im * diff.im;
//     }
//     return zsl.float.sqrt(norm);
// }

// fn is_hermitian(a: []const zsl.cf64, size: u32) bool {
//     for (0..size) |i| {
//         for (i + 1..size) |j| {
//             if (!std.math.approxEqRel(f64, a[i * size + j].re, a[j * size + i].re, 1e-9) or
//                 !std.math.approxEqRel(f64, a[i * size + j].im, -a[j * size + i].im, 1e-9))
//             {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

// fn print_complex_matrix(desc: []const u8, m: u32, n: u32, a: []zsl.cf64, lda: u32, order: zsl.Order) void {
//     std.debug.print("\n{s}\n", .{desc});
//     if (order == .row_major) {
//         var i: u32 = 0;
//         while (i < m) : (i += 1) {
//             var j: u32 = 0;
//             while (j < n) : (j += 1) {
//                 std.debug.print("{d:.4} + {d:.4}i  ", .{ a[i * lda + j].re, a[i * lda + j].im });
//             }
//             std.debug.print("\n", .{});
//         }
//     } else {
//         var i: u32 = 0;
//         while (i < m) : (i += 1) {
//             var j: u32 = 0;
//             while (j < n) : (j += 1) {
//                 std.debug.print("{d:.4} + {d:.4}i  ", .{ a[i + j * lda].re, a[i + j * lda].im });
//             }
//             std.debug.print("\n", .{});
//         }
//     }
// }

fn randomPermutation(rand: std.Random, data: []usize) void {
    // Initialize with identity permutation
    var i: usize = 0;
    while (i < data.len) : (i += 1) {
        data[i] = i;
    }

    // Shuffle using Fisher-Yates algorithm
    i = data.len - 1;
    while (i > 0) : (i -= 1) {
        const j = rand.intRangeAtMost(usize, 0, i);
        const temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }
}

fn randomMatrix(comptime M: type, allocator: std.mem.Allocator, rand: std.Random, rows: usize, cols: usize) !M {
    switch (comptime zsl.meta.matrixType(M)) {
        .general_dense => {
            var result: M = try .init(allocator, rows, cols);

            var i: usize = 0;
            while (i < rows) : (i += 1) {
                var j: usize = 0;
                while (j < cols) : (j += 1) {
                    result.set(
                        i,
                        j,
                        zsl.numeric.cast(
                            zsl.meta.Numeric(M),
                            if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                            else
                                rand.float(f64),
                        ),
                    ) catch unreachable;
                }
            }

            return result;
        },
        .symmetric_dense, .hermitian_dense => {
            var result: M = try .init(allocator, rows);

            var i: usize = 0;
            while (i < rows) : (i += 1) {
                var j: usize = i;
                while (j < rows) : (j += 1) {
                    result.set(
                        i,
                        j,
                        zsl.numeric.cast(
                            zsl.meta.Numeric(M),
                            if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                zsl.cf64{ .re = rand.float(f64), .im = if ((comptime zsl.meta.isHermitianMatrix(M)) and i == j) 0.0 else rand.float(f64) }
                            else
                                rand.float(f64),
                        ),
                    ) catch unreachable;
                }
            }

            return result;
        },
        .triangular_dense => {
            var result: M = try M.init(allocator, rows, cols);

            if (comptime zsl.meta.uploOf(M) == .upper) {
                var i: usize = 0;
                while (i < zsl.int.min(rows, cols)) : (i += 1) {
                    if (comptime zsl.meta.diagOf(M) == .non_unit) {
                        result.set(
                            i,
                            i,
                            zsl.numeric.cast(
                                zsl.meta.Numeric(M),
                                if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                    zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                                else
                                    rand.float(f64),
                            ),
                        ) catch unreachable;
                    }

                    var j: usize = i + 1;
                    while (j < cols) : (j += 1) {
                        result.set(
                            i,
                            j,
                            zsl.numeric.cast(
                                zsl.meta.Numeric(M),
                                if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                    zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                                else
                                    rand.float(f64),
                            ),
                        ) catch unreachable;
                    }
                }
            } else {
                var i: usize = 0;
                while (i < rows) : (i += 1) {
                    var j: usize = 0;
                    while (j < i and j < cols) : (j += 1) {
                        result.set(
                            i,
                            j,
                            zsl.numeric.cast(
                                zsl.meta.Numeric(M),
                                if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                    zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                                else
                                    rand.float(f64),
                            ),
                        ) catch unreachable;
                    }

                    if ((comptime zsl.meta.diagOf(M) == .non_unit) and i < cols) {
                        result.set(
                            i,
                            i,
                            zsl.numeric.cast(
                                zsl.meta.Numeric(M),
                                if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                                    zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                                else
                                    rand.float(f64),
                            ),
                        ) catch unreachable;
                    }
                }
            }

            return result;
        },
        .general_sparse => {
            const nnz: usize = (rows * cols) / 100;

            var builder: zsl.matrix.builder.Sparse(zsl.meta.Numeric(M)) = try .init(allocator, rows, cols, nnz);
            errdefer builder.deinit(allocator);

            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                builder.appendAssumeCapacity(
                    rand.intRangeAtMost(usize, 0, rows - 1),
                    rand.intRangeAtMost(usize, 0, cols - 1),
                    zsl.numeric.cast(
                        zsl.meta.Numeric(M),
                        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                            zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                        else
                            rand.float(f64),
                    ),
                );
            }

            return try builder.compile(allocator, zsl.meta.layoutOf(M));
        },
        .symmetric_sparse, .hermitian_sparse => {
            const nnz: usize = (rows * cols) / 100;

            var builder: zsl.matrix.builder.Sparse(zsl.meta.Numeric(M)) = try .init(allocator, rows, rows, nnz);
            errdefer builder.deinit(allocator);

            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                const r = rand.intRangeAtMost(usize, 0, rows - 1);
                const c = rand.intRangeAtMost(usize, 0, cols - 1);

                builder.appendAssumeCapacity(
                    r,
                    c,
                    zsl.numeric.cast(
                        zsl.meta.Numeric(M),
                        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                            zsl.cf64{ .re = rand.float(f64), .im = if ((comptime zsl.meta.isHermitianMatrix(M)) and r == c) 0.0 else rand.float(f64) }
                        else
                            rand.float(f64),
                    ),
                );
            }

            return if (comptime zsl.meta.isSymmetricSparseMatrix(M))
                builder.compileSymmetric(allocator, zsl.meta.uploOf(M), zsl.meta.layoutOf(M))
            else
                builder.compileHermitian(allocator, zsl.meta.uploOf(M), zsl.meta.layoutOf(M));
        },
        .triangular_sparse => {
            const nnz: usize = (rows * cols) / 100;

            var builder: zsl.matrix.builder.Sparse(zsl.meta.Numeric(M)) = try .init(allocator, rows, cols, nnz);
            errdefer builder.deinit(allocator);

            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                const r = rand.intRangeAtMost(usize, 0, rows - 1);
                const c = rand.intRangeAtMost(usize, 0, cols - 1);

                builder.appendAssumeCapacity(
                    r,
                    c,
                    zsl.numeric.cast(
                        zsl.meta.Numeric(M),
                        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                            zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                        else
                            rand.float(f64),
                    ),
                );
            }

            return builder.compileTriangular(allocator, zsl.meta.uploOf(M), zsl.meta.diagOf(M), zsl.meta.layoutOf(M));
        },
        .diagonal => {
            var result: M = try .init(allocator, rows, cols);
            errdefer result.deinit(allocator);

            var i: usize = 0;
            while (i < zsl.int.min(rows, cols)) : (i += 1) {
                result.set(
                    i,
                    i,
                    zsl.numeric.cast(
                        zsl.meta.Numeric(M),
                        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(M)))
                            zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                        else
                            rand.float(f64),
                    ),
                ) catch unreachable;
            }

            return result;
        },
        .permutation => {
            var result: M = try .init(allocator, rows);
            errdefer result.deinit(allocator);

            randomPermutation(rand, result.data[0..rows]);

            return result;
        },
        else => unreachable,
    }
}

fn printMatrix(desc: []const u8, A: anytype) void {
    std.debug.print("\nMatrix {s}:\n\n", .{desc});

    var i: u32 = 0;
    while (i < A.rows) : (i += 1) {
        std.debug.print("\t", .{});

        var j: u32 = 0;
        while (j < A.cols) : (j += 1) {
            // if (comptime zsl.meta.isComplex(zsl.meta.Numeric(@TypeOf(A)))) {
            //     std.debug.print("{d:7.4} + {d:7.4}i    ", .{ (A.get(i, j) catch unreachable).re, (A.get(i, j) catch unreachable).im });
            // } else {
            //     std.debug.print("{d:5.4}    ", .{A.get(i, j) catch unreachable});
            // }
            std.debug.print("{}    ", .{A.get(i, j) catch unreachable});
        }
        std.debug.print("\n", .{});
    }
    std.debug.print("\n", .{});
}

fn randomVector(comptime V: type, allocator: std.mem.Allocator, rand: std.Random, len: usize) !V {
    switch (comptime zsl.meta.vectorType(V)) {
        .dense => {
            var result: V = try .init(allocator, len);

            var i: usize = 0;
            while (i < len) : (i += 1) {
                result.set(
                    i,
                    zsl.numeric.cast(
                        zsl.meta.Numeric(V),
                        if (comptime zsl.meta.isComplex(zsl.meta.Numeric(V)))
                            zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                        else
                            rand.float(f64),
                    ),
                ) catch unreachable;
            }

            return result;
        },
        .sparse => {
            const nnz: usize = zsl.int.max(1, rand.intRangeAtMost(usize, len / 10, len / 2));

            var result: V = try .init(allocator, len, nnz);
            errdefer result.deinit(allocator);

            // generate random indices
            var used: std.AutoHashMap(usize, void) = .init(allocator);
            defer used.deinit();
            var count: usize = 0;
            while (count < nnz) : (count += 1) {
                const i = rand.intRangeAtMost(usize, 0, len - 1);
                if (!used.contains(i)) {
                    try used.put(i, {});
                    try result.set(
                        allocator,
                        i,
                        zsl.numeric.cast(
                            zsl.meta.Numeric(V),
                            if (comptime zsl.meta.isComplex(zsl.meta.Numeric(V)))
                                zsl.cf64{ .re = rand.float(f64), .im = rand.float(f64) }
                            else
                                rand.float(f64),
                        ),
                    );
                } else {
                    count -= 1; // try again
                }
            }

            return result;
        },
        else => unreachable,
    }
}

fn printVector(desc: []const u8, v: anytype) void {
    std.debug.print("\nVector {s}:\n", .{desc});

    var i: usize = 0;
    while (i < v.len) : (i += 1) {
        std.debug.print("{}\n", .{v.get(i) catch unreachable});
    }
    std.debug.print("\n", .{});
}
