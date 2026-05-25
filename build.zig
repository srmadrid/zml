const std = @import("std");
const builtin = @import("builtin");

const IntMode = enum {
    default,
    wrap,
    saturate,
};

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    const options = b.addOptions();

    const opt_int_mode = b.option(IntMode, "int_mode", "Integer operation mode") orelse IntMode.wrap;
    options.addOption(IntMode, "int_mode", opt_int_mode);

    const opt_max_threads = b.option(usize, "max_threads", "Maximum number of threads") orelse 64;
    options.addOption(usize, "max_threads", opt_max_threads);

    const opt_max_dimensions = b.option(usize, "max_dimensions", "Maximum number of dimensions for dense and strided arrays") orelse 8;
    options.addOption(usize, "max_dimensions", opt_max_dimensions);

    const opt_link_cblas = b.option([]const u8, "link_cblas", "Link CBLAS implementation; required if calling any function from linalg.cblas");
    options.addOption(?[]const u8, "link_cblas", opt_link_cblas);

    const opt_link_lapacke = b.option([]const u8, "link_lapacke", "Link LAPACKE implementation; required if calling any function from linalg.lapacke");
    options.addOption(?[]const u8, "link_lapacke", opt_link_lapacke);

    inline for (1..4) |cache_level| {
        const opt_l_size = b.option(usize, std.fmt.comptimePrint("l{d}_size", .{cache_level}), "Override total L1 data cache size in bytes") orelse
            if (target.result.cpu.arch == builtin.cpu.arch and
                target.result.os.tag == builtin.os.tag)
                getCacheSize(b, cache_level)
            else blk: {
                const default = getDefaultCacheSize(cache_level);

                std.debug.print(
                    "warning: cross-compiling ({s}-{s} → {s}-{s}), cannot auto-detect L{d} size; " ++
                        "pass -Dl{d}-size=<bytes> to tune, defaulting to {d}KB\n",
                    .{
                        @tagName(builtin.cpu.arch),       @tagName(builtin.os.tag),
                        @tagName(target.result.cpu.arch), @tagName(target.result.os.tag),
                        cache_level,                      cache_level,
                        default,
                    },
                );

                break :blk default * 1024;
            };
        options.addOption(usize, std.fmt.comptimePrint("l{d}_size", .{cache_level}), opt_l_size);
    }

    const module = b.addModule("zsl", .{
        .root_source_file = b.path("src/zsl.zig"),
        .target = target,
        .optimize = optimize,
    });
    module.addOptions("options", options);

    // Executable (for testing)
    const exe = b.addExecutable(.{
        .name = "main",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    exe.root_module.addImport("zsl", module);

    if (opt_link_cblas != null) {
        exe.root_module.linkSystemLibrary(opt_link_cblas.?, .{});
    }

    if (opt_link_lapacke != null) {
        exe.root_module.linkSystemLibrary(opt_link_lapacke.?, .{});
    }

    b.installArtifact(exe);

    const run_cmd = b.addRunArtifact(exe);
    run_cmd.step.dependOn(b.getInstallStep());
    const run_step = b.step("run", "Run the executable");
    run_step.dependOn(&run_cmd.step);

    // Compile only CBLAS
    const cblas_lib = b.addLibrary(.{
        .linkage = .dynamic,
        .name = "blas",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/cblas.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    cblas_lib.root_module.addImport("zsl", module);

    const cblas_install = b.addInstallArtifact(cblas_lib, .{});

    const cblas_step = b.step("cblas", "Compile CBLAS library");
    cblas_step.dependOn(&cblas_install.step);

    // Tests
    const opt_verbose_tests = b.option(bool, "verbose_tests", "Enable verbose output for tests") orelse false;
    options.addOption(bool, "verbose_tests", opt_verbose_tests);

    const lib_unit_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("test/zsl.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    lib_unit_tests.root_module.addImport("zsl", module);

    if (opt_link_cblas != null) {
        lib_unit_tests.root_module.linkSystemLibrary(opt_link_cblas.?, .{});
    }

    if (opt_link_lapacke != null) {
        lib_unit_tests.root_module.linkSystemLibrary(opt_link_lapacke.?, .{});
    }

    const run_lib_unit_tests = b.addRunArtifact(lib_unit_tests);
    const test_step = b.step("test", "Run unit tests");
    test_step.dependOn(&run_lib_unit_tests.step);

    // Documentation
    const lib = b.addLibrary(.{
        .name = "zsl",
        .root_module = module,
    });

    const install_docs = b.addInstallDirectory(.{
        .source_dir = lib.getEmittedDocs(),
        .install_dir = .prefix,
        .install_subdir = "docs",
    });

    const docs_step = b.step("docs", "Update documentation in the `docs/` directory");
    docs_step.dependOn(&install_docs.step);

    // Steps
    const check_step = b.step("check", "Check if the code compiles; this is for ZLS");
    check_step.dependOn(&exe.step);
}

/// Default cache size in KB.
fn getDefaultCacheSize(cache_level: usize) usize {
    return switch (cache_level) {
        1 => 32,
        2 => 256,
        3 => 0,
        else => unreachable,
    };
}

/// Attempts to get the system's l{cache_level} cache size, with cache_level 1,
/// 2 or 3, and defaults to 32KB, 256KB or 0B, respectively, if it fails.
fn getCacheSize(b: *std.Build, cache_level: usize) usize {
    if (cache_level == 0 or cache_level > 3)
        return 0;

    const default = getDefaultCacheSize(cache_level);

    return detectCacheSize(b, cache_level) catch |err| {
        std.debug.print(
            "warning: L{d} cache detection failed ({s}), defaulting to {d}KB\n",
            .{ cache_level, @errorName(err), default },
        );

        return default * 1024;
    };
}

fn detectCacheSize(b: *std.Build, cache_level: usize) !usize {
    return switch (builtin.os.tag) {
        .linux => linuxCacheSize(b, cache_level),
        .macos, .ios, .tvos => darwinCacheSize(b, cache_level),
        .windows => windowsCacheSize(b, cache_level),
        else => error.UnsupportedPlatform,
    };
}

/// Walks /sys/devices/system/cpu/cpu0/cache/index* looking for a Data or
/// Unified cache at level 1.
fn linuxCacheSize(b: *std.Build, cache_level: usize) !usize {
    var path_buf: [128]u8 = undefined;

    for (0..16) |cache_idx| {
        // Check cache level
        const level_path = try std.fmt.bufPrint(
            &path_buf,
            "/sys/devices/system/cpu/cpu0/cache/index{d}/level",
            .{cache_idx},
        );

        const level = readSysFsInt(b.graph.io, level_path) catch continue;

        if (level != cache_level)
            continue;

        // Check cache type (Data or Unified)
        const type_path = try std.fmt.bufPrint(
            &path_buf,
            "/sys/devices/system/cpu/cpu0/cache/index{d}/type",
            .{cache_idx},
        );

        const cache_type = try readSysFsStr(b.graph.io, type_path, &path_buf);

        if (std.mem.eql(u8, cache_type, "Instruction"))
            continue;

        // Read size
        const size_path = try std.fmt.bufPrint(
            &path_buf,
            "/sys/devices/system/cpu/cpu0/cache/index{d}/size",
            .{cache_idx},
        );

        const size_str = try readSysFsStr(b.graph.io, size_path, &path_buf);

        return parseKernelCacheSize(size_str);
    }

    return switch (cache_level) {
        1 => error.L1CacheEntryNotFound,
        2 => error.L2CacheEntryNotFound,
        3 => error.L3CacheEntryNotFound,
        else => unreachable,
    };
}

/// Reads a small sysfs text file and trims whitespace.
fn readSysFsStr(io: std.Io, path: []const u8, tmp: []u8) ![]const u8 {
    const file = try std.Io.Dir.openFileAbsolute(io, path, .{});
    defer file.close(io);
    const n = try file.readPositionalAll(io, tmp, 0);
    return std.mem.trim(u8, tmp[0..n], &std.ascii.whitespace);
}

fn readSysFsInt(io: std.Io, path: []const u8) !usize {
    var tmp: [32]u8 = undefined;
    const s = try readSysFsStr(io, path, &tmp);
    return std.fmt.parseInt(usize, s, 10);
}

/// Parses strings like "32K", "512K", "8M" into bytes.
fn parseKernelCacheSize(s: []const u8) !usize {
    if (s.len == 0) return error.EmptyCacheSizeString;
    const suffix = s[s.len - 1];
    const multiplier: usize = switch (suffix) {
        'K', 'k' => 1024,
        'M', 'm' => 1024 * 1024,
        'G', 'g' => 1024 * 1024 * 1024,
        else => 1,
    };
    const digits = if (multiplier != 1) s[0 .. s.len - 1] else s;
    return (try std.fmt.parseInt(usize, digits, 10)) * multiplier;
}

/// `sysctl -n hw.l{n}dcachesize` returns the size in bytes as a plain integer.
fn darwinCacheSize(b: *std.Build, cache_level: usize) !usize {
    const io = b.graph.io;

    const key: []const u8 = switch (cache_level) {
        1 => "hw.l1dcachesize",
        2 => "hw.l2cachesize",
        3 => "hw.l3cachesize",
        else => unreachable,
    };

    var child = try std.process.spawn(io, .{
        .argv = &.{ "sysctl", "-n", key },
        .stdin = .ignore,
        .stdout = .pipe,
        .stderr = .ignore,
    });

    var out_buf: [64]u8 = undefined;
    var reader_buf: [256]u8 = undefined;
    var freader = child.stdout.?.reader(io, &reader_buf);
    const n = try freader.interface.readSliceShort(&out_buf);

    const term = try child.wait(io);
    if (term != .exited or term.exited != 0)
        return error.SysctlFailed;

    const trimmed = std.mem.trim(u8, out_buf[0..n], &std.ascii.whitespace);
    return std.fmt.parseInt(usize, trimmed, 10);
}

/// Win32_CacheMemory.Level: 3 = L1, 4 = L2, 5 = L3. MaxCacheSize is in KB.
fn windowsCacheSize(b: *std.Build, cache_level: usize) !usize {
    const io = b.graph.io;

    const wmi_level: usize = switch (cache_level) {
        1 => 3,
        2 => 4,
        3 => 5,
        else => return error.UnsupportedCacheLevel,
    };

    var level_buf: [2]u8 = undefined;
    const level_str = std.fmt.bufPrint(&level_buf, "{d}", .{wmi_level}) catch unreachable;

    const script = try std.mem.concat(b.allocator, u8, &.{
        "$c = Get-CimInstance Win32_CacheMemory | Where-Object { $_.Level -eq ",
        level_str,
        " } | Select-Object -First 1; if ($c) { $c.MaxCacheSize } else { exit 1 }",
    });
    defer b.allocator.free(script);

    var child = try std.process.spawn(io, .{
        .argv = &.{ "powershell", "-NoProfile", "-Command", script },
        .stdin = .ignore,
        .stdout = .pipe,
        .stderr = .ignore,
        .cwd = .{ .inherit = {} },
    });

    var out_buf: [64]u8 = undefined;
    var reader_buf: [256]u8 = undefined;
    var freader = child.stdout.?.reader(io, &reader_buf);
    const n = try freader.interface.readSliceShort(&out_buf);

    const term = try child.wait(io);
    if (term != .exited or term.exited != 0)
        return error.PowerShellQueryFailed;

    const trimmed = std.mem.trim(u8, out_buf[0..n], &std.ascii.whitespace);
    const kib = try std.fmt.parseInt(usize, trimmed, 10);
    return kib * 1024;
}
