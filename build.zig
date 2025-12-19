const std = @import("std");

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
    const opt_max_dimensions = b.option(u32, "max_dimensions", "Maximum number of dimensions for `Array`s") orelse 8;
    options.addOption(u32, "max_dimensions", opt_max_dimensions);

    // Option to provide BLAS and LAPACK implementations
    const opt_link_cblas = b.option([]const u8, "link_cblas", "Link CBLAS implementation");
    options.addOption(?[]const u8, "link_cblas", opt_link_cblas);
    const opt_link_lapacke = b.option([]const u8, "link_lapacke", "Link LAPACKE implementation");
    options.addOption(?[]const u8, "link_lapacke", opt_link_lapacke);

    const module = b.addModule("zml", .{
        .root_source_file = b.path("src/zml.zig"),
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

    exe.root_module.addImport("zml", module);

    if (opt_link_cblas != null or opt_link_lapacke != null) {
        exe.linkLibC();
    }
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

    cblas_lib.root_module.addImport("zml", module);

    const cblas_install = b.addInstallArtifact(cblas_lib, .{});

    const cblas_step = b.step("cblas", "Compile CBLAS library");
    cblas_step.dependOn(&cblas_install.step);

    // Tests
    const opt_verbose_tests = b.option(bool, "verbose_tests", "Enable verbose output for tests") orelse false;
    options.addOption(bool, "verbose_tests", opt_verbose_tests);

    const lib_unit_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("test/zml.zig"),
            .target = target,
            .optimize = optimize,
        }),
    });

    lib_unit_tests.root_module.addImport("zml", module);

    if (opt_link_cblas != null or opt_link_lapacke != null) {
        lib_unit_tests.linkLibC();
    }
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
        .name = "zml",
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
