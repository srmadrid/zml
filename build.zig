const std = @import("std");

pub fn build(b: *std.Build) void {
    const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    // Option to provide BLAS and LAPACK implementations
    const options = b.addOptions();
    const opt_link_cblas = b.option([]const u8, "link_cblas", "Link CBLAS implementation");
    options.addOption(?[]const u8, "link_cblas", opt_link_cblas);
    const opt_link_lapacke = b.option([]const u8, "link_lapacke", "Link LAPACKE implementation");
    options.addOption(?[]const u8, "link_apacke", opt_link_lapacke);

    const module = b.addModule("zml", .{
        .root_source_file = b.path("src/zml.zig"),
    });
    module.addOptions("options", options);

    // Executable (for testing)
    const exe = b.addExecutable(.{
        .name = "main",
        .root_source_file = b.path("src/main.zig"),
        .target = target,
        .optimize = optimize,
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

    // Tests
    const lib_unit_tests = b.addTest(.{
        .root_source_file = b.path("test/zml.zig"),
        .target = target,
        .optimize = optimize,
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

    // Examples
    const example_step = b.step("examples", "Build examples");
    for ([_][]const u8{}) |example_name| {
        const example = b.addExecutable(.{
            .name = example_name,
            .root_source_file = b.path(b.fmt("examples/{s}.zig", .{example_name})),
            .target = target,
            .optimize = optimize,
        });
        const install_example = b.addInstallArtifact(example, .{});
        example.root_module.addImport("zml", module);
        example_step.dependOn(&example.step);
        example_step.dependOn(&install_example.step);
    }

    // Steps
    const check_step = b.step("check", "Check if the code compiles; this is for ZLS");
    check_step.dependOn(&exe.step);
}
