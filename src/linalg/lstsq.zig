const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const linalg = @import("../linalg.zig");

pub fn LstsqMethod(A: type) type {
    if (!meta.isMatrix(A) or meta.isBuilderMatrix(A))
        @compileError("zsl.linalg.LstsqMethod: A must be a non-builder matrix type, got \n\tA = " ++ @typeName(A) ++ "\n");

    return union(enum) {};
}
