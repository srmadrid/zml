//! Namespace for symmetric matrix types.

const dense = @import("symmetric/dense.zig");
pub const Dense = dense.Dense;
const sparse = @import("symmetric/sparse.zig");
pub const Sparse = sparse.Sparse;
const block = @import("symmetric/block.zig");
pub const Block = block.Block;
