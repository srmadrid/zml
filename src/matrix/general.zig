//! Namespace for general matrix types.

const dense = @import("general/dense.zig");
pub const Dense = dense.Dense;
const sparse = @import("general/sparse.zig");
pub const Sparse = sparse.Sparse;
const block = @import("general/block.zig");
pub const Block = block.Block;
