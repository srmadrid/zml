//! Namespace for matrix builder types.

const sparse = @import("builder/sparse.zig");
pub const Sparse = sparse.Sparse;
const block = @import("builder/block.zig");
pub const Block = block.Block;
