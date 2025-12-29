//! Namespace for symmetric matrix types.

const dense = @import("symmetric/dense.zig");
pub const Dense = dense.Dense;
const banded = @import("symmetric/banded.zig");
pub const Banded = banded.Banded;
const tridiagonal = @import("symmetric/tridiagonal.zig");
pub const Tridiagonal = tridiagonal.Tridiagonal;
const sparse = @import("symmetric/sparse.zig");
pub const Sparse = sparse.Sparse;
const block = @import("symmetric/block.zig");
pub const Block = block.Block;
