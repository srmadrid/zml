//! Namespace for hermitian matrix types.

const dense = @import("hermitian/dense.zig");
pub const Dense = dense.Dense;
const banded = @import("hermitian/banded.zig");
pub const Banded = banded.Banded;
const tridiagonal = @import("hermitian/tridiagonal.zig");
pub const Tridiagonal = tridiagonal.Tridiagonal;
const sparse = @import("hermitian/sparse.zig");
pub const Sparse = sparse.Sparse;
const block = @import("hermitian/block.zig");
pub const Block = block.Block;
