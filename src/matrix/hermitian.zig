//! Namespace for hermitian matrix types.

const dense = @import("hermitian/dense.zig");
pub const Dense = dense.Dense;
const sparse = @import("hermitian/sparse.zig");
pub const Sparse = sparse.Sparse;
