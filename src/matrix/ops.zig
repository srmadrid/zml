const std = @import("std");

// Maybe in place versions (_) also require allocator if the coerced type of the inputs is arbitrary precision

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const EnsureMatrix = types.EnsureMatrix;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const matrix = @import("../array.zig");

const general = @import("general.zig");
const symmetric = @import("symmetric.zig");
const hermitian = @import("hermitian.zig");
const triangular = @import("triangular.zig");
const diagonal = @import("diagonal.zig");
const banded = @import("banded.zig");
const tridiagonal = @import("tridiagonal.zig");
const sparse = @import("sparse.zig");
