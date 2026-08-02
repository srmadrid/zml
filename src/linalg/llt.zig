pub const Static = @import("llt/static.zig").Static;
pub const Dense = @import("llt/dense.zig").Dense;
// pub const Sparse = @import("llt/sparse.zig");

pub const Factor = @import("llt/ops.zig").Factor;
pub const factor = @import("llt/ops.zig").factor;
pub const factorAlloc = @import("llt/ops.zig").factorAlloc;
pub const factorInto = @import("llt/ops.zig").factorInto;
pub const factorIntoUnchecked = @import("llt/ops.zig").factorIntoUnchecked;
