// Elementwise generic functions
pub const Apply2 = @import("ops/apply2.zig").Apply2;
pub const apply2 = @import("ops/apply2.zig").apply2;
pub const apply2Alloc = @import("ops/apply2Alloc.zig").apply2Alloc;
pub const apply2Into = @import("ops/apply2Into.zig").apply2Into;

// Arithmetic operations
pub const Add = @import("ops/add.zig").Add;
pub const add = @import("ops/add.zig").add;
pub const addAlloc = @import("ops/addAlloc.zig").addAlloc;
pub const addInto = @import("ops/addInto.zig").addInto;
pub const Sub = @import("ops/sub.zig").Sub;
pub const sub = @import("ops/sub.zig").sub;
pub const subAlloc = @import("ops/subAlloc.zig").subAlloc;
pub const subInto = @import("ops/subInto.zig").subInto;
pub const Mul = @import("ops/mul.zig").Mul;
pub const mul = @import("ops/mul.zig").mul;
pub const mulAlloc = @import("ops/mulAlloc.zig").mulAlloc;
pub const mulInto = @import("ops/mulInto.zig").mulInto;
pub const Div = @import("ops/div.zig").Div;
pub const div = @import("ops/div.zig").div;
pub const divAlloc = @import("ops/divAlloc.zig").divAlloc;
pub const divInto = @import("ops/divInto.zig").divInto;
