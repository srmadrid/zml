// Elementwise generic functions
pub const Apply2 = @import("ops/apply2.zig").Apply2;
pub const apply2 = @import("ops/apply2.zig").apply2;
pub const apply2Unchecked = @import("ops/apply2.zig").apply2Unchecked;
pub const apply2Alloc = @import("ops/apply2.zig").apply2Alloc;
pub const apply2Into = @import("ops/apply2.zig").apply2Into;
pub const apply2IntoUnchecked = @import("ops/apply2.zig").apply2IntoUnchecked;

// Arithmetic operations
pub const Add = @import("ops/add.zig").Add;
pub const add = @import("ops/add.zig").add;
pub const addUnchecked = @import("ops/add.zig").addUnchecked;
pub const addAlloc = @import("ops/add.zig").addAlloc;
pub const addInto = @import("ops/add.zig").addInto;
pub const addIntoUnchecked = @import("ops/add.zig").addIntoUnchecked;
pub const Sub = @import("ops/sub.zig").Sub;
pub const sub = @import("ops/sub.zig").sub;
pub const subUnchecked = @import("ops/sub.zig").subUnchecked;
pub const subAlloc = @import("ops/sub.zig").subAlloc;
pub const subInto = @import("ops/sub.zig").subInto;
pub const subIntoUnchecked = @import("ops/sub.zig").subIntoUnchecked;
pub const Mul = @import("ops/mul.zig").Mul;
pub const mul = @import("ops/mul.zig").mul;
pub const mulAlloc = @import("ops/mul.zig").mulAlloc;
pub const mulInto = @import("ops/mul.zig").mulInto;
pub const mulIntoUnchecked = @import("ops/mul.zig").mulIntoUnchecked;
pub const Div = @import("ops/div.zig").Div;
pub const div = @import("ops/div.zig").div;
pub const divAlloc = @import("ops/div.zig").divAlloc;
pub const divInto = @import("ops/div.zig").divInto;
pub const divIntoUnchecked = @import("ops/div.zig").divIntoUnchecked;
