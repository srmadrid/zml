//! Namespace for permutation matrix types.

pub const Static = @import("permutation/static.zig").Static;
pub const Sparse = @import("permutation/sparse.zig").Sparse;

pub const Direction = enum {
    forward,
    backward,

    pub fn invert(self: Direction) Direction {
        return switch (self) {
            .forward => .backward,
            .backward => .forward,
        };
    }
};
