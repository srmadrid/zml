const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");

pub fn Complex(comptime T: type) type {
    if (T != rational.Rational and T != real.Real)
        @compileError("Unsupported type for Complex: " ++ @typeName(T));

    return struct {
        re: T,
        im: T,
        flags: Flags,

        pub const empty: Complex(T) = .{
            .re = .empty,
            .im = .empty,
            .flags = .{ .owns_data = false, .writable = false },
        };

        pub fn deinit(self: *Complex(T), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                self.re.deinit(allocator);
                self.im.deinit(allocator);
            }

            self.* = undefined;
        }

        pub fn copy(self: Complex(T), allocator: std.mem.Allocator) !Complex(T) {
            var re: T = try self.re.copy(allocator);
            errdefer re.deinit(allocator);
            const im: T = try self.im.copy(allocator);

            return .{
                .re = re,
                .im = im,
                .flags = .{ .owns_data = true, .writable = true },
            };
        }
    };
}

// Arithmetic operations
// pub const add = @import("complex/add.zig").add;
// pub const add_ = @import("complex/add_.zig").add_;
// pub const sub = @import("complex/sub.zig").sub;
// pub const sub_ = @import("complex/sub_.zig").sub_;
// pub const mul = @import("complex/mul.zig").mul;
// pub const mul_ = @import("complex/mul_.zig").mul_;
// pub const div = @import("complex/div.zig").div;
// pub const div_ = @import("complex/div_.zig").div_;

// Comparison operations
// pub const eq = @import("complex/eq.zig").eq;
// pub const ne = @import("complex/ne.zig").ne;

// Basic operations
// pub const abs = @import("complex/abs.zig").abs;
pub const neg = @import("complex/neg.zig").neg;

pub const Flags = packed struct {
    owns_data: bool = true,
    writable: bool = true,
};
