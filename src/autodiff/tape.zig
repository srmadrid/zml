const std = @import("std");

const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const autodiff = @import("../autodiff.zig");

pub fn Tape(N: type) type {
    if (comptime !meta.isNumeric(N))
        @compileError("zsl.autodiff.Tape: N must be a numeric type, got \n\tN: " ++ @typeName(N) ++ "\n");

    return struct {
        nodes: [*]Node,
        len: usize,
        _nlen: usize,

        pub const Node = struct {
            op: autodiff.Op,
            left: usize,
            right: usize,
            val: N,
            grad: N = numeric.cast(N, 0),
        };

        pub const Numeric = N;

        pub fn init(allocator: std.mem.Allocator, capacity: usize) !autodiff.Tape(N) {
            return .{
                .nodes = (try allocator.alloc(Node, capacity)).ptr,
                .len = 0,
                ._nlen = capacity,
            };
        }

        pub fn deinit(self: *autodiff.Tape(N), allocator: std.mem.Allocator) void {
            allocator.free(self.nodes[0..self._nlen]);

            self.* = undefined;
        }

        pub fn pushAssumeCapacity(self: *autodiff.Tape(N), node: autodiff.Tape(N).Node) usize {
            const id = self.len;
            self.nodes[id] = node;
            self.len += 1;
            return id;
        }

        pub fn clear(self: *autodiff.Tape(N)) void {
            self.len = 0;
        }
    };
}
