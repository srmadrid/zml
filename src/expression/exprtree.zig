const std = @import("std");
const zml = @import("../zml.zig");
const Symbol = zml.Symbol;

/// A tree containing an expression.
pub const ExpressionTree = struct {
    root: Node,

    pub const Node = struct {
        value: *Symbol,
        parent: ?*Node,
        children: std.ArrayListUnmanaged(*Node),

        /// Initialize with the input symbol.
        fn init(allocator: std.mem.Allocator, value: *Symbol, parent: ?*Node) !Node {
            return Node{
                .value = value,
                .parent = parent,
                .children = try std.ArrayListUnmanaged(*Node).initCapacity(allocator, 2),
            };
        }

        /// Release all the memory of the node and its children.
        fn deinit(self: *Node, allocator: std.mem.Allocator) void {
            for (self.children.items) |child| {
                child.deinit(allocator);
            }

            self.children.deinit(allocator);
        }
    };

    /// Initialize the tree with the value as the root.
    pub fn init(allocator: std.mem.Allocator, value: Symbol) !ExpressionTree {
        return ExpressionTree{
            .root = try Node.init(allocator, value, null),
        };
    }

    /// Deinitialize the tree and release all of its memory.
    pub fn deinit(self: *ExpressionTree, allocator: std.mem.Allocator) void {
        if (self.root) {
            self.root.deinit(allocator);
        }
    }
};
