const std = @import("std");
const zml = @import("../zml.zig");
const Symbol = zml.Symbol;

/// A
pub const ExpressionTree = struct {
    root: Node,

    pub const Node = struct {
        value: Symbol,
        parent: ?*Node,
        left: ?*Node,
        right: ?*Node,

        fn init(value: Symbol) Node {
            return Node{
                .value = value,
                .parent = null,
                .left = null,
                .right = null,
            };
        }
    };

    pub fn init(value: Symbol) ExpressionTree {
        return ExpressionTree{
            .root = Node.init(value),
        };
    }
};
