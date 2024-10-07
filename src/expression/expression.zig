const std = @import("std");
const zml = @import("../zml.zig");
const Symbol = zml.Symbol;
const SymbolType = @import("../symbol.zig").SymbolType;
const ExpressionTree = @import("exprtree.zig").ExpressionTree;

/// A mathematical expression.
pub const Expression = struct {
    /// String representation of the expression.
    string: []const u8,
    /// Tree representation of the expression.
    tree: ExpressionTree,
    /// Dependencies.
    dependencies: std.ArrayListUnmanaged(*Symbol),
    /// Symbols used only inside the expression.
    intdependencies: std.ArrayListUnmanaged(Symbol),
    /// Allocator for the expression.
    allocator: std.mem.Allocator,

    /// Initialize an expression.
    pub fn init(allocator: std.mem.Allocator, expression: []const u8, dependencies: []const *Symbol) !Expression {
        const dp = std.ArrayListUnmanaged(*Symbol).initBuffer(dependencies);
        const idp = try std.ArrayListUnmanaged(Symbol).initCapacity(allocator, 2);

        return Expression{
            .string = expression,
            .tree = ExpressionTree{ .root = null },
            .dependencies = dp,
            .intdependencies = idp,
            .allocator = allocator,
        };
    }

    /// Tokenizes an input string into an ArrayList of Symbols.
    pub fn tokenize(allocator: std.mem.Allocator, expression: []const u8) !std.ArrayListUnmanaged(Token) {
        var array = try std.ArrayListUnmanaged(Token).initCapacity(allocator, 2);
        var stack = try std.ArrayListUnmanaged(TokenType).initCapacity(allocator, 2);

        var i: usize = 0;
        while (i < expression.len) {
            std.debug.print("{d}, ", .{i});
            // \{x\in\mathbb{R}\mid x > 0, \sin(x) = 0\}
            const stacklen = stack.items.len;

            if (expression[i] == ' ') {
                i += 1;
                continue;
            }

            if (expression[i] == ',') {
                try array.append(allocator, Token{ .string = expression[i .. i + 1], .type = TokenType.Comma });
                i += 1;
                continue;
            }

            if (expression[i] == '(') {
                try array.append(allocator, Token{ .string = expression[i .. i + 1], .type = TokenType.OpeningParenthesis });
                try stack.append(allocator, TokenType.OpeningParenthesis);
                i += 1;
                continue;
            }

            var j: usize = i;
            if (expression[i] == '\\') {
                if (i + 1 >= expression.len) {
                    array.deinit(allocator);
                    stack.deinit(allocator);
                    return error.Untokenizable;
                }

                if (expression[i + 1] == '{') {
                    try array.append(allocator, Token{ .string = expression[i .. i + 2], .type = TokenType.OpeningBrackets });
                    i += 2;
                    continue;
                }

                if (expression[i + 1] == '}') {
                    try array.append(allocator, Token{ .string = expression[i .. i + 2], .type = TokenType.ClosingBrackets });
                    i += 2;
                    continue;
                }

                j += 1;
            }

            while (j < expression.len and expression[j] != ' ' and expression[j] != ',' and expression[j] != '\\') {
                if (expression[j] == '(') {
                    try stack.append(allocator, TokenType.OpeningParenthesis);
                } else if (expression[j] == '{') {
                    try stack.append(allocator, TokenType.OpeningBrackets);
                } else if (expression[j] == '}') {
                    if (stack.pop() != TokenType.OpeningBrackets) {
                        array.deinit(allocator);
                        stack.deinit(allocator);
                        return error.Untokenizable;
                    }
                } else if (expression[j] == ')') {
                    if (stack.pop() != TokenType.OpeningParenthesis) {
                        array.deinit(allocator);
                        stack.deinit(allocator);
                        return error.Untokenizable;
                    }
                }

                j += 1;
            }

            try array.append(allocator, Token{ .string = expression[i..j], .type = TokenType.Symbol });

            if (stack.items.len != stacklen) {
                array.deinit(allocator);
                stack.deinit(allocator);
                return error.Untokenizable;
            }

            i = j;
        }

        std.debug.print("\n", .{});
        stack.deinit(allocator);
        return array;
    }
};

pub const Token = struct {
    string: []const u8,
    type: TokenType,
};

pub const TokenType = enum {
    Symbol,
    OpeningParenthesis,
    ClosingParenthesis,
    OpeningBrackets,
    ClosingBrackets,
    Comma,
};

test {
    _ = @import("exprtree.zig");
    std.testing.refAllDeclsRecursive(@This());
}
