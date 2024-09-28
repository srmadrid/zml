const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const SymbolType = @import("symbol.zig").SymbolType;
const Expression = zml.Expression;
const Variable = zml.Variable;

/// A mathematical set.
pub const Set = struct {
    /// Symbol of the set.
    symbol: Symbol,
    /// Variable from which the set is created. The domain of the set is
    /// contained here, and if it is equal to the set itself, and no condition
    /// is passed, then the set is only declared, i.e., it is equivalent to
    /// just saying ``Let $G$ be a set'' without explicitly defining $G$ itself.
    variable: Variable,
    /// Conditions that the elements satisfy.
    condition: Expression,

    /// Initialization functions
    pub const init = struct {
        /// Initialize a set using set-builder style.
        pub fn builder(allocator: std.mem.Allocator, name: []const u8, variable: Variable, conditions: []const Expression) !Set {
            // Some function that joins predicates with ANDs.
            const expr = Expression.combineAND(allocator, conditions);

            const symbol = try Symbol.init(name, SymbolType.Set, &[_]*Symbol{&variable.symbol}, allocator);
            symbol.dependents.append(allocator, &variable.symbol);

            return Set{
                .symbol = symbol,
                .variable = variable,
                .condition = expr,
            };
        }
    };
};
