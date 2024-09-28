const std = @import("std");
const zml = @import("zml.zig");
const Symbol = zml.Symbol;
const SymbolType = @import("symbol.zig").SymbolType;
const Set = zml.Set;
const Expression = zml.Expression;

/// A mathematical variable.
pub const Variable = struct {
    /// Symbol of the variable.
    symbol: Symbol,
    /// Set in which the variable lives.
    domain: *Set,

    /// Initialize a variable.
    pub fn init(allocator: std.mem.Allocator, name: []const u8, domain: *Set) !Variable {
        const symbol = try Symbol.init(name, SymbolType.Variable, &[_]*Symbol{&domain.symbol}, allocator);
        symbol.dependents.append(allocator, &domain.symbol);

        return Variable{
            .symbol = symbol,
            .domain = domain,
        };
    }

    /// Deinitialize a variable.
    pub fn deinit(self: *Variable) !void {
        try self.symbol.deinit();
    }
};
