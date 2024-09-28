const std = @import("std");

/// The basic type for ZML.
pub const Symbol = struct {
    /// String representation of the symbol.
    name: []const u8,
    /// Unique identification number of the symbol.
    id: usize,
    /// Type of the symbol.
    type: SymbolType,
    /// Dependencies.
    dependencies: std.ArrayListUnmanaged(*Symbol),
    /// Dependent symbols.
    dependents: std.ArrayListUnmanaged(*Symbol),
    /// Allocator for the symbol.
    allocator: std.mem.Allocator,

    /// Initialize a symbol.
    pub fn init(name: []const u8, symbol_type: SymbolType, dependencies: []const *Symbol, allocator: std.mem.Allocator) !Symbol {
        const _dependencies = try std.ArrayListUnmanaged(*Symbol).initCapacity(allocator, @min(2, dependencies.len));
        const _dependents = try std.ArrayListUnmanaged(*Symbol).initCapacity(allocator, 2);

        return Symbol{
            .name = name,
            .id = generateID(),
            .type = symbol_type,
            .dependencies = _dependencies,
            .dependents = _dependents,
            .allocator = allocator,
        };
    }

    /// Deinitialize a symbol.
    pub fn deinit(self: *Symbol) !void {
        if (self.dependents.items.len != 0) {
            return error.Error;
        }

        self.dependencies.deinit(self.allocator);
        self.dependents.deinit(self.allocator);
    }

    /// Returns a pointer to the parent object.
    pub fn object(self: *Symbol) anyopaque {
        return @fieldParentPtr("symbol", self);
    }
};

/// Generate a random unused ID.
pub fn generateID() usize {
    const number = std.crypto.random.int(usize);
    // Somehow check if ID is in use
    return number;
}

/// Types of symbols.
pub const SymbolType = enum {
    /// An element of a set. Fixed.
    Element,
    /// A variable of a set. Can be any element of a set.
    Variable,
    /// Mathematical set.
    Set,
};
