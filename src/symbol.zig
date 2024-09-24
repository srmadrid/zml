const std = @import("std");

/// The basic type for ZML.
pub const Symbol = struct {
    /// String representation of the symbol.
    name: []const u8,
    /// Unique identification number of the symbol.
    id: usize,
    /// Type of the symbol.
    type: SymbolType,
    /// Pointer to the mathematical object.
    object: *anyopaque,
    /// Dependent symbols.
    dependents: std.ArrayListUnmanaged(*Symbol),
    /// Dependencies.
    dependencies: std.ArrayListUnmanaged(*Symbol),
    /// Allocator for the symbol.
    allocator: ?std.mem.Allocator,
};

/// Types of symbols.
pub const SymbolType = enum {
    /// An element of a set. Fixed.
    Element,
    /// A variable of a set. Can be any element of a set.
    Variable,
    /// Mathematical set.
    Set,
};
