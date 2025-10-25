const std = @import("std");

const types = @import("types.zig");
const ops = @import("ops.zig");
const constants = @import("constants.zig");
const int = @import("int.zig");
const float = @import("float.zig");
const integer = @import("integer.zig");
const Integer = integer.Integer;
const complex = @import("complex.zig");
const Complex = complex.Complex;

pub const Rational = struct {
    num: Integer,
    den: Integer,
    flags: Flags,

    pub const empty: Rational = .{
        .num = .empty,
        .den = .empty,
        .flags = .{ .owns_data = false, .writable = false },
    };

    /// Initializes a new `Rational` with the specified numerator and
    /// denominator sizes.
    ///
    /// Sizes of zero are allowed, resulting in a numerator or denominator with
    /// no allocated limbs.
    ///
    /// Parameters
    /// ----------
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations.
    ///
    /// `numsize` (`u32`):
    /// The size of the numerator in limbs.
    ///
    /// `densize` (`u32`):
    /// The size of the denominator in limbs.
    ///
    /// Returns
    /// -------
    /// `Rational`:
    /// The newly initialized `Rational`.
    ///
    /// Errors
    /// ------
    /// `std.mem.Allocator.Error.OutOfMemory`:
    /// If memory allocation fails.
    pub fn init(allocator: std.mem.Allocator, numsize: u32, densize: u32) !Rational {
        var num: Integer = try Integer.init(allocator, numsize);
        errdefer num.deinit(allocator);

        var den: Integer = try Integer.init(allocator, densize);
        errdefer den.deinit(allocator);

        return .{
            .num = num,
            .den = den,
            .flags = .{ .owns_data = true, .writable = true },
        };
    }

    /// Initializes a new `Rational` with the specified numerator and
    /// denominator.
    ///
    /// Signature
    /// ---------
    /// ```zig
    /// fn initSet(allocator: std.mem.Allocator, numerator: N, denominator: D) !Rational
    /// ```
    ///
    /// Parameters
    /// ----------
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations.
    ///
    /// `numerator` (`anytype`):
    /// The value to set the numerator to.
    ///
    /// `denominator` (`anytype`):
    /// The value to set the denominator to.
    ///
    /// Returns
    /// -------
    /// `Rational`:
    /// The newly initialized `Rational`.
    ///
    /// Errors
    /// ------
    /// `std.mem.Allocator.Error.OutOfMemory`:
    /// If memory allocation fails.
    ///
    /// `rational.Error.ZeroDenominator`:
    /// If the denominator is zero.
    pub fn initSet(allocator: std.mem.Allocator, numerator: anytype, denominator: anytype) !Rational {
        // Edit to properly handle floats, cfloats, complexes and other rationals correctly.
        // If two complex numbers are given, perform complex division before, or just ignore im part (integer just ignores it, but noe division is ever involved there)?
        if (ops.eq(denominator, 0, .{}) catch unreachable)
            return Error.ZeroDenominator;

        var num: Integer = try Integer.initSet(allocator, numerator);
        errdefer num.deinit(allocator);

        var den: Integer = try Integer.initSet(allocator, denominator);
        errdefer den.deinit(allocator);

        var r: Rational = .{
            .num = num,
            .den = den,
            .flags = .{ .owns_data = true, .writable = true },
        };

        try r.reduce(allocator);

        return r;
    }

    /// Deinitializes the `Rational`, freeing any allocated memory.
    ///
    /// If the `Rational` does not own its data, no memory is freed and this
    /// becomes a no-op.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*Rational`):
    /// A pointer to the `Rational` to deinitialize.
    ///
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory deallocation. Must be the same allocator
    /// used to initialize `self`.
    ///
    /// Returns
    /// -------
    /// `void`
    pub fn deinit(self: *Rational, allocator: std.mem.Allocator) void {
        if (self.flags.owns_data) {
            self.num.deinit(allocator);
            self.den.deinit(allocator);
        }

        self.* = undefined;
    }

    /// Reduces the `Rational` to its simplest form, ensuring that the numerator
    /// and denominator are coprime.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*Rational`):
    /// A pointer to the `Rational` to reduce.
    ///
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations. Must be the same allocator
    /// used to initialize `self`.
    ///
    /// Returns
    /// -------
    /// `void`
    ///
    /// Errors
    /// ------
    /// `std.mem.Allocator.Error.OutOfMemory`:
    /// If memory allocation fails.
    ///
    /// `rational.Error.ZeroDenominator`:
    /// If the denominator is zero.
    pub fn reduce(self: *Rational, allocator: std.mem.Allocator) !void {
        if (!self.flags.writable)
            return Error.NotWritable;

        var g: Integer = try integer.gcd(allocator, self.num, self.den);
        defer g.deinit(allocator);

        try integer.div_(allocator, &self.num, self.num, g);
        try integer.div_(allocator, &self.den, self.den, g);
    }

    /// Creates a copy of the `Rational`.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to copy.
    ///
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations.
    ///
    /// Returns
    /// -------
    /// `Rational`:
    /// The newly created copy of the `Rational`.
    ///
    /// Errors
    /// ------
    /// `std.mem.Allocator.Error.OutOfMemory`:
    /// If memory allocation fails.
    pub fn copy(self: *const Rational, allocator: std.mem.Allocator) !Rational {
        var num: Integer = try self.num.copy(allocator);
        errdefer num.deinit(allocator);
        const den: Integer = try self.den.copy(allocator);

        return .{
            .num = num,
            .den = den,
            .flags = .{ .owns_data = true, .writable = true },
        };
    }

    /// Converts the `Rational` to an integer type `Int`, performing truncation
    /// if necessary.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to convert.
    ///
    /// `Int` (`type`):
    /// The integer type to convert to.
    ///
    /// Returns
    /// -------
    /// `Int`:
    /// The converted integer value.
    pub fn toInt(self: *const Rational, comptime Int: type) Int {
        comptime if (types.numericType(Int) != .int)
            @compileError("rational.toInt requires Int to be an int type, got " ++ @typeName(Int));

        const num_f: f128 = self.num.toFloat(f128);
        const den_f: f128 = self.den.toFloat(f128);

        const q: f128 = num_f / den_f;
        if (@typeInfo(Int).int.signedness == .unsigned) {
            if (q < 0)
                return 0;

            if (float.gt(q, int.maxVal(Int)))
                return int.maxVal(Int);
        } else {
            if (float.lt(q, int.minVal(Int)))
                return int.minVal(Int);

            if (float.gt(q, int.maxVal(Int)))
                return int.maxVal(Int);
        }

        return types.scast(Int, q);
    }

    /// Converts the `Rational` to a floating-point type `Float`.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to convert.
    ///
    /// `Float` (`type`):
    /// The floating-point type to convert to.
    ///
    /// Returns
    /// -------
    /// `Float`:
    /// The converted floating-point value.
    pub fn toFloat(self: *const Rational, comptime Float: type) Float {
        comptime if (types.numericType(Float) != .float)
            @compileError("rational.toFloat requires Float to be a float type, got " ++ @typeName(Float));

        const num_f: f128 = self.num.toFloat(f128);
        const den_f: f128 = self.den.toFloat(f128);

        return types.scast(Float, num_f / den_f);
    }

    /// Returns a view of the `Rational` as a complex number with zero imaginary
    /// part.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to view as a complex number.
    ///
    /// Returns
    /// -------
    /// `Complex(Rational)`:
    /// The `Rational` viewed as a complex number.
    pub fn asComplex(self: *const Rational) Complex(Rational) {
        var re: Rational = self.*;
        re.flags.owns_data = false;

        return .{
            .re = re,
            .im = constants.zero(Rational, .{}) catch unreachable,
            .flags = .{ .owns_data = false, .writable = false },
        };
    }

    /// Converts the `Rational` into a complex number with zero imaginary part.
    /// After the conversion, the new complex number takes ownership of the data
    /// and the original `Rational` is invalidated.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*Rational`):
    /// A pointer to the `Rational` to convert.
    ///
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations.
    ///
    /// Returns
    /// -------
    /// `Complex(Rational)`:
    /// The newly created complex number.
    ///
    /// Errors
    /// ------
    /// `std.mem.Allocator.Error.OutOfMemory`:
    /// If memory allocation fails.
    ///
    /// `rational.Error.DataNotOwned`:
    /// If the `Rational` does not own its data.
    pub fn toComplex(self: *Rational, allocator: std.mem.Allocator) !Complex(Rational) {
        if (!self.flags.owns_data)
            return Error.DataNotOwned;

        const result: Complex(Rational) = .{
            .re = self,
            .im = try constants.zero(Integer, .{ .allocator = allocator }),
            .flags = .{ .owns_data = true, .writable = self.flags.writable },
        };

        self.* = undefined;

        return result;
    }

    /// Creates a copy of the `Rational` as a complex number with zero imaginary
    /// part.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to copy.
    ///
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations.
    ///
    /// Returns
    /// -------
    /// `Complex(Rational)`:
    /// The newly created complex number.
    ///
    /// Errors
    /// ------
    /// `std.mem.Allocator.Error.OutOfMemory`:
    /// If memory allocation fails.
    pub fn copyToComplex(self: *const Rational, allocator: std.mem.Allocator) !Complex(Rational) {
        var re: Rational = try self.copy(allocator);
        errdefer re.deinit(allocator);

        return .{
            .re = re,
            .im = try constants.zero(Rational, .{ .allocator = allocator }),
            .flags = .{ .owns_data = true, .writable = true },
        };
    }
};

// Arithmetic operations
pub const add = @import("rational/add.zig").add;
pub const add_ = @import("rational/add_.zig").add_;
pub const sub = @import("rational/sub.zig").sub;
pub const sub_ = @import("rational/sub_.zig").sub_;
pub const mul = @import("rational/mul.zig").mul;
pub const mul_ = @import("rational/mul_.zig").mul_;
pub const div = @import("rational/div.zig").div;
pub const div_ = @import("rational/div_.zig").div_;

// Comparison operations
// pub const cmp = @import("rational/cmp.zig").cmp;
// pub const eq = @import("rational/eq.zig").eq;
// pub const ne = @import("rational/ne.zig").ne;
// pub const lt = @import("rational/lt.zig").lt;
// pub const le = @import("rational/le.zig").le;
// pub const gt = @import("rational/gt.zig").gt;
// pub const ge = @import("rational/ge.zig").ge;

// Basic operations
pub const abs = @import("rational/abs.zig").abs;
pub const neg = @import("rational/neg.zig").neg;

pub const Error = error{
    ZeroDenominator,
    ZeroDivision,
    NotWritable,
};

pub const Flags = packed struct {
    owns_data: bool = true,
    writable: bool = true,
};
