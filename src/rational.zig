//! Namespace for rational operations.

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

pub var default_accuracy: u32 = 50;
pub var default_internal_accuracy: u32 = 60;

/// Arbotrary-precision rational type, represented as a fraction of two
/// arbitrary-precision integers (`Integer`).
pub const Rational = struct {
    num: Integer,
    den: Integer,
    flags: Flags,

    /// Type signature
    pub const is_rational = {};
    pub const is_real_type = {};
    pub const is_signed = {};
    pub const is_allocated = true;

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
    /// The value to set the numerator to. Must be a numeric type or a string.
    ///
    /// `denominator` (`anytype`):
    /// The value to set the denominator to. Must be a numeric type or a string.
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
        var r: Rational = try Rational.init(allocator, 0, 0);
        errdefer r.deinit(allocator);

        try r.set(allocator, numerator, denominator);
        return r;
    }

    /// Deinitializes the `Rational`, freeing any allocated memory and
    /// invalidating it.
    ///
    /// If the `Rational` does not own its data, no memory is freed and this
    /// only invalidates it.
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

    /// Sets the value of the `Rational`. If either `numerator` or `denominator`
    /// is not an integer type, a full division is performed to set the value.
    ///
    /// Signature
    /// ---------
    /// ```zig
    /// fn set(allocator: std.mem.Allocator, numerator: N, denominator: D) !void
    /// ```
    ///
    /// Parameters
    /// ----------
    /// `self` (`*Rational`):
    /// A pointer to the `Rational` to set.
    ///
    /// `allocator` (`std.mem.Allocator`):
    /// The allocator to use for memory allocations. Must be the same allocator
    /// used to initialize `self`.
    ///
    /// `numerator` (`anytype`):
    /// The value to set the numerator to. Must be a numeric type or a string.
    /// Complex values use their real part.
    ///
    /// `denominator` (`anytype`):
    /// The value to set the denominator to. Must be a numeric type or a string.
    /// Complex values use their real part.
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
    /// `rational.Error.NotWritable`:
    /// If the `Rational` is not writable.
    ///
    /// `rational.Error.NotFinite`:
    /// If the value is a float or complex number that is not finite.
    ///
    /// `rational.Error.ZeroDenominator`:
    /// If the denominator is zero.
    pub fn set(self: *Rational, allocator: std.mem.Allocator, numerator: anytype, denominator: anytype) !void {
        const N: type = @TypeOf(numerator);
        const D: type = @TypeOf(denominator);

        if (!self.flags.writable)
            return Error.NotWritable;

        if (comptime types.isNumeric(N) and types.isNumeric(D)) {
            switch (comptime types.numericType(N)) {
                .bool => switch (comptime types.numericType(D)) {
                    .bool => {
                        if (!denominator)
                            return Error.ZeroDenominator;

                        try self.num.reserve(allocator, 1);
                        try self.den.reserve(allocator, 1);

                        if (numerator) {
                            self.num.limbs[0] = 1;
                            self.num.size = 1;
                            self.num.positive = true;
                        } else {
                            self.num.size = 0;
                            self.num.positive = true;
                        }

                        self.den.limbs[0] = 1;
                        self.den.size = 1;
                    },
                    .int => {
                        if (denominator == 0)
                            return Error.ZeroDenominator;

                        if (!numerator) {
                            try self.num.reserve(allocator, 1);
                            self.num.size = 0;
                            self.num.positive = true;

                            try self.den.reserve(allocator, 1);
                            self.den.limbs[0] = 1;
                            self.den.size = 1;
                        } else {
                            var dvalue = @import("int/asInteger.zig").asInteger(denominator);
                            dvalue[0].limbs = &dvalue[1];

                            try self.num.set(allocator, 1);
                            try self.den.set(allocator, dvalue[0]);
                        }
                    },
                    .float => {
                        if (denominator == 0.0)
                            return Error.ZeroDenominator;

                        if (!numerator) {
                            try self.num.reserve(allocator, 1);
                            self.num.size = 0;
                            self.num.positive = true;

                            try self.den.reserve(allocator, 1);
                            self.den.limbs[0] = 1;
                            self.den.size = 1;
                        } else {
                            var dvalue = try @import("float/asRational.zig").asRational(denominator);
                            dvalue[0].num.limbs = &dvalue[1][0];
                            dvalue[0].den.limbs = &dvalue[1][1];

                            try self.num.set(allocator, dvalue[0].den);
                            try self.den.set(allocator, dvalue[0].num);
                        }
                    },
                    .cfloat => return self.set(allocator, numerator, denominator.re),
                    .integer => {
                        if (denominator.size == 0)
                            return Error.ZeroDenominator;

                        if (!numerator) {
                            try self.num.reserve(allocator, 1);
                            self.num.size = 0;
                            self.num.positive = true;

                            try self.den.set(allocator, denominator);
                        } else {
                            try self.num.set(allocator, 1);
                            try self.den.set(allocator, denominator);
                        }
                    },
                    .rational => {
                        if (denominator.num.size == 0)
                            return Error.ZeroDenominator;

                        if (!numerator) {
                            try self.num.reserve(allocator, 1);
                            self.num.size = 0;
                            self.num.positive = true;

                            try self.den.set(allocator, denominator);
                        } else {
                            try self.num.set(allocator, denominator.den);
                            try self.den.set(allocator, denominator.num);
                        }
                    },
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator, denominator.re),
                },
                .int => switch (comptime types.numericType(D)) {
                    .bool => {
                        if (!denominator)
                            return Error.ZeroDenominator;

                        var nvalue = @import("int/asInteger.zig").asInteger(numerator);
                        nvalue[0].limbs = &nvalue[1];

                        try self.num.set(allocator, nvalue[0]);
                        try self.den.set(allocator, 1);
                    },
                    .int => {
                        if (denominator == 0)
                            return Error.ZeroDenominator;

                        var nvalue = @import("int/asInteger.zig").asInteger(numerator);
                        nvalue[0].limbs = &nvalue[1];

                        var dvalue = @import("int/asInteger.zig").asInteger(denominator);
                        dvalue[0].limbs = &dvalue[1];

                        try self.num.set(allocator, nvalue[0]);
                        try self.den.set(allocator, dvalue[0]);
                    },
                    .float => {
                        if (denominator == 0.0)
                            return Error.ZeroDenominator;

                        var nvalue = @import("int/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1];

                        var dvalue = try @import("float/asRational.zig").asRational(denominator);
                        dvalue[0].num.limbs = &dvalue[1][0];
                        dvalue[0].den.limbs = &dvalue[1][1];

                        try div_(allocator, self, nvalue[0], dvalue[0]);
                    },
                    .cfloat => return self.set(allocator, numerator, denominator.re),
                    .integer => {
                        if (denominator.size == 0)
                            return Error.ZeroDenominator;

                        var nvalue = @import("int/asInteger.zig").asInteger(numerator);
                        nvalue[0].limbs = &nvalue[1];

                        try self.num.set(allocator, nvalue[0]);
                        try self.den.set(allocator, denominator);
                    },
                    .rational => {
                        if (denominator.num.size == 0)
                            return Error.ZeroDenominator;

                        var nvalue = @import("int/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1];

                        try div_(allocator, self, nvalue[0], denominator);
                    },
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator, denominator.re),
                },
                .float => switch (comptime types.numericType(D)) {
                    .bool => {
                        if (!denominator)
                            return Error.ZeroDenominator;

                        var nvalue = try @import("float/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1][0];
                        nvalue[0].den.limbs = &nvalue[1][1];

                        try self.num.set(allocator, nvalue[0].num);
                        try self.den.set(allocator, nvalue[0].den);
                    },
                    .int => {
                        if (denominator == 0)
                            return Error.ZeroDenominator;

                        var nvalue = try @import("float/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1][0];
                        nvalue[0].den.limbs = &nvalue[1][1];

                        var dvalue = @import("int/asRational.zig").asRational(denominator);
                        dvalue[0].num.limbs = &dvalue[1];

                        try div_(allocator, self, nvalue[0], dvalue[0]);
                    },
                    .float => {
                        if (denominator == 0.0)
                            return Error.ZeroDenominator;

                        var nvalue = try @import("float/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1][0];
                        nvalue[0].den.limbs = &nvalue[1][1];

                        var dvalue = try @import("float/asRational.zig").asRational(denominator);
                        dvalue[0].num.limbs = &dvalue[1][0];
                        dvalue[0].den.limbs = &dvalue[1][1];

                        try div_(allocator, self, nvalue[0], dvalue[0]);
                    },
                    .cfloat => return self.set(allocator, numerator, denominator.re),
                    .integer => {
                        if (denominator.size == 0)
                            return Error.ZeroDenominator;

                        var nvalue = try @import("float/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1][0];
                        nvalue[0].den.limbs = &nvalue[1][1];

                        try div_(allocator, self, nvalue[0], denominator.asRational());
                    },
                    .rational => {
                        if (denominator.num.size == 0)
                            return Error.ZeroDenominator;

                        var nvalue = try @import("float/asRational.zig").asRational(numerator);
                        nvalue[0].num.limbs = &nvalue[1][0];
                        nvalue[0].den.limbs = &nvalue[1][1];

                        try div_(allocator, self, nvalue[0], denominator);
                    },
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator, denominator.re),
                },
                .cfloat => switch (comptime types.numericType(D)) {
                    .bool => return self.set(allocator, numerator.re, denominator),
                    .int => return self.set(allocator, numerator.re, denominator),
                    .float => return self.set(allocator, numerator.re, denominator),
                    .cfloat => return self.set(allocator, numerator.re, denominator.re),
                    .integer => return self.set(allocator, numerator.re, denominator),
                    .rational => return self.set(allocator, numerator.re, denominator),
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator.re, denominator.re),
                },
                .integer => switch (comptime types.numericType(D)) {
                    .bool => {
                        if (!denominator)
                            return Error.ZeroDenominator;

                        try self.num.set(allocator, numerator);
                        try self.den.set(allocator, 1);
                    },
                    .int => {
                        if (denominator == 0)
                            return Error.ZeroDenominator;

                        try self.num.set(allocator, numerator);
                        try self.den.set(allocator, denominator);
                    },
                    .float => {
                        if (denominator == 0.0)
                            return Error.ZeroDenominator;

                        var dvalue = try @import("float/asRational.zig").asRational(denominator);
                        dvalue[0].num.limbs = &dvalue[1][0];
                        dvalue[0].den.limbs = &dvalue[1][1];

                        try div_(allocator, self, numerator.asRational(), dvalue[0]);
                    },
                    .cfloat => return self.set(allocator, numerator, denominator.re),
                    .integer => {
                        if (denominator.size == 0)
                            return Error.ZeroDenominator;

                        try self.num.set(allocator, numerator);
                        try self.den.set(allocator, denominator);
                    },
                    .rational => {
                        if (denominator.num.size == 0)
                            return Error.ZeroDenominator;

                        try div_(allocator, self, numerator.asRational(), denominator);
                    },
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator, denominator.re),
                },
                .rational => switch (comptime types.numericType(D)) {
                    .bool => {
                        if (!denominator)
                            return Error.ZeroDenominator;

                        try self.num.set(allocator, numerator.num);
                        try self.den.set(allocator, numerator.den);
                    },
                    .int => {
                        if (denominator == 0)
                            return Error.ZeroDenominator;

                        var dvalue = @import("int/asRational.zig").asRational(denominator);
                        dvalue[0].num.limbs = &dvalue[1];

                        try div_(allocator, self, numerator, dvalue[0]);
                    },
                    .float => {
                        if (denominator == 0.0)
                            return Error.ZeroDenominator;

                        var dvalue = try @import("float/asRational.zig").asRational(denominator);
                        dvalue[0].num.limbs = &dvalue[1][0];
                        dvalue[0].den.limbs = &dvalue[1][1];

                        try div_(allocator, self, numerator, dvalue[0]);
                    },
                    .cfloat => return self.set(allocator, numerator, denominator.re),
                    .integer => {
                        if (denominator.size == 0)
                            return Error.ZeroDenominator;

                        try div_(allocator, self, numerator, denominator.asRational());
                    },
                    .rational => {
                        if (denominator.num.size == 0)
                            return Error.ZeroDenominator;

                        try div_(allocator, self, numerator, denominator);
                    },
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator, denominator.re),
                },
                .real => @compileError("Real type not supported yet"),
                .complex => switch (comptime types.numericType(D)) {
                    .bool => return self.set(allocator, numerator.re, denominator),
                    .int => return self.set(allocator, numerator.re, denominator),
                    .float => return self.set(allocator, numerator.re, denominator),
                    .cfloat => return self.set(allocator, numerator.re, denominator.re),
                    .integer => return self.set(allocator, numerator.re, denominator),
                    .rational => return self.set(allocator, numerator.re, denominator),
                    .real => @compileError("Real type not supported yet"),
                    .complex => return self.set(allocator, numerator.re, denominator.re),
                },
            }
        } else if (comptime N == []const u8 or N == []u8) {
            @compileError("String type not supported yet");
        } else {
            @compileError("Value must be a numeric type or a string");
        }
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

    /// Converts the `Rational` to an int type `Int`, performing truncation if
    /// necessary.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to convert.
    ///
    /// `Int` (`type`):
    /// The int type to convert to.
    ///
    /// Returns
    /// -------
    /// `Int`:
    /// The converted int value.
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

    /// Converts the `Rational` to a float type `Float`.
    ///
    /// Parameters
    /// ----------
    /// `self` (`*const Rational`):
    /// A pointer to the `Rational` to convert.
    ///
    /// `Float` (`type`):
    /// The float type to convert to.
    ///
    /// Returns
    /// -------
    /// `Float`:
    /// The converted float value.
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
            .flags = .{ .owns_data = false, .writable = true },
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
            .flags = .{ .owns_data = true, .writable = true },
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
pub const cmp = @import("rational/cmp.zig").cmp;
pub const eq = @import("rational/eq.zig").eq;
pub const ne = @import("rational/ne.zig").ne;
pub const lt = @import("rational/lt.zig").lt;
pub const le = @import("rational/le.zig").le;
pub const gt = @import("rational/gt.zig").gt;
pub const ge = @import("rational/ge.zig").ge;

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
