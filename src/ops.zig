const std = @import("std");

const types = @import("types.zig");
const cast = types.cast;
const scast = types.scast;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const MulCoerce = types.MulCoerce;
const EnsureArray = types.EnsureArray;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;

const int = @import("int.zig");
const float = @import("float.zig");
const cfloat = @import("cfloat.zig");
const integer = @import("integer.zig");
const rational = @import("rational.zig");
const real = @import("real.zig");
const complex = @import("complex.zig");
const vector = @import("vector.zig");
const matrix = @import("matrix.zig");
const array = @import("array.zig");

// Utility operations
pub const copysign = @import("ops/copysign.zig").copysign;
pub const re = @import("ops/re.zig").re;
pub const re_ = @import("ops/re_.zig").re_;
pub const im = @import("ops/im.zig").im;
pub const im_ = @import("ops/im_.zig").im_;
pub const conj = @import("ops/conj.zig").conj;
pub const conj_ = @import("ops/conj_.zig").conj_;

// Arithmetic operations
pub const add = @import("ops/add.zig").add;
pub const add_ = @import("ops/add_.zig").add_;
pub const sub = @import("ops/sub.zig").sub;
pub const sub_ = @import("ops/sub_.zig").sub_;
pub const mul = @import("ops/mul.zig").mul;
pub const mul_ = @import("ops/mul_.zig").mul_;
pub const div = @import("ops/div.zig").div;
pub const div_ = @import("ops/div_.zig").div_;

// Comparison operations
pub const eq = @import("ops/eq.zig").eq;
pub const eq_ = @import("ops/eq_.zig").eq_;
pub const ne = @import("ops/ne.zig").ne;
pub const ne_ = @import("ops/ne_.zig").ne_;
pub const lt = @import("ops/lt.zig").lt;
pub const lt_ = @import("ops/lt_.zig").lt_;
pub const le = @import("ops/le.zig").le;
pub const le_ = @import("ops/le_.zig").le_;
pub const gt = @import("ops/gt.zig").gt;
pub const gt_ = @import("ops/gt_.zig").gt_;
pub const ge = @import("ops/ge.zig").ge;
pub const ge_ = @import("ops/ge_.zig").ge_;
pub const max = @import("ops/max.zig").max;
pub const max_ = @import("ops/max_.zig").max_;
pub const min = @import("ops/min.zig").min;
pub const min_ = @import("ops/min_.zig").min_;

// Basic operations
pub const abs = @import("ops/abs.zig").abs;
pub const abs_ = @import("ops/abs_.zig").abs_;
pub const abs1 = @import("ops/abs1.zig").abs1;
pub const abs1_ = @import("ops/abs1_.zig").abs1_;
pub const abs2 = @import("ops/abs2.zig").abs2;
pub const abs2_ = @import("ops/abs2_.zig").abs2_;
pub const neg = @import("ops/neg.zig").neg;
pub const neg_ = @import("ops/neg_.zig").neg_;

// Exponential functions
pub const exp = @import("ops/exp.zig").exp;
pub const exp_ = @import("ops/exp_.zig").exp_;
pub const exp10 = @import("ops/exp10.zig").exp10;
pub const exp10_ = @import("ops/exp10_.zig").exp10_;
pub const exp2 = @import("ops/exp2.zig").exp2;
pub const exp2_ = @import("ops/exp2_.zig").exp2_;
pub const log = @import("ops/log.zig").log;
pub const log_ = @import("ops/log_.zig").log_;
pub const log10 = @import("ops/log10.zig").log10;
pub const log10_ = @import("ops/log10_.zig").log10_;
pub const log2 = @import("ops/log2.zig").log2;
pub const log2_ = @import("ops/log2_.zig").log2_;

// Power functions
pub const pow = @import("ops/pow.zig").pow;
pub const pow_ = @import("ops/pow_.zig").pow_;
pub const sqrt = @import("ops/sqrt.zig").sqrt;
pub const sqrt_ = @import("ops/sqrt_.zig").sqrt_;
pub const cbrt = @import("ops/cbrt.zig").cbrt;
pub const cbrt_ = @import("ops/cbrt_.zig").cbrt_;
pub const hypot = @import("ops/hypot.zig").hypot;
pub const hypot_ = @import("ops/hypot_.zig").hypot_;

// Trigonometric functions
pub const sin = @import("ops/sin.zig").sin;
pub const sin_ = @import("ops/sin_.zig").sin_;
pub const cos = @import("ops/cos.zig").cos;
pub const cos_ = @import("ops/cos_.zig").cos_;
pub const tan = @import("ops/tan.zig").tan;
pub const tan_ = @import("ops/tan_.zig").tan_;
pub const asin = @import("ops/asin.zig").asin;
pub const asin_ = @import("ops/asin_.zig").asin_;
pub const acos = @import("ops/acos.zig").acos;
pub const acos_ = @import("ops/acos_.zig").acos_;
pub const atan = @import("ops/atan.zig").atan;
pub const atan_ = @import("ops/atan_.zig").atan_;
pub const atan2 = @import("ops/atan2.zig").atan2;
pub const atan2_ = @import("ops/atan2_.zig").atan2_;

// Hyperbolic functions
pub const sinh = @import("ops/sinh.zig").sinh;
pub const sinh_ = @import("ops/sinh_.zig").sinh_;
pub const cosh = @import("ops/cosh.zig").cosh;
pub const cosh_ = @import("ops/cosh_.zig").cosh_;
pub const tanh = @import("ops/tanh.zig").tanh;
pub const tanh_ = @import("ops/tanh_.zig").tanh_;
pub const asinh = @import("ops/asinh.zig").asinh;
pub const asinh_ = @import("ops/asinh_.zig").asinh_;
pub const acosh = @import("ops/acosh.zig").acosh;
pub const acosh_ = @import("ops/acosh_.zig").acosh_;
pub const atanh = @import("ops/atanh.zig").atanh;
pub const atanh_ = @import("ops/atanh_.zig").atanh_;

// Error and gamma functions
pub const erf = @import("ops/erf.zig").erf;
pub const erf_ = @import("ops/erf_.zig").erf_;
pub const erfc = @import("ops/erfc.zig").erfc;
pub const erfc_ = @import("ops/erfc_.zig").erfc_;
pub const gamma = @import("ops/gamma.zig").gamma;
pub const gamma_ = @import("ops/gamma_.zig").gamma_;
pub const lgamma = @import("ops/lgamma.zig").lgamma;
pub const lgamma_ = @import("ops/lgamma_.zig").lgamma_;

// Nearest integer operations
pub inline fn ceil(
    x: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Numeric(@TypeOf(x))) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X)) {
        comptime if (types.isArbitraryPrecision(Numeric(X))) {
            if (types.numericType(Numeric(X)) == .int) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{
                        .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                    },
                );
            }
        } else {
            types.validateContext(
                @TypeOf(ctx),
                .{
                    .array_allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );
        };

        return array.ceil(
            ctx.array_allocator,
            x,
            types.stripStruct(ctx, &.{"array_allocator"}),
        );
    } else if (comptime types.isMatrix(X)) {
        @compileError("zml.ceil not defined for matrices, convert to array first");
    } else {
        switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.ceil not defined for " ++ @typeName(X)),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return x;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                return float.ceil(x);
            },
            .cfloat => @compileError("zml.ceil not defined for " ++ @typeName(X)),
            else => @compileError("zml.ceil not implemented for " ++ @typeName(X) ++ " yet"),
        }
    }
}

pub inline fn ceil_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.ceil_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    if (comptime types.isArray(O)) {
        comptime if (types.isArbitraryPrecision(Numeric(O))) {
            if (types.isArbitraryPrecision(Numeric(X))) {
                if (Numeric(O) == Numeric(X)) {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                } else {
                    types.validateContext(
                        @TypeOf(ctx),
                        .{
                            .allocator = .{ .type = std.mem.Allocator, .required = true },
                            .internal_allocator = .{ .type = std.mem.Allocator, .required = true },
                        },
                    );
                }
            } else {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            }
        } else {
            if (types.isArbitraryPrecision(Numeric(X))) {
                types.validateContext(
                    @TypeOf(ctx),
                    .{ .internal_allocator = .{ .type = std.mem.Allocator, .required = true } },
                );
            } else {
                types.validateContext(@TypeOf(ctx), .{});
            }
        };

        return array.ceil_(o, x, ctx);
    } else if (comptime types.isArray(X)) {
        @compileError("zml.ceil_: o must be an array if x is an array, got " ++ @typeName(O) ++ " and " ++ @typeName(X));
    } else if (comptime types.isMatrix(O) or types.isMatrix(X)) {
        @compileError("zml.ceil_ not defined for matrices, convert to array first");
    } else {
        switch (comptime types.numericType(X)) {
            .bool => @compileError("zml.ceil_ not defined for " ++ @typeName(X) ++ " input type"),
            .int => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = x;
            },
            .float => {
                comptime types.validateContext(@TypeOf(ctx), .{});

                o.* = scast(O, float.ceil(x));
            },
            .cfloat => @compileError("zml.ceil_ not defined for " ++ @typeName(X) ++ " input type"),
            else => @compileError("zml.ceil_ not implemented for " ++ @typeName(O) ++ " and " ++ @typeName(X) ++ " output and input types yet"),
        }
    }
}

pub inline fn init(
    comptime T: type,
    ctx: anytype,
) !T {
    if (comptime !types.isNumeric(T))
        @compileError("zml.init requires T to be a numeric type, got " ++ @typeName(T));

    switch (comptime types.numericType(T)) {
        .bool => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return false;
        },
        .int => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .float => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return 0;
        },
        .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return .{ .re = 0, .im = 0 };
        },
        else => @compileError("zml.init not implemented for " ++ @typeName(T) ++ " yet"),
    }
}

pub const set = @import("ops/set.zig").set;

pub inline fn copy(
    x: anytype,
    ctx: anytype,
) !@TypeOf(x) {
    const X: type = @TypeOf(x);

    if (comptime types.isArray(X))
        @compileError("zml.copy not implemented for arrayss yet.");

    if (comptime types.isMatrix(X))
        @compileError("zml.copy not implemented for matrices yet.");

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            return x;
        },
        else => @compileError("zml.copy not implemented for " ++ @typeName(X) ++ " yet"),
    }
}

pub inline fn deinit(
    x: anytype,
    ctx: anytype,
) void {
    comptime if (!types.isPointer(@TypeOf(x)) or types.isConstPointer(@TypeOf(x)))
        @compileError("zml.deinit requires x a mutable pointer, got " ++ @typeName(@TypeOf(x)));

    const X: type = Child(@TypeOf(x));

    if (comptime types.isArray(X))
        @compileError("zml.deinit not implemented for arrays yet.");

    if (comptime types.isMatrix(X))
        @compileError("zml.deinit not implemented for matrices yet.");

    switch (comptime types.numericType(X)) {
        .bool, .int, .float, .cfloat => {
            comptime types.validateContext(@TypeOf(ctx), .{});

            // No deinitialization needed for fixed precision types, this is a no-op.
        },
        .integer => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );

            x.deinit(ctx.allocator);
        },
        .rational => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );

            x.deinit(ctx.allocator);
        },
        .real => @compileError("zml.deinit not implemented for " ++ @typeName(X) ++ " yet"),
        .complex => {
            comptime types.validateContext(
                @TypeOf(ctx),
                .{
                    .allocator = .{ .type = std.mem.Allocator, .required = true },
                },
            );

            x.deinit(ctx.allocator);
        },
        .expression => @compileError("zml.deinit not implemented for " ++ @typeName(X) ++ " yet"),
    }
}
