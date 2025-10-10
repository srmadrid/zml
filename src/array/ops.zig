const std = @import("std");

const types = @import("../types.zig");
const Coerce = types.Coerce;
const Scalar = types.Scalar;
const Numeric = types.Numeric;
const Child = types.Child;
const EnsureFloat = types.EnsureFloat;
const EnsureArray = types.EnsureArray;
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;

const int = @import("../int.zig");
const ops = @import("../ops.zig");

const array = @import("../array.zig");

const dense = @import("dense.zig");
const strided = @import("strided.zig");
const sparse = @import("sparse.zig");

pub const apply1 = @import("ops/apply1.zig").apply1;
pub const apply1_ = @import("ops/apply1_.zig").apply1_;
pub const apply2 = @import("ops/apply2.zig").apply2;
pub const apply2_ = @import("ops/apply2_.zig").apply2_;

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
    allocator: std.mem.Allocator,
    x: anytype,
    ctx: anytype,
) !EnsureArray(@TypeOf(x), Numeric(@TypeOf(x))) {
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(X)) {
        types.validateContext(
            @TypeOf(ctx),
            .{ .allocator = .{ .type = std.mem.Allocator, .required = true } },
        );
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1(
        allocator,
        x,
        ops.ceil,
        ctx,
    );
}

pub inline fn ceil_(
    o: anytype,
    x: anytype,
    ctx: anytype,
) !void {
    const O: type = Numeric(Child(@TypeOf(o)));
    const X: type = Numeric(@TypeOf(x));

    comptime if (types.isArbitraryPrecision(O) or types.isArbitraryPrecision(X)) {
        @compileError("arbitrary precision types not implemented yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    return apply1_(o, x, ops.ceil_, ctx);
}
