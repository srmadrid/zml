const std = @import("std");

const types = @import("../types.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const complex = @import("../complex.zig");

pub fn check_aliasing(o: anytype, x: anytype) bool {
    switch (comptime types.numericType(types.Scalar(types.Child(@TypeOf(o))))) {
        .rational => switch (comptime types.numericType(@TypeOf(x))) {
            .bool, .int, .float, .cfloat => return false,
            .integer => return o.re.num.limbs == x.limbs or o.re.den.limbs == x.limbs or
                o.im.num.limbs == x.limbs or o.im.den.limbs == x.limbs,
            .rational => return o.re.num.limbs == x.num.limbs or o.re.num.limbs == x.den.limbs or
                o.re.den.limbs == x.num.limbs or o.re.den.limbs == x.den.limbs or
                o.im.num.limbs == x.num.limbs or o.im.num.limbs == x.den.limbs or
                o.im.den.limbs == x.num.limbs or o.im.den.limbs == x.den.limbs,
            .real => @compileError("zml.complex.check_aliasing not implemented yet for complex and real types"),
            .complex => switch (comptime types.numericType(types.Scalar(types.Child(@TypeOf(x))))) {
                .rational => return o.re.num.limbs == x.re.num.limbs or o.re.num.limbs == x.re.den.limbs or
                    o.re.num.limbs == x.im.num.limbs or o.re.num.limbs == x.im.den.limbs or
                    o.re.den.limbs == x.re.num.limbs or o.re.den.limbs == x.re.den.limbs or
                    o.re.den.limbs == x.im.num.limbs or o.re.den.limbs == x.im.den.limbs or
                    o.im.num.limbs == x.re.num.limbs or o.im.num.limbs == x.re.den.limbs or
                    o.im.num.limbs == x.im.num.limbs or o.im.num.limbs == x.im.den.limbs or
                    o.im.den.limbs == x.re.num.limbs or o.im.den.limbs == x.re.den.limbs or
                    o.im.den.limbs == x.im.num.limbs or o.im.den.limbs == x.im.den.limbs,
                .real => @compileError("zml.complex.check_aliasing not implemented yet for complex and real types"),
                else => unreachable,
            },
            else => unreachable,
        },
        .real => @compileError("zml.complex.check_aliasing not implemented yet for real output types"),
        else => unreachable,
    }
}
