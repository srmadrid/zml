const types = @import("../../types.zig");

const dense = @import("../dense.zig");
const strided = @import("../strided.zig");
// const sparse = @import("../sparse.zig");

///
pub fn apply2_(
    o: anytype,
    x: anytype,
    y: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("apply2_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O))
        @compileError("apply2_: o must be an array, got " ++ @typeName(O));

    comptime if (!types.isArray(X) and !types.isArray(Y))
        @compileError("apply2_: at least one of x or y must be an array, got " ++ @typeName(X) ++ " and " ++ @typeName(Y));

    comptime if (@typeInfo(@TypeOf(op_)) != .@"fn" or (@typeInfo(@TypeOf(op_)).@"fn".params.len != 3 and @typeInfo(@TypeOf(op_)).@"fn".params.len != 4))
        @compileError("apply2_: op must be a function of three arguments, or a function of four arguments with the fourth argument being a context, got " ++ @typeName(@TypeOf(op_)));

    if (comptime !types.isArray(X)) {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(Y)) {
                .dense => return dense.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(Y)) {
                .dense => return strided.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else if (comptime !types.isArray(Y)) {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(X)) {
                .dense => return dense.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(X)) {
                .dense => return strided.apply2_(o, x, y, op_, ctx),
                .strided => return strided.apply2_(o, x, y, op_, ctx),
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(X)) {
                .dense => switch (comptime types.arrayType(Y)) {
                    .dense => return dense.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .strided => switch (comptime types.arrayType(Y)) {
                    .dense => return strided.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(X)) {
                .dense => switch (comptime types.arrayType(Y)) {
                    .dense => return strided.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .strided => switch (comptime types.arrayType(Y)) {
                    .dense => return strided.apply2_(o, x, y, op_, ctx),
                    .strided => return strided.apply2_(o, x, y, op_, ctx),
                    .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                    .numeric => unreachable,
                },
                .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply2_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    }
}
