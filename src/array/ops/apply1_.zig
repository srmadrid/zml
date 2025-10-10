const types = @import("../../types.zig");

const dense = @import("../dense.zig");
const strided = @import("../strided.zig");
// const sparse = @import("../sparse.zig");

///
pub fn apply1_(
    o: anytype,
    x: anytype,
    comptime op_: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("apply1_: o must be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isArray(O))
        @compileError("apply1_: o must be an array, got " ++ @typeName(O));

    comptime if (@typeInfo(@TypeOf(op_)) != .@"fn" or (@typeInfo(@TypeOf(op_)).@"fn".params.len != 2 and @typeInfo(@TypeOf(op_)).@"fn".params.len != 3))
        @compileError("apply1_: op_ must be a function of two arguments, or a function of three arguments with the third argument being a context, got " ++ @typeName(@TypeOf(op_)));

    if (comptime !types.isArray(X)) {
        switch (comptime types.arrayType(O)) {
            .dense => return dense.apply1_(o, x, op_, ctx),
            .strided => return strided.apply1_(o, x, op_, ctx),
            .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    } else {
        switch (comptime types.arrayType(O)) {
            .dense => switch (comptime types.arrayType(X)) {
                .dense => return dense.apply1_(o, x, op_, ctx), // dense dense apply1_
                .strided => return strided.apply1_(o, x, op_, ctx), // dense strided apply1_
                .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .strided => switch (comptime types.arrayType(X)) {
                .dense => return strided.apply1_(o, x, op_, ctx), // strided dense apply1_
                .strided => return strided.apply1_(o, x, op_, ctx), // strided strided apply1_
                .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
                .numeric => unreachable,
            },
            .sparse => @compileError("apply1_ not implemented for sparse arrays yet"),
            .numeric => unreachable,
        }
    }
}
