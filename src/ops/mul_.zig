const std = @import("std");

const types = @import("../types.zig");
const ops = @import("../ops.zig");
const int = @import("../int.zig");
const float = @import("../float.zig");
const cfloat = @import("../cfloat.zig");
const integer = @import("../integer.zig");
const rational = @import("../rational.zig");
const real = @import("../real.zig");
const complex = @import("../complex.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");
const array = @import("../array.zig");

/// Performs in-place multiplication between two operands of compatible types.
///
/// The `mul_` routine computes the product `x * y` and stores the result
/// directly into `o`, automatically validating the provided context. The
/// operation is performed in the coerced precision of the operands, and the
/// result is then cast to the output type. It supports both fixed-precision and
/// arbitrary-precision arithmetic, as well as structured data domains. The
/// supported type combinations are:
/// - **Numeric = Numeric * Numeric**: scalar multiplication.
/// - **Numeric = Vector * Vector**: vector dot product.
/// - **Vector = Numeric * Vector** and **Vector = Vector * Numeric**:
///   element-wise multiplication.
/// - **Vector = Vector * Matrix** and **Vector = Matrix * Vector**:
///   matrix-vector multiplication, with the vector treated as a row or column
///   vector, respectively.
/// - **Matrix = Numeric * Matrix** and **Matrix = Matrix * Numeric**:
///   element-wise multiplication.
/// - **Matrix = Matrix * Matrix**: matrix-matrix multiplication.
/// - **Array = Numeric * Array**, **Array = Array * Numeric**, and **Array =
///   Array * Array**: broadcasted element-wise multiplication.
/// - **Expression = Any * Expression**, **Expression = Expression * Any**, and
///   **Expression = Expression * Expression**: symbolic multiplication.
///
/// Signature
/// ---------
/// ```zig
/// fn mul_(o: *O, x: X, y: Y, ctx: anytype) !void
/// ```
///
/// Parameters
/// ----------
/// `o` (`anytype`):
/// The output pointer where the result will be stored. For arbitrary-precision
/// or structured types, `o` must point to a properly initialized value.
///
/// `x` (`anytype`):
/// The left operand.
///
/// `y` (`anytype`):
/// The right operand.
///
/// `ctx` (`anytype`):
/// A context struct providing necessary resources and configuration for the
/// operation. The required fields depend on the output and operand types. If
/// the context is missing required fields or contains unnecessary or wrongly
/// typed fields, the compiler will emit a detailed error message describing the
/// expected structure.
///
/// Returns
/// -------
/// `void`
///
/// Errors
/// ------
/// ``:
///
/// Notes
/// -----
/// When the output and either of the input types are the same, aliasing is
/// allowed.
pub inline fn mul_(
    o: anytype,
    x: anytype,
    y: anytype,
    ctx: anytype,
) !void {
    comptime var O: type = @TypeOf(o);
    const X: type = @TypeOf(x);
    const Y: type = @TypeOf(y);

    comptime if (!types.isPointer(O) or types.isConstPointer(O))
        @compileError("zml.mul_ requires the output to be a mutable pointer, got " ++ @typeName(O));

    O = types.Child(O);

    comptime if (!types.isExpression(O) and !types.isExpression(X) and !types.isExpression(Y) and
        !types.isArray(O) and !types.isArray(X) and !types.isArray(Y) and
        !types.isMatrix(O) and !types.isMatrix(X) and !types.isMatrix(Y) and
        !types.isVector(O) and !types.isVector(X) and !types.isVector(Y) and
        !types.isNumeric(O) and !types.isNumeric(X) and !types.isNumeric(Y))
        @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

    const C: type = types.Coerce(O, types.Coerce(X, Y));

    switch (comptime types.domainType(O)) {
        .expression => @compileError("zml.mul_ not implemented yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        .array => switch (comptime types.domainType(X)) {
            .array, .numeric => switch (comptime types.domainType(Y)) {
                .array, .numeric => { // array = array * array, array = numeric * array, array = array * numeric
                    comptime if (types.isNumeric(X) and types.isNumeric(Y))
                        @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types");

                    comptime switch (types.numericType(types.Numeric(O))) {
                        .bool, .int, .float, .cfloat => switch (types.numericType(types.Numeric(C))) {
                            .bool => @compileError("zml.add_ not defined for " ++ @typeName(X) ++ " and " ++ @typeName(Y)),
                            .int => {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .mode = .{ .type = int.Mode, .required = false, .default = .default },
                                    },
                                );
                            },
                            .float, .cfloat => {
                                types.validateContext(@TypeOf(ctx), .{});
                            },
                            .integer, .rational, .real, .complex => {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );
                            },
                        },
                        .integer, .rational, .real, .complex => switch (types.numericType(types.Numeric(C))) {
                            .bool, .float, .cfloat => {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );
                            },
                            .int => {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .mode = .{ .type = int.Mode, .required = false, .default = .default },
                                    },
                                );
                            },
                            .integer, .rational, .real, .complex => {
                                types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .element_allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                    },
                                );
                            },
                        },
                    };

                    return array.mul_(
                        o,
                        x,
                        y,
                        ctx,
                    );
                },
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .matrix => switch (comptime types.domainType(X)) {
            .matrix => switch (comptime types.domainType(Y)) {
                .matrix => @compileError("matrix.mul_ not implemented yet"),
                .numeric => @compileError("matrix.mul_ not implemented yet"),
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            .numeric => switch (comptime types.domainType(Y)) {
                .matrix => @compileError("matrix.mul_ not implemented yet"),
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .vector => switch (comptime types.domainType(X)) {
            .vector => switch (comptime types.domainType(Y)) {
                .matrix => @compileError("vector.mul_ not implemented yet"),
                .numeric => @compileError("vector.mul_ not implemented yet"),
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            .numeric => switch (comptime types.domainType(Y)) {
                .matrix => @compileError("vector.mul_ not implemented yet"),
                .vector => @compileError("vector.mul_ not implemented yet"),
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
        .numeric => switch (comptime types.domainType(X)) {
            .vector => switch (comptime types.domainType(Y)) {
                .vector => @compileError("vector.mul_ not implemented yet"),
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            .numeric => switch (comptime types.domainType(Y)) {
                .numeric => { // numeric = numeric * numeric
                    switch (comptime types.numericType(O)) {
                        .bool, .int, .float, .cfloat => switch (comptime types.numericType(C)) {
                            .bool => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .int => {
                                const spec =
                                    .{
                                        .mode = .{ .type = int.Mode, .required = false, .default = .default },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                ops.set(
                                    o,
                                    int.mul(
                                        x,
                                        y,
                                        types.getFieldOrDefault(ctx, spec, "mode"),
                                    ),
                                    .{},
                                ) catch unreachable;
                            },
                            .float => {
                                comptime types.validateContext(@TypeOf(ctx), .{});

                                ops.set(
                                    o,
                                    float.mul(x, y),
                                    .{},
                                ) catch unreachable;
                            },
                            .cfloat => {
                                comptime types.validateContext(@TypeOf(ctx), .{});

                                ops.set(
                                    o,
                                    cfloat.mul(x, y),
                                    .{},
                                ) catch unreachable;
                            },
                            .integer => {
                                const spec =
                                    .{
                                        .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try integer.mul_(
                                        ctx.buffer_allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    ops.set(
                                        o,
                                        buffer.*,
                                        .{},
                                    ) catch unreachable;
                                } else {
                                    var result: integer.Integer = try integer.mul(
                                        ctx.buffer_allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(ctx.buffer_allocator);

                                    ops.set(
                                        o,
                                        result,
                                        .{},
                                    ) catch unreachable;
                                }
                            },
                            .rational => {
                                const spec =
                                    .{
                                        .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try rational.mul_(
                                        ctx.buffer_allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    ops.set(
                                        o,
                                        buffer.*,
                                        .{},
                                    ) catch unreachable;
                                } else {
                                    var result: rational.Rational = try rational.mul(
                                        ctx.buffer_allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(ctx.buffer_allocator);

                                    ops.set(
                                        o,
                                        result,
                                        .{},
                                    ) catch unreachable;
                                }
                            },
                            .real => @compileError("zml.mul_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .complex => {
                                const spec =
                                    .{
                                        .buffer_allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer = .{ .type = ?*types.Coerce(X, Y), .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try complex.mul_(
                                        ctx.buffer_allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    ops.set(
                                        o,
                                        buffer.*,
                                        .{},
                                    ) catch unreachable;
                                } else {
                                    var result: types.Coerce(X, Y) = try complex.mul(
                                        ctx.buffer_allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(ctx.buffer_allocator);

                                    ops.set(
                                        o,
                                        result,
                                        .{},
                                    ) catch unreachable;
                                }
                            },
                        },
                        .integer => switch (comptime types.numericType(C)) {
                            .bool => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .int => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .mode = .{ .type = int.Mode, .required = false, .default = .default },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                try ops.set(
                                    o,
                                    int.mul(
                                        x,
                                        y,
                                        types.getFieldOrDefault(ctx, spec, "mode"),
                                    ),
                                    types.stripStruct(ctx, &.{"mode"}),
                                );
                            },
                            .float => {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    float.mul(x, y),
                                    ctx,
                                );
                            },
                            .cfloat => {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    cfloat.mul(x, y),
                                    ctx,
                                );
                            },
                            .integer => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (@import("../integer/check_aliasing.zig").check_aliasing(o, x) or
                                    @import("../integer/check_aliasing.zig").check_aliasing(o, y))
                                {
                                    // Aliasing -> use buffer if provided
                                    if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                        try integer.mul_(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            buffer,
                                            x,
                                            y,
                                        );

                                        try ops.set(
                                            o,
                                            buffer.*,
                                            .{ .allocator = ctx.allocator },
                                        );
                                    } else {
                                        var tx = try @import("../integer/check_aliasing_alloc.zig").check_aliasing_alloc(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            o,
                                            x,
                                        );
                                        defer ops.deinit(
                                            &tx,
                                            if (comptime types.isArbitraryPrecision(X))
                                                .{
                                                    .allocator = if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                        buffer_allocator
                                                    else
                                                        ctx.allocator,
                                                }
                                            else
                                                .{},
                                        );
                                        var ty = try @import("../integer/check_aliasing_alloc.zig").check_aliasing_alloc(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            o,
                                            y,
                                        );
                                        defer ops.deinit(
                                            &ty,
                                            if (comptime types.isArbitraryPrecision(Y))
                                                .{
                                                    .allocator = if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                        buffer_allocator
                                                    else
                                                        ctx.allocator,
                                                }
                                            else
                                                .{},
                                        );

                                        try integer.mul_(
                                            ctx.allocator,
                                            o,
                                            tx,
                                            ty,
                                        );
                                    }
                                } else {
                                    // No aliasing -> no need for buffer
                                    try integer.mul_(
                                        ctx.allocator,
                                        o,
                                        x,
                                        y,
                                    );
                                }
                            },
                            .rational => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try rational.mul_(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    try ops.set(
                                        o,
                                        buffer.*,
                                        .{ .allocator = ctx.allocator },
                                    );
                                } else {
                                    var result: rational.Rational = try rational.mul(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                    );

                                    try ops.set(
                                        o,
                                        result,
                                        .{ .allocator = ctx.allocator },
                                    );
                                }
                            },
                            .real => @compileError("zml.mul_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .complex => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*types.Coerce(X, Y), .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try complex.mul_(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    try ops.set(
                                        o,
                                        buffer.*,
                                        .{ .allocator = ctx.allocator },
                                    );
                                } else {
                                    var result: types.Coerce(X, Y) = try complex.mul(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                    );

                                    try ops.set(
                                        o,
                                        result,
                                        .{ .allocator = ctx.allocator },
                                    );
                                }
                            },
                        },
                        .rational => switch (comptime types.numericType(C)) {
                            .bool => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .int => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .mode = .{ .type = int.Mode, .required = false, .default = .default },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                try ops.set(
                                    o,
                                    int.mul(
                                        x,
                                        y,
                                        types.getFieldOrDefault(ctx, spec, "mode"),
                                    ),
                                    types.stripStruct(ctx, &.{"mode"}),
                                );
                            },
                            .float => {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    float.mul(x, y),
                                    ctx,
                                );
                            },
                            .cfloat => {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    cfloat.mul(x, y),
                                    ctx,
                                );
                            },
                            .integer => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try integer.mul_(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    try ops.set(
                                        o,
                                        buffer.*,
                                        .{ .allocator = ctx.allocator },
                                    );
                                } else {
                                    var result: integer.Integer = try integer.mul(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                    );

                                    try ops.set(
                                        o,
                                        result,
                                        .{ .allocator = ctx.allocator },
                                    );
                                }
                            },
                            .rational => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (@import("../rational/check_aliasing.zig").check_aliasing(o, x) or
                                    @import("../rational/check_aliasing.zig").check_aliasing(o, y))
                                {
                                    // Aliasing -> use buffer if provided
                                    if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                        try rational.mul_(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            buffer,
                                            x,
                                            y,
                                        );

                                        try ops.set(
                                            o,
                                            buffer.*,
                                            .{ .allocator = ctx.allocator },
                                        );
                                    } else {
                                        var tx = try @import("../rational/check_aliasing_alloc.zig").check_aliasing_alloc(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            o,
                                            x,
                                        );
                                        defer ops.deinit(
                                            &tx,
                                            if (comptime types.isArbitraryPrecision(X))
                                                .{
                                                    .allocator = if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                        buffer_allocator
                                                    else
                                                        ctx.allocator,
                                                }
                                            else
                                                .{},
                                        );
                                        var ty = try @import("../rational/check_aliasing_alloc.zig").check_aliasing_alloc(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            o,
                                            y,
                                        );
                                        defer ops.deinit(
                                            &ty,
                                            if (comptime types.isArbitraryPrecision(Y))
                                                .{
                                                    .allocator = if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                        buffer_allocator
                                                    else
                                                        ctx.allocator,
                                                }
                                            else
                                                .{},
                                        );

                                        try rational.mul_(
                                            ctx.allocator,
                                            o,
                                            tx,
                                            ty,
                                        );
                                    }
                                } else {
                                    // No aliasing -> no need for buffer
                                    try rational.mul_(
                                        ctx.allocator,
                                        o,
                                        x,
                                        y,
                                    );
                                }
                            },
                            .real => @compileError("zml.mul_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .complex => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*types.Coerce(X, Y), .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try complex.mul_(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    try ops.set(
                                        o,
                                        buffer.*,
                                        .{ .allocator = ctx.allocator },
                                    );
                                } else {
                                    var result: types.Coerce(X, Y) = try complex.mul(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                    );

                                    try ops.set(
                                        o,
                                        result,
                                        .{ .allocator = ctx.allocator },
                                    );
                                }
                            },
                        },
                        .real => @compileError("zml.mul_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                        .complex => switch (comptime types.numericType(C)) {
                            .bool => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .int => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .mode = .{ .type = int.Mode, .required = false, .default = .default },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                try ops.set(
                                    o,
                                    int.mul(
                                        x,
                                        y,
                                        types.getFieldOrDefault(ctx, "mode", int.Mode, .default),
                                    ),
                                    types.stripStruct(ctx, &.{"mode"}),
                                );
                            },
                            .float => {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    float.mul(x, y),
                                    ctx,
                                );
                            },
                            .cfloat => {
                                comptime types.validateContext(
                                    @TypeOf(ctx),
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                    },
                                );

                                try ops.set(
                                    o,
                                    cfloat.mul(x, y),
                                    ctx,
                                );
                            },
                            .integer => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*integer.Integer, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try integer.mul_(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    try ops.set(
                                        o,
                                        buffer.*,
                                        .{ .allocator = ctx.allocator },
                                    );
                                } else {
                                    var result: integer.Integer = try integer.mul(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                    );

                                    try ops.set(
                                        o,
                                        result,
                                        .{ .allocator = ctx.allocator },
                                    );
                                }
                            },
                            .rational => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*rational.Rational, .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                    try rational.mul_(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        buffer,
                                        x,
                                        y,
                                    );

                                    try ops.set(
                                        o,
                                        buffer.*,
                                        .{ .allocator = ctx.allocator },
                                    );
                                } else {
                                    var result: rational.Rational = try rational.mul(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                        x,
                                        y,
                                    );
                                    defer result.deinit(
                                        if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                            buffer_allocator
                                        else
                                            ctx.allocator,
                                    );

                                    try ops.set(
                                        o,
                                        result,
                                        .{ .allocator = ctx.allocator },
                                    );
                                }
                            },
                            .real => @compileError("zml.mul_ not implemeneted yet for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
                            .complex => {
                                const spec =
                                    .{
                                        .allocator = .{ .type = std.mem.Allocator, .required = true },
                                        .buffer_allocator = .{ .type = ?std.mem.Allocator, .required = false, .default = null },
                                        .buffer = .{ .type = ?*types.Coerce(X, Y), .required = false, .default = null },
                                    };

                                comptime types.validateContext(@TypeOf(ctx), spec);

                                if (@import("../complex/check_aliasing.zig").check_aliasing(o, x) or
                                    @import("../complex/check_aliasing.zig").check_aliasing(o, y))
                                {
                                    // Aliasing -> use buffer if provided
                                    if (types.getFieldOrDefault(ctx, spec, "buffer")) |buffer| {
                                        try complex.mul_(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            buffer,
                                            x,
                                            y,
                                        );

                                        try ops.set(
                                            o,
                                            buffer.*,
                                            .{ .allocator = ctx.allocator },
                                        );
                                    } else {
                                        var tx = try @import("../complex/check_aliasing_alloc.zig").check_aliasing_alloc(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            o,
                                            x,
                                        );
                                        defer ops.deinit(
                                            &tx,
                                            if (comptime X == integer.Integer or X == rational.Rational or X == real.Real or
                                                X == complex.Complex(rational.Rational) or X == complex.Complex(real.Real))
                                                .{
                                                    .allocator = if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                        buffer_allocator
                                                    else
                                                        ctx.allocator,
                                                }
                                            else
                                                .{},
                                        );
                                        var ty = try @import("../complex/check_aliasing_alloc.zig").check_aliasing_alloc(
                                            if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                buffer_allocator
                                            else
                                                ctx.allocator,
                                            o,
                                            y,
                                        );
                                        defer ops.deinit(
                                            &ty,
                                            if (comptime Y == integer.Integer or Y == rational.Rational or Y == real.Real or
                                                Y == complex.Complex(rational.Rational) or Y == complex.Complex(real.Real))
                                                .{
                                                    .allocator = if (types.getFieldOrDefault(ctx, spec, "buffer_allocator")) |buffer_allocator|
                                                        buffer_allocator
                                                    else
                                                        ctx.allocator,
                                                }
                                            else
                                                .{},
                                        );

                                        try complex.mul_(
                                            ctx.allocator,
                                            o,
                                            x,
                                            y,
                                        );
                                    }
                                } else {
                                    // No aliasing -> no need for buffer
                                    try complex.mul_(
                                        ctx.allocator,
                                        o,
                                        x,
                                        y,
                                    );
                                }
                            },
                        },
                    }
                },
                else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
            },
            else => @compileError("zml.mul_ not defined for " ++ @typeName(O) ++ ", " ++ @typeName(X) ++ " and " ++ @typeName(Y) ++ " types"),
        },
    }
}
