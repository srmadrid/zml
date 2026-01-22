const std = @import("std");

const types = @import("../types.zig");

/// Validates a context struct against a specification struct.
///
/// This function checks that the context struct matches the specification
/// struct, ensuring that all required fields are present and have the correct
/// typesn, and no unexpected fields are present. If any discrepancies are
/// found, a compile-time error is raised with a detailed message.
///
/// ## Arguments
/// * `ctx` (`comptime type`): The type of the context struct to validate. Must
///   be `void` or a struct type.
/// * `spec` (`anytype`): The specification struct that defines the expected
///   fields. Must have the following structure:
///   ```zig
///   .{
///       .field_name1 = .{ .type = type1, .required = true },
///       .field_name2 = .{ .type = type2, .required = false, .default = default_value },
///       ...
///   }
///   ```
///   where `field_name` is the name of the field, `type` is the expected type
///   of the field, `required` is a boolean indicating whether the field is
///   required or not, and `default` is the default value for optional fields.
///   Additionally, a `description` field can be included for documentation
///   purposes.
///
/// ## Returns
/// `void`
pub fn validateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);
    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (comptime ctxinfo != .@"struct" and ctxinfo != .void)
        @compileError("zml.types.validateContext: Ctx must be void or a struct, got\n\tCtx: " ++ @typeName(Ctx) ++ "\n");

    if (comptime specinfo.@"struct".fields.len == 0 and (Ctx == void or ctxinfo.@"struct".fields.len == 0))
        return; // No context needed, nothing to validate

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = field_spec.required;
        const expected_type = field_spec.type;

        if (comptime @hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);

            if (comptime !typesMatch(actual_type, expected_type))
                @compileError(formatSpecCtxMismatch(Ctx, spec));
        } else {
            if (comptime required)
                @compileError(formatSpecCtxMismatch(Ctx, spec));
        }
    }

    // Check for unexpected fields in context
    inline for (ctxinfo.@"struct".fields) |ctx_field| {
        if (!@hasField(SpecType, ctx_field.name)) {
            @compileError(formatSpecCtxMismatch(Ctx, spec));
        }
    }
}

fn formatSpecCtxMismatch(
    comptime Ctx: type,
    comptime spec: anytype,
) []const u8 {
    return std.fmt.comptimePrint(
        "Context struct must have the following structure:\n{s}\nGot:\n{s}",
        .{ formatSpec(spec), formatCtx(Ctx) },
    );
}

/// Validates that a context struct contains at least the required fields
/// defined in a specification struct.
///
/// This function differs from `validateContext` in that it only checks for the
/// presence and types of required fields, and the types of present optional
/// fields, allowing the context struct to have additional fields not specified
/// in the spec. If any required fields are missing or have incorrect types, or
/// if present optional fields have incorrect types, a compile-time error is
/// raised with a detailed message.
///
/// ## Arguments
/// * `ctx` (`comptime type`): The type of the context struct to validate. Must
///   be `void` or a struct type.
/// * `spec` (`anytype`): The specification struct that defines the expected
///   fields. Must have the following structure:
///   ```zig
///   .{
///       .field_name1 = .{ .type = type1, .required = true },
///       .field_name2 = .{ .type = type2, .required = false, .default = default_value },
///       ...
///   }
///   ```
///   where `field_name` is the name of the field, `type` is the expected type
///   of the field, `required` is a boolean indicating whether the field is
///   required or not, and `default` is the default value for optional fields.
///   Additionally, a `description` field can be included for documentation
///   purposes.
///
/// ## Returns
/// `void`
pub fn partialValidateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);
    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (comptime ctxinfo != .@"struct" and ctxinfo != .void)
        @compileError("zml.types.partialValidateContext: Ctx must be void or a struct, got\n\tCtx: " ++ @typeName(Ctx) ++ "\n");

    if (comptime specinfo.@"struct".fields.len == 0 and (Ctx == void or ctxinfo.@"struct".fields.len == 0))
        return; // No context needed, nothing to validate

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = field_spec.required;
        const expected_type = field_spec.type;

        if (comptime @hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);

            if (comptime !typesMatch(actual_type, expected_type))
                @compileError(formatPartialSpecCtxMismatch(Ctx, spec));
        } else {
            if (comptime required)
                @compileError(formatPartialSpecCtxMismatch(Ctx, spec));
        }
    }
}

fn formatPartialSpecCtxMismatch(
    comptime Ctx: type,
    comptime spec: anytype,
) []const u8 {
    return std.fmt.comptimePrint(
        "Context struct must at least have the following structure (additional fields are allowed):\n{s}\nGot:\n{s}",
        .{ formatSpec(spec), formatCtx(Ctx) },
    );
}

pub fn ctxHasField(
    comptime T: type,
    comptime field_name: []const u8,
    comptime FieldType: type,
) bool {
    if (@hasField(T, field_name)) {
        if (comptime typesMatch(@FieldType(T, field_name), FieldType))
            return true;
    }

    return false;
}

/// Retrieves a field from a context struct, or returns the default value from
/// a specification struct if the field is not present.
///
/// ## Arguments
/// * `ctx` (`anytype`): The context struct instance from which to retrieve the
///   field.
/// * `spec` (`anytype`): The specification struct that defines the expected
///   fields. Must have the following structure:
///   ```zig
///   .{
///       .field_name1 = .{ .type = type1, .required = true },
///       .field_name2 = .{ .type = type2, .required = false, .default = default_value },
///       ...
///   }
///   ```
///   where `field_name` is the name of the field, `type` is the expected type
///   of the field, `required` is a boolean indicating whether the field is
///   required or not, and `default` is the default value for optional fields.
///   Additionally, a `description` field can be included for documentation
///   purposes.
/// * field_name (`[]const u8`): The name of the field to retrieve.
///
/// ## Returns
/// `@field(spec, field_name).type`: The value of the field from the context
/// struct if it exists and has the correct type, or the default value from the
/// specification struct otherwise.
pub fn getFieldOrDefault(ctx: anytype, comptime spec: anytype, comptime field_name: []const u8) @field(spec, field_name).type {
    const Ctx = @TypeOf(ctx);
    const FieldType = @field(spec, field_name).type;

    if (@hasField(Ctx, field_name)) {
        const actual_type = @FieldType(Ctx, field_name);

        if (comptime !typesMatch(actual_type, FieldType))
            @compileError("zml.types.getFieldOrDefault: field '" ++ field_name ++ "' has type " ++ @typeName(@FieldType(Ctx, field_name)) ++ ", expected " ++ @typeName(FieldType));

        return @field(ctx, field_name);
    }

    return @field(spec, field_name).default;
}

fn typesMatch(
    comptime actual_type: type,
    comptime expected_type: type,
) bool {
    const actual_info = @typeInfo(actual_type);
    const expected_info = @typeInfo(expected_type);

    if (expected_info == .@"fn") {
        if (comptime actual_info != .@"fn")
            return false;

        // Compare function signatures
        if (comptime actual_info.@"fn".params.len != expected_info.@"fn".params.len)
            return false;

        inline for (actual_info.@"fn".params, expected_info.@"fn".params) |actual_param, expected_param| {
            if (comptime expected_param.is_generic)
                continue;

            if (comptime actual_param.type.? != expected_param.type.?)
                return false;
        }

        // Compare return types
        if (comptime actual_info.@"fn".return_type == null) // Cannot be resolved, assume match
            return true;

        const actual_ret_info = @typeInfo(actual_info.@"fn".return_type.?);
        const expected_ret_info = @typeInfo(expected_info.@"fn".return_type.?);
        if (comptime expected_ret_info == .error_union) {
            if (comptime actual_ret_info != .error_union)
                return false;

            return actual_ret_info.error_union.payload == expected_ret_info.error_union.payload;
        }

        if (comptime actual_ret_info == .error_union)
            return false;

        return actual_info.@"fn".return_type.? == expected_info.@"fn".return_type.?;
    }

    if (comptime actual_type == @TypeOf(.enum_literal)) // Special case for enum literals
        return expected_info == .@"enum" or (expected_info == .optional and @typeInfo(expected_info.optional.child) == .@"enum");

    if (comptime expected_info == .optional) // Handle optional types
        return if (comptime actual_info == .optional)
            actual_info.optional.child == expected_info.optional.child or actual_type == @TypeOf(null)
        else
            actual_type == expected_info.optional.child or actual_type == @TypeOf(null);

    return actual_type == expected_type;
}

fn formatSpec(
    comptime spec: anytype,
) []const u8 {
    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    comptime var spec_str: []const u8 = "";
    if (specinfo.@"struct".fields.len == 0) {
        spec_str = ".{}\n";
    } else {
        comptime var max_name_len: usize = 0;
        comptime var max_type_len: usize = 0;
        inline for (specinfo.@"struct".fields) |field| {
            const field_type = @typeName(@field(@field(spec, field.name), "type"));
            if (field.name.len > max_name_len) max_name_len = field.name.len;
            if (field_type.len > max_type_len) max_type_len = field_type.len;
        }

        spec_str = spec_str ++ ".{\n";
        inline for (specinfo.@"struct".fields, 0..) |field, i| {
            // Check for description field
            if (@hasField(@TypeOf(@field(spec, field.name)), "description")) {
                if (i != 0)
                    spec_str = spec_str ++ "\n";

                spec_str = spec_str ++ std.fmt.comptimePrint(
                    "    // {s}\n",
                    .{@field(@field(spec, field.name), "description")},
                );
            }

            const field_type = @typeName(@field(@field(spec, field.name), "type"));
            const required = @hasField(@TypeOf(@field(spec, field.name)), "required") and @field(@field(spec, field.name), "required");
            spec_str = spec_str ++ std.fmt.comptimePrint(
                "    .{s}: {s}",
                .{
                    field.name,
                    field_type,
                },
            );

            for (0..((max_name_len + max_type_len) - field.name.len - field_type.len)) |_| {
                spec_str = spec_str ++ " ";
            }

            spec_str = spec_str ++ std.fmt.comptimePrint(
                "    ({s})\n",
                .{
                    if (required)
                        "required"
                    else
                        std.fmt.comptimePrint(
                            "optional, default = {}",
                            .{@field(@field(spec, field.name), "default")},
                        ),
                },
            );
        }
        spec_str = spec_str ++ "}\n";
    }

    return spec_str;
}

fn formatCtx(
    comptime Ctx: type,
) []const u8 {
    const ctxinfo = @typeInfo(Ctx);

    comptime var ctx_str: []const u8 = "";
    if (ctxinfo.@"struct".fields.len == 0) {
        ctx_str = ".{}\n";
    } else {
        ctx_str = ctx_str ++ ".{\n";
        inline for (ctxinfo.@"struct".fields) |field| {
            const field_type = @typeName(@FieldType(Ctx, field.name));
            ctx_str = ctx_str ++ std.fmt.comptimePrint(
                "    .{s}: {s}\n",
                .{
                    field.name,
                    field_type,
                },
            );
        }
        ctx_str = ctx_str ++ "}\n";
    }

    return ctx_str;
}

pub fn MixStructFields(comptime S1: type, comptime S2: type) type {
    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("zml.types.MixStructFields: both S1 and S2 must be struct types, got\n\tS1: " ++
            @typeName(S1) ++ "\n\tS2: " ++ @typeName(S2) ++ "\n");

    for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("zml.type.MixStructFields: field '" ++ field.name ++ "' already exists in the second struct");
    }

    return @Type(.{
        .@"struct" = .{
            .layout = .auto,
            .fields = info1.@"struct".fields ++ info2.@"struct".fields,
            .decls = &.{},
            .is_tuple = false,
        },
    });
}

/// Merges the fields of two struct instances into a new struct instance.
/// Structs cannot have overlapping field names.
///
/// ## Arguments
/// * `s1` (`anytype`): The first struct.
/// * `s2` (`anytype`): The second struct.
///
/// ## Returns
/// `types.MixStructFields(@TypeOf(s1), @TypeOf(s2))`: A new struct instance
/// containing all fields from `s1` and `s2`.
///
/// ## Examples
/// ```zig
/// const s1 = .{ .a = 1, .b = 2 };
/// const s2 = .{ .c = 3, .d = 4 };
///
/// const mixed = mixStructFields(s1, s2);
/// // mixed is .{ .a = 1, .b = 2, .c = 3, .d = 4 }
/// ```
pub fn mixStructFields(s1: anytype, s2: anytype) MixStructFields(@TypeOf(s1), @TypeOf(s2)) {
    const S1 = @TypeOf(s1);
    const S2 = @TypeOf(s2);

    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("zml.types.mixStructFields: both s1 and s2 must be structs, got\n\ts1: " ++
            @typeName(S1) ++ "\n\ts2: " ++ @typeName(S2) ++ "\n");

    var result: MixStructFields(S1, S2) = undefined;
    inline for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("zml.types.mixStructFields: field '" ++ field.name ++ "' already exists in the second struct");

        @field(result, field.name) = @field(s1, field.name);
    }

    inline for (info2.@"struct".fields) |field| {
        @field(result, field.name) = @field(s2, field.name);
    }

    return result;
}

pub fn StripStructFields(comptime S: type, comptime fields_to_remove: []const []const u8) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("zml.types.StripStructFields: S must be a struct type, got\n\tS: " ++ @typeName(S) ++ "\n");

    // Calculate how many fields will remain
    comptime var remaining_count: usize = 0;
    inline for (info.@"struct".fields) |field| {
        comptime var should_remove: bool = false;
        inline for (fields_to_remove) |field_to_remove| {
            if (std.mem.eql(u8, field.name, field_to_remove)) {
                should_remove = true;
                break;
            }
        }

        if (!should_remove) {
            remaining_count += 1;
        }
    }

    // Create new fields array
    comptime var new_fields: [remaining_count]std.builtin.Type.StructField = undefined;
    comptime var new_field_index: usize = 0;
    inline for (info.@"struct".fields) |field| {
        var should_remove = false;
        inline for (fields_to_remove) |field_to_remove| {
            if (std.mem.eql(u8, field.name, field_to_remove)) {
                should_remove = true;
                break;
            }
        }

        if (!should_remove) {
            new_fields[new_field_index] = field;
            new_field_index += 1;
        }
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

/// Removes specified fields from a struct instance, returning a new struct
/// instance without those fields.
///
/// ## Arguments
/// * `s` (`anytype`): The struct instance from which to remove fields.
/// * fields_to_remove (`[]const []const u8`): An array of field names to
///   remove. If a specified field does not exist in `s`, it is ignored.
///
/// ## Returns
/// `types.StripStructFields(@TypeOf(s), fields_to_remove)`: A new struct
/// instance without the specified fields.
///
/// ## Examples
/// ```zig
/// const s = .{ .a = 1, .b = 2, .c = 3, .d = 4 };
///
/// const stripped = zml.types.stripStructFields(s, &.{ "b", "d" });
/// // stripped is .{ .a = 1, .c = 3 }
/// ```
pub fn stripStructFields(s: anytype, comptime fields_to_remove: []const []const u8) StripStructFields(@TypeOf(s), fields_to_remove) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);

    if (info != .@"struct")
        @compileError("zml.types.stripStructFields: s must be a struct, got\n\ts: " ++ @typeName(S) ++ "\n");

    var result: StripStructFields(S, fields_to_remove) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}

pub fn RenameStructFields(comptime S: type, comptime fields_to_rename: anytype) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("zml.types.RenameStructFields: S must be a struct type, got\n\tS: " ++ @typeName(S) ++ "\n");

    const F: type = @TypeOf(fields_to_rename);
    const finfo = @typeInfo(F);

    if (finfo != .@"struct")
        @compileError("zml.types.RenameStructFields: fields_to_rename must be a struct, got\n\tfields_to_rename: " ++ @typeName(F) ++ "\n");

    // Check that all new names do not exist in the original struct
    inline for (finfo.@"struct".fields) |field| {
        if (@hasField(S, @field(fields_to_rename, field.name)))
            @compileError("zml.types.RenameStructFields: new field name '" ++ @field(fields_to_rename, field.name) ++ "' already exists in the original struct");
    }

    // Create new fields array
    comptime var new_fields: [info.@"struct".fields.len]std.builtin.Type.StructField = undefined;
    comptime var new_field_index: usize = 0;
    inline for (info.@"struct".fields) |field| {
        comptime var renamed = false;
        inline for (finfo.@"struct".fields) |ffield| {
            if (comptime std.mem.eql(u8, field.name, ffield.name)) {
                new_fields[new_field_index] = .{
                    .name = @as([:0]const u8, @field(fields_to_rename, ffield.name)),
                    .type = field.type,
                    .default_value_ptr = field.default_value_ptr,
                    .is_comptime = field.is_comptime,
                    .alignment = field.alignment,
                };
                renamed = true;
                break;
            }
        }

        if (!renamed) {
            new_fields[new_field_index] = field;
        }
        new_field_index += 1;
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

/// Renames fields of a struct according to a mapping provided in another
/// struct.
///
/// If a field specified in `fields_to_rename` does not exist in `s`, it is
/// ignored. Fields in `s` that are not specified in `fields_to_rename` retain
/// their original names.
///
/// ## Arguments
/// * `s` (`anytype`): The struct instance whose fields are to be renamed.
/// * `fields_to_rename` (`anytype`): A struct instance that defines the mapping
///   of old field names to new field names. Each field in this struct should
///   have the name of the field to be renamed in `s`, and its value should be a
///   string representing the new namen i.e.:
///   ```zig
///   .{
///       .old_field_name1 = "new_field_name1",
///       .old_field_name2 = "new_field_name2",
///       ...
///   }
///   ```
///
/// ## Returns
/// `types.RenameStructFields(@TypeOf(s), fields_to_rename)`: A new struct
/// instance with the specified fields renamed.
///
/// ## Examples
/// ```zig
/// const s = .{ .a = 1, .b = 2 };
/// const fields_to_rename = .{
///     .a = "alpha",
///     .b = "beta",
/// };
///
/// const renamed = zml.types.renameStructFields(s, fields_to_rename);
/// // renamed is .{ .alpha = 1, .beta = 2 }
/// ```
pub fn renameStructFields(s: anytype, comptime fields_to_rename: anytype) RenameStructFields(@TypeOf(s), fields_to_rename) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("zml.types.renameStructFields: s must be a struct, got\n\ts: " ++ @typeName(S) ++ "\n");

    const F: type = @TypeOf(fields_to_rename);
    const finfo = @typeInfo(F);

    var result: RenameStructFields(S, fields_to_rename) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        comptime var renamed = false;
        inline for (finfo.@"struct".fields) |ffield| {
            if (comptime std.mem.eql(u8, field.name, @field(fields_to_rename, ffield.name))) {
                @field(result, field.name) = @field(s, ffield.name);
                renamed = true;
                break;
            }
        }

        if (!renamed) {
            @field(result, field.name) = @field(s, field.name);
        }
    }

    return result;
}

pub fn KeepStructFields(comptime S: type, comptime fields_to_keep: []const []const u8) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("zml.types.KeepStructFields: S must be a struct type, got\n\tS: " ++ @typeName(S) ++ "\n");

    // Create new fields array
    comptime var temp_new_fields: [fields_to_keep.len]std.builtin.Type.StructField = undefined;
    comptime var real_num_fields: comptime_int = 0;
    inline for (fields_to_keep) |field_to_keep| {
        comptime var field_info: std.builtin.Type.StructField = undefined;
        inline for (info.@"struct".fields) |field| {
            if (comptime std.mem.eql(u8, field.name, field_to_keep)) {
                field_info = field;
                temp_new_fields[real_num_fields] = field_info;
                real_num_fields += 1;
                break;
            }
        }
    }

    comptime var new_fields: [real_num_fields]std.builtin.Type.StructField = undefined;
    inline for (0..real_num_fields) |i| {
        new_fields[i] = temp_new_fields[i];
    }

    return @Type(.{ .@"struct" = .{
        .layout = .auto,
        .fields = &new_fields,
        .decls = &.{},
        .is_tuple = false,
    } });
}

/// Creates a new struct type by keeping only the specified fields from the
/// original struct type.
///
/// If a field specified in `fields_to_keep` does not exist in `s`, it is
/// ignored. Only the fields listed in `fields_to_keep` will be present in the
/// new struct type.
///
/// ## Arguments
/// * `s` (`anytype`): The struct instance from which to keep fields.
/// * `fields_to_keep` (`[]const []const u8`): An array of field names to keep
///   in the new struct type.
///
/// ## Returns
/// `types.KeepStructFields(@TypeOf(s), fields_to_keep)`: A new struct instance
/// containing only the specified fields.
///
/// ## Examples
/// ```zig
/// const s = .{ .a = 1, .b = 2, .c = 3, .d = 4 };
///
/// const kept = zml.types.keepStructFields(s, &.{ "a", "c" });
/// // kept is .{ .a = 1, .c = 3 }
/// ```
pub fn keepStructFields(s: anytype, comptime fields_to_keep: []const []const u8) KeepStructFields(@TypeOf(s), fields_to_keep) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("zml.types.keepStructFields: s must be a struct, got\n\ts: " ++ @typeName(S) ++ "\n");

    var result: KeepStructFields(S, fields_to_keep) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}
