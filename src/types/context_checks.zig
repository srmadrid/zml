const std = @import("std");

const types = @import("../types.zig");

/// Validates a context struct against a specification struct.
///
/// This function checks that the context struct matches the specification
/// struct, ensuring that all required fields are present and have the correct
/// types. It also checks for unexpected fields in the context struct.
///
/// Parameters
/// ----------
/// ctx (`anytype`): The context struct to validate. Must be a struct type.
/// spec (`anytype`): The specification struct that defines the expected fields
/// and their types. Must have the following structure:
/// ```zig
/// .{
///     .field_name = .{ .type = type, .required = bool },
///     ...
/// }
/// ```
/// where `field_name` is the name of the field, `type` is the expected type of
/// the field, and `required` is a boolean indicating whether the field is
/// required or not. If `required` is `false`, then another field named
/// `default` must be present, specifying the default value for the field.
///
/// Returns
/// -------
/// `void`: If the context struct is valid according to the specification.
///
/// Notes
/// -----
/// If the specification struct has no fields, the context struct is allowed to
/// be `void`.
pub fn validateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (specinfo.@"struct".fields.len == 0 and Ctx == void)
        return; // No context needed, nothing to validate

    if (ctxinfo != .@"struct")
        @compileError("Expected struct for context, got " ++ @typeName(Ctx));

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = @hasField(@TypeOf(field_spec), "required") and field_spec.required;
        const expected_type = field_spec.type;

        if (@hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);
            const type_info = @typeInfo(expected_type);
            const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
                break :blk type_info == .@"enum" or (type_info == .optional and @typeInfo(type_info.optional.child) == .@"enum");
            } else if (type_info == .optional)
                actual_type == type_info.optional.child or actual_type == @TypeOf(null)
            else
                actual_type == expected_type;

            if (!types_match)
                @compileError(formatSpecCtxMismatch(Ctx, spec));
        } else if (required) {
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
    const ctxinfo = @typeInfo(Ctx);

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
        inline for (specinfo.@"struct".fields) |field| {
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

    return std.fmt.comptimePrint(
        "Context struct must have the following structure:\n{s}\nGot:\n{s}",
        .{ spec_str, ctx_str },
    );
}

/// Validates that a context struct contains all required fields as defined
/// by a specification struct.
///
/// This function differs from `validateContext` in that it does not raise
/// a compile error if the context struct contains unexpected fields.
///
/// Parameters
/// ----------
/// ctx (`anytype`): The context struct to validate. Must be a struct type.
///
/// spec (`anytype`): The specification struct that defines the expected fields
/// and their types. Must have the following structure:
/// ```zig
/// .{
///     .field_name = .{ .type = type, .required = bool },
///     ...
/// }
/// ```
/// where `field_name` is the name of the field, `type` is the expected type of
/// the field, and `required` is a boolean indicating whether the field is
/// required or not.
///
/// Returns
/// -------
/// `void`: If the context struct is valid according to the specification.
pub fn partialValidateContext(comptime Ctx: type, comptime spec: anytype) void {
    const ctxinfo = @typeInfo(Ctx);

    const SpecType = @TypeOf(spec);
    const specinfo = @typeInfo(SpecType);

    if (specinfo.@"struct".fields.len == 0 and Ctx == void)
        return; // No context needed, nothing to validate

    if (ctxinfo != .@"struct")
        @compileError("Expected struct for context, got " ++ @typeName(Ctx));

    // Check that all required fields exist and types match
    inline for (specinfo.@"struct".fields) |field| {
        const field_spec = @field(spec, field.name);
        const required = @hasField(@TypeOf(field_spec), "required") and field_spec.required;
        const expected_type = field_spec.type;

        if (@hasField(Ctx, field.name)) {
            const actual_type = @FieldType(Ctx, field.name);
            const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
                const type_info = @typeInfo(expected_type);
                break :blk type_info == .@"enum" or (type_info == .optional and @typeInfo(type_info.optional.child) == .@"enum");
            } else actual_type == expected_type;

            if (!types_match)
                @compileError(formatCtxRequiredFields(Ctx, spec));
        } else if (required) {
            @compileError(formatCtxRequiredFields(Ctx, spec));
        }
    }
}

fn formatCtxRequiredFields(
    comptime Ctx: type,
    comptime spec: anytype,
) []const u8 {
    const ctxinfo = @typeInfo(Ctx);

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
        inline for (specinfo.@"struct".fields) |field| {
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
                    if (required) "required" else "optional",
                },
            );
        }
        spec_str = spec_str ++ "}\n";
    }

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

    return std.fmt.comptimePrint(
        "Context struct must contain at least the following fields:\n{s}\nGot:\n{s}",
        .{ spec_str, ctx_str },
    );
}

pub fn ctxHasField(
    comptime T: type,
    comptime field_name: []const u8,
    comptime FieldType: type,
) bool {
    const has_field: bool = @hasField(T, field_name);

    if (has_field) {
        if (@FieldType(T, field_name) == FieldType or
            (@FieldType(T, field_name) == @TypeOf(.enum_literal) and
                (@typeInfo(FieldType) == .@"enum" or
                    (@typeInfo(FieldType) == .optional and @typeInfo(@typeInfo(FieldType).optional.child) == .@"enum"))) or
            (@FieldType(T, field_name) == @TypeOf(null) and
                @typeInfo(FieldType) == .optional))
            return true;
    }

    return false;
}

pub fn getFieldOrDefault(ctx: anytype, comptime spec: anytype, comptime field_name: []const u8) @field(spec, field_name).type {
    const T = @TypeOf(ctx);
    const FieldType = @field(spec, field_name).type;

    if (@hasField(T, field_name)) {
        const actual_type = @FieldType(T, field_name);
        const type_info = @typeInfo(FieldType);
        const expected_type = FieldType;
        const types_match = if (actual_type == @TypeOf(.enum_literal)) blk: { // Special case for enum literals
            break :blk type_info == .@"enum" or (type_info == .optional and @typeInfo(type_info.optional.child) == .@"enum");
        } else if (type_info == .optional)
            actual_type == type_info.optional.child or actual_type == @TypeOf(null)
        else
            actual_type == expected_type;

        if (!types_match)
            @compileError("Field '" ++ field_name ++ "' has type " ++ @typeName(@FieldType(T, field_name)) ++ ", expected " ++ @typeName(FieldType));

        return @field(ctx, field_name);
    }

    return @field(spec, field_name).default;
}

pub fn MixStructs(comptime S1: type, comptime S2: type) type {
    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("MixStructs: both types must be structs");

    for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("Field '" ++ field.name ++ "' already exists in the second struct");
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

pub fn mixStructs(s1: anytype, s2: anytype) MixStructs(@TypeOf(s1), @TypeOf(s2)) {
    const S1 = @TypeOf(s1);
    const S2 = @TypeOf(s2);

    const info1 = @typeInfo(S1);
    const info2 = @typeInfo(S2);

    if (info1 != .@"struct" or info2 != .@"struct")
        @compileError("mixStructs: both types must be structs");

    var result: MixStructs(S1, S2) = undefined;
    inline for (info1.@"struct".fields) |field| {
        if (@hasField(S2, field.name))
            @compileError("Field '" ++ field.name ++ "' already exists in the second struct");

        @field(result, field.name) = @field(s1, field.name);
    }
    inline for (info2.@"struct".fields) |field| {
        @field(result, field.name) = @field(s2, field.name);
    }

    return result;
}

pub fn StripStruct(comptime S: type, comptime fields_to_remove: []const []const u8) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

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

pub fn stripStruct(s: anytype, comptime fields_to_remove: []const []const u8) StripStruct(@TypeOf(s), fields_to_remove) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    var result: StripStruct(S, fields_to_remove) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}

pub fn RenameStructFields(comptime S: type, comptime fields_to_rename: anytype) type {
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    const F: type = @TypeOf(fields_to_rename);
    const finfo = @typeInfo(F);

    if (finfo != .@"struct")
        @compileError("fields_to_rename must be a struct");

    // Check that all new names do not exist in the original struct
    inline for (finfo.@"struct".fields) |field| {
        if (@hasField(S, @field(fields_to_rename, field.name)))
            @compileError("Field '" ++ @field(fields_to_rename, field.name) ++ "' already exists in the struct");
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

/// Renames fields of a struct according to a mapping provided in another struct.
///
/// If a field specified in `fields_to_rename` does not exist in `s`, it is
/// ignored. Fields in `s` that are not specified in `fields_to_rename` retain
/// their original names.
///
/// Parameters
/// ----------
/// `s` (`anytype`): The struct instance whose fields are to be renamed.
/// Must be a struct type.
///
/// `fields_to_rename` (`anytype`): A struct instance that defines the mapping
/// of old field names to new field names. Each field in this struct should
/// have the name of the field to be renamed in `s`, and its value
/// should be a string representing the new name. I.e.:
/// ```zig
/// .{
///     .old_field_name1 = "new_field_name1",
///     .old_field_name2 = "new_field_name2",
///     ...
/// }
/// ```
///
/// Returns
/// -------
/// The new struct instance with fields renamed according to the mapping.
pub fn renameStructFields(s: anytype, comptime fields_to_rename: anytype) RenameStructFields(@TypeOf(s), fields_to_rename) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

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
        @compileError("Type must be a struct");

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
/// The fields specified in `fields_to_keep` need not exist in the original
/// struct; any non-existing fields will simply be ignored.
///
/// Parameters
/// ----------
/// `s` (`anytype`): The struct instance from which to keep fields. Must be a
/// struct type.
///
/// `fields_to_keep` (`[]const []const u8`): An array of field names to keep
/// in the new struct type.
///
/// Returns
/// -------
/// The new struct instance containing only the specified fields.
pub fn keepStructFields(s: anytype, comptime fields_to_keep: []const []const u8) KeepStructFields(@TypeOf(s), fields_to_keep) {
    const S = @TypeOf(s);
    const info = @typeInfo(S);
    if (info != .@"struct")
        @compileError("Type must be a struct");

    var result: KeepStructFields(S, fields_to_keep) = undefined;
    inline for (@typeInfo(@TypeOf(result)).@"struct".fields) |field| {
        @field(result, field.name) = @field(s, field.name);
    }

    return result;
}
