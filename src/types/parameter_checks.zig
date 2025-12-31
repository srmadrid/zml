//! Namespace for "anytype" parameter checks used in parameter type validation.

const std = @import("std");

const types = @import("../types.zig");

const max_domains_checked = 32;
const max_parameters = 32;

const CheckedDomain = struct {
    first_level: FirstLevel, // The first level type (domain, numeric, vector, matrix_kind, matrix_storage, matrix, array)
    element_type: ?[]const u8, // The element type expression, if any (e.g., "float" in "matrix.dense(float)")

    const FirstLevel = union(enum) {
        domain: types.Domain,
        numeric: []const u8, // The numeric type expression as a string since we allow more complex expressions here
        vector: types.VectorType,
        matrix_kind: types.MatrixKind,
        matrix_storage: types.MatrixStorage,
        matrix: types.MatrixType,
        array: types.ArrayType,

        pub fn fromAny(comptime s: anytype) FirstLevel {
            return switch (@TypeOf(s)) {
                types.Domain => .{ .domain = s },
                types.NumericType => .{ .numeric = s },
                types.VectorType => .{ .vector = s },
                types.MatrixKind => .{ .matrix_kind = s },
                types.MatrixStorage => .{ .matrix_storage = s },
                types.MatrixType => .{ .matrix = s },
                types.ArrayType => .{ .array = s },
                else => unreachable,
            };
        }

        pub fn toString(self: FirstLevel) []const u8 {
            return switch (self) {
                .domain => self.domain.toString(),
                .numeric => "numeric",
                .vector => self.vector.toString(),
                .matrix_kind => self.matrix_kind.toString(),
                .matrix_storage => self.matrix_storage.toString(),
                .matrix => self.matrix.toString(),
                .array => self.array.toString(),
            };
        }
    };
};

fn mixSubdomains(comptime base: types.Domain, kind: anytype, storage: anytype) switch (base) {
    .matrix => types.MatrixType,
    else => unreachable,
} {
    switch (base) {
        .matrix => {
            const K = @TypeOf(kind);
            const S = @TypeOf(storage);

            if (K != types.MatrixKind or S != types.MatrixStorage) {
                @compileError("For base domain 'matrix', kind and storage must be of types MatrixKind and MatrixStorage, respectively.");
            }

            switch (kind) {
                .general => switch (storage) {
                    .dense => return .dense_general,
                    .banded => return .banded_general,
                    .tridiagonal => return .tridiagonal_general,
                    .sparse => return .sparse_general,
                    .block => return .block_general,
                },
                .symmetric => switch (storage) {
                    .dense => return .dense_symmetric,
                    .banded => return .banded_symmetric,
                    .tridiagonal => return .tridiagonal_symmetric,
                    .sparse => return .sparse_symmetric,
                    .block => return .block_symmetric,
                },
                .hermitian => switch (storage) {
                    .dense => return .dense_hermitian,
                    .banded => return .banded_hermitian,
                    .tridiagonal => return .tridiagonal_hermitian,
                    .sparse => return .sparse_hermitian,
                    .block => return .block_hermitian,
                },
                .triangular => switch (storage) {
                    .dense => return .dense_triangular,
                    .banded => unreachable,
                    .tridiagonal => unreachable,
                    .sparse => return .sparse_triangular,
                    .block => unreachable,
                },
                .diagonal => return .diagonal,
                .permutation => return .permutation,
            }
        },
        else => unreachable,
    }
}

const CategoryOperator = enum {
    /// bool, int, integer
    integral,
    /// float, dyadic, cfloat, rational, real, complex
    nonintegral,
    /// bool, int, float, dyadic, cfloat
    fixed,
    /// integer, rational, real, complex
    arbitrary,
    /// bool, int, float, dyadic, integer, rational, real
    real,
    /// cfloat, complex
    complex,
    /// signed int, float, dyadic, cfloat, integer, rational, real, complex
    signed,
    /// bool, unsigned int
    unsigned,

    pub fn match(self: CategoryOperator) fn (type) bool {
        return switch (self) {
            .integral => return types.isIntegral,
            .nonintegral => return types.isNonIntegral,
            .fixed => return types.isFixedPrecision,
            .arbitrary => return types.isArbitraryPrecision,
            .real => return types.isReal,
            .complex => return types.isComplex,
            .signed => return types.isSigned,
            .unsigned => return types.isUnsigned,
        };
    }
};

/// Given a type signature (as a string) and a type, checks whether the type
/// matches the signature, triggering a compile-time error if not.
///
/// Parameters
/// ----------
/// comptime mode (`enum { silent, compile_error }`): The mode of operation. If
/// set to `silent`, the function returns a boolean indicating whether the type
/// matches the signature. If set to `compile_error`, a compile-time error is
/// triggered if the type does not match the signature.
///
/// comptime type_signature (`[]const u8`): The type signature to check against.
/// The type signature must start with:
///
/// - "": pass by value
/// - "*": mutable one-item pointer
/// - "* const": constant one-item pointer
/// - "[*]": mutable many-item pointer
/// - "[*] const": constant many-item pointer
/// - "[]": mutable slice
/// - "[] const": constant slice
///
/// Then, it must specify the type:
///
/// - "numeric": any supported numeric type. For specific numeric types, use the
///   names defined in the `NumericType` enum.
/// - "vector": any vector type. For specific vector types, use "vector.<storage>",
///   where `<storage>` is one of:
///
///     - "dense": dense vector
///     - "sparse": sparse vector
///
/// - "matrix": any matrix type. For specific matrix types, use "matrix.<kind>.<storage>",
///   where `<kind>` is one of:
///
///     - "general": general matrix
///     - "symmetric": symmetric matrix
///     - "hermitian": hermitian matrix
///     - "triangular": triangular matrix
///     - "diagonal": diagonal matrix (<storage> is ignored)
///     - "permutation": permutation matrix (<storage> is ignored)
///
///   and `<storage>` is one of:
///
///     - "dense": dense matrix
///     - "banded": banded matrix
///     - "tridiagonal": tridiagonal matrix
///     - "sparse": sparse matrix
///     - "block": block sparse matrix
///   Setting only the kind (e.g., "matrix.symmetric") or only the storage
///   (e.g., "matrix.sparse") is also accepted, allowing for more general checks,
///   otherwise the order must be `<kind>.<storage>`.
/// - "array": any array type. For specific array types, use "array.<storage>",
///   where `<storage>` is one of:
///     - "dense": dense array
///     - "strided": strided array
///     - "sparse": sparse array
/// - "expression": just use `Expression` as the type, why are you even checking this?
/// - "any": any type is accepted
///
/// In any container type (vector, matrix, array, etc.), after the type name and
/// without spaces, you can add a numeric type expression in angle brackets to
/// specify the element type. For example, "matrix.dense(float)" specifies a
/// dense matrix with float elements. If no element type is specified, any
/// element type is accepted. Inside the angle brackets, the same rules as for
/// numeric types apply, including the more granular constraints described below.
/// For more granular control over numeric types, constraint operators can be
/// used. Any operator <op> must be placed immediately before the type it
/// applies to, and with no spaces, like <op><type>. We define three kinds of
/// operators, principal operators (which can only be used once or twice at the
/// beginning of the numeric type expression), additional operators (which can
/// be used multiple times and in any order after the principal operators), and
/// final operators (which must be at the end of the expression and apply to
/// everything before them). The principal operators are:
///
/// - Ordering operators (only for numeric types, following the order in the
///   `NumericType` enum):
///     - "<": only for numeric types, accepts types lower than the specified type.
///     - "<=": only for numeric types, accepts types lower than or equal to the
///     specified type.
///     - ">": only for numeric types, accepts types higher than the specified type.
///     - ">=": only for numeric types, accepts types higher than or equal to the
///     specified type.
/// - "=": accepts only the specified type. It is equivalent to not using any operator
///
/// but may speed up compilation in some cases.
/// If two ordering operators are used, the first must be a lower bound ("<" or "<=")
/// and the second an upper bound (">" or ">="). For example, "<=dyadic >=int" accepts
/// any numeric type between int and float, inclusive (in this case, ints, floats and
/// dyadics). The additional operators can also be used for different domains. However,
/// if used in this way and a numeric type is specified as a base type (not inside a container),
/// the numeric part must be at the end of the signature. For example,
/// "matrix.sparse(@real) <op> matrix.block(@fixed @real) <op> >=int <=dyadic" is valid,
/// but ">=int <=dyadic <op> matrix.sparse(@real) <op> matrix.block(@fixed @real)" is not.
/// The additional operators are:
///
/// - "!": exclusion operator, for excluding types. For example, "numeric !float !complex"
/// accepts any numeric type except float and complex types. The principal part
/// must be "numeric" or a range when using exclusion operators on numeric types, or
/// "any" or a general container type (vector, matrix, array) when excluding container types,
/// for example, "matrix !matrix.sparse" accepts any matrix type except sparse matrices, or
/// "any !numeric" accepts any type except numeric types.
/// - "|": inclusion operator, for including only specific types. The principal
/// part must be a specific type or a range when using inclusion operators.
/// For instance, ">=integer |int" accepts anything from integer upwards, but
/// also int (which is below integer). For container types this works similarly,
/// e.g., "matrix |vector" accepts any matrix or vector type.
///
/// The final operators are all category operators, and are prefixed with "@":
///
/// - "@integral": accepts only integral types (bool, int, and integer).
/// - "@nonintegral": accepts only non-integral types (rational, real, float, dyadic, cfloat, complex, and expression).
/// - "@fixed": accepts only fixed precision types (bool, int, float, dyadic, and cfloat).
/// - "@arbitrary": accepts only arbitrary precision types (integer, rational,
/// real, and complex).
/// - "@real": accepts only real types (bool, int, integer, rational, real, float).
/// - "@complex": accepts only complex types (cfloat, complex).
/// - "@signed": accepts only signed types (int, integer, rational, real, float, cfloat,
/// complex).
/// - "@unsigned": accepts only unsigned types (bool and uint).
///
/// Final operators can be combined. For example, "numeric @fixed @real" accepts
/// any fixed precision real type (bool, int, float, dyadic). They also do not need
/// the principal part, so "@real" alone is valid and accepts any real numeric type.
/// However, ain important constraint is that final operators apply to everything else
/// within the same numeric type expression. For example, "@real !float |complex"
/// still excludes complex types, because the category operator is applied last.
/// The type signature is case-insensitive and ignores repeated spaces. Below are
/// examples of valid type signatures. We divide by the nature of the examples:
///
/// - Basic signatures (pointer prefixes):
///     - "any": any non-pointer type
///     - "numeric": any numeric type
///     - "* numeric": mutable one-item pointer to any numeric type
///     - "* const numeric": constant one-item pointer to any numeric type
///     - "[] numeric": mutable slice of any numeric type
///     - "* any": mutable one-item pointer to any type
/// - Basic numeric expressions:
///     - "=float": only float, equivalent to just "float"
///     - ">=int": int or any type higher than int
///     - ">=int <=dyadic": any type between int and dyadic, inclusive
///     - "numeric !float !complex": any numeric type except float and complex types
///     - ">=integer |int": anything from integer upwards, but also int
///     - "numeric @real": any real numeric type (not cfloat or complex)
///     - ">=int @signed": any signed type equal to or higher than int
///     - "numeric @fixed @real": any fixed precision real type (bool, int, float, dyadic)
/// - Basic container types:
///     - "vector": any vector type holding any element type
///     - "matrix.symmetric.": any symmetric matrix type holding any element type
///     - "matrix.general.dense": a `matrix.general.Dense` holding any element type
/// - Container types with numeric expressions:
///     - "vector(float)": any vector type holding float elements
///     - "matrix(@complex)": any matrix type holding complex element types
/// - Complex combinations:
///     - "[] matrix.sparse(@real)": mutable slice of sparse matrices holding real element types
///     - "* const matrix.symmetric(rational |integer)": constant one-item pointer to symmetric matrices
///     holding rational or integer element types
///     - "matrix.general.dense<numeric !float @signed>": a `matrix.general.Dense`
///     holding any signed numeric element type except float
///     - "matrix |vector": any matrix or vector type holding any element type
///
/// comptime T (`type`): The type to check against the signature.
///
/// comptime fn_name (`[]const u8`): The name of the function from which this
/// check is being called. Used for error messages.
///
/// comptime param_name (`[]const u8`): The name of the parameter being checked.
/// Used for error messages.
///
/// Returns
/// -------
/// `void`
pub fn checkParameterType(
    comptime mode: enum { silent, compile_error }, // silent -> non-matching types are ignored, compile_error -> trigger compile error on non-matching types. Errors in the signature itself always trigger compile errors.
    comptime type_signature: []const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
) switch (mode) {
    .silent => bool,
    .compile_error => void,
} {
    comptime var VT: type = T;

    // Preprocess type signature: lowercase and trim leading spaces
    comptime var signature_array: [type_signature.len + 1]u8 = undefined;
    for (type_signature, 0..) |c, i|
        signature_array[i] = std.ascii.toLower(c);
    signature_array[type_signature.len] = 0;
    comptime var signature: []const u8 = signature_array[0..signature_array.len];

    // Pointer or slice prefix matching
    signature = std.mem.trimStart(u8, signature, " ");

    comptime var pointer: ?std.builtin.Type.Pointer.Size = null;
    comptime var constant: bool = false;
    switch (signature[0]) {
        '*' => {
            pointer = .one;

            if (signature.len >= 7 and
                std.mem.eql(u8, signature[1..7], " const"))
            {
                // constant one-item pointer
                if (!types.isPointer(VT)) // any non-const pointer can be coerced to const
                    @compileError(
                        std.fmt.comptimePrint(
                            "{s} requires parameter '{s}' to be a one-item pointer, but got '{s}'.",
                            .{ fn_name, param_name, @typeName(T) },
                        ),
                    );

                constant = true;
                VT = types.Child(VT);
                signature = signature[7..];
            } else {
                // mutable one-item pointer
                if (!types.isPointer(VT) or types.isConstPointer(VT))
                    @compileError(
                        std.fmt.comptimePrint(
                            "{s} requires parameter '{s}' to be a mutable one-item pointer, but got '{s}'.",
                            .{ fn_name, param_name, @typeName(T) },
                        ),
                    );

                VT = types.Child(VT);
                signature = signature[1..];
            }
        },
        '[' => {
            switch (signature[1]) {
                '*' => {
                    pointer = .many;

                    if (signature.len >= 8 and
                        std.mem.eql(u8, signature[2..8], "] const"))
                    {
                        // constant many-item pointer
                        if (!types.isManyPointer(VT)) // any non-const pointer can be coerced to const
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a many-item pointer, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = types.Child(VT);
                        signature = signature[8..];
                    } else {
                        // mutable many-item pointer
                        if (!types.isManyPointer(VT) or types.isConstPointer(VT))
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a mutable many-item pointer, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = types.Child(VT);
                        signature = signature[2..];
                    }
                },
                ']' => {
                    pointer = .slice;

                    if (signature.len >= 7 and
                        std.mem.eql(u8, signature[2..7], " const"))
                    {
                        // constant slice
                        if (!types.isSlice(VT))
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a slice, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = types.Child(VT);
                        signature = signature[7..];
                    } else {
                        // mutable slice
                        if (!types.isSlice(VT) or types.isConstPointer(VT))
                            @compileError(
                                std.fmt.comptimePrint(
                                    "{s} requires parameter '{s}' to be a mutable slice, but got '{s}'.",
                                    .{ fn_name, param_name, @typeName(T) },
                                ),
                            );

                        VT = types.Child(VT);
                        signature = signature[2..];
                    }
                },
                else => @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized pointer type in type signature '{s}' in parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                ),
            }
        },
        else => {},
    }

    signature = std.mem.trimStart(u8, signature, " ");

    // Actual type checking
    comptime var matched_any_positive = false; // Must match at least one positive type ("" or "|")
    comptime var matched_any_negative = false; // Must not match any negative type ("!")
    comptime var finished = false;
    comptime var checked = 0; // Number of principal checks done
    comptime var domains_checked: [max_domains_checked]?CheckedDomain = .{null} ** max_domains_checked; // For compilation error messages
    comptime var positive_operators: [max_domains_checked]bool = .{false} ** max_domains_checked;
    while (!finished) : (checked += 1) {
        signature = std.mem.trimStart(u8, signature, " ");

        if (signature.len == 0 or signature[0] == 0) {
            finished = true;
            break;
        }

        // Check for "|" and "!" operators
        comptime var positive_operator = true;
        if (signature.len >= 1 and signature[0] == '|') {
            if (checked == 0) {
                @compileError(
                    std.fmt.comptimePrint(
                        "The '|' operator cannot be used as the first operator in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );
            }

            positive_operator = true;
            signature = signature[1..];
        } else if (signature.len >= 1 and signature[0] == '!') {
            if (checked == 0) {
                @compileError(
                    std.fmt.comptimePrint(
                        "The '!' operator cannot be used as the first operator in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );
            }

            positive_operator = false;
            signature = signature[1..];
        }

        signature = std.mem.trimStart(u8, signature, " ");

        positive_operators[checked] = positive_operator;

        // Check for "any". Only valid if nothing else has been checked yet
        if (std.mem.eql(u8, signature[0..3], "any")) {
            if (checked == 0) {
                // "any" matches everything supported by the library, ignore rest of signature
                if (!types.isSupportedType(VT))
                    @compileError(
                        std.fmt.comptimePrint(
                            "{s} requires parameter '{s}' to be any supported type, but got '{s}'.",
                            .{ fn_name, param_name, @typeName(T) },
                        ),
                    )
                else
                    return; // Matched anything, ignore rest of signature
            } else {
                @compileError(
                    std.fmt.comptimePrint(
                        "The 'any' type can only be used as the first principal type, in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );
            }
        }

        comptime var matched_domain: types.Domain = .numeric;

        // Check for domain types
        inline for (std.meta.fields(types.Domain)) |domain| {
            if (domain.name.len <= signature.len and
                std.mem.eql(u8, signature[0..domain.name.len], domain.name))
            {
                matched_domain = @enumFromInt(domain.value);

                if (matched_domain == .numeric) {
                    // Numeric types are handled later
                    continue;
                }

                // Consume matched part
                domains_checked[checked] = .{ .first_level = .{ .domain = matched_domain }, .element_type = null };
                signature = signature[domain.name.len..];

                break;
            }
        }

        if (matched_domain == .numeric) {
            // If numeric domain, either we matched numeric or nothing yet
            const numeric_constraint = checkNumericConstraints(
                type_signature,
                signature,
                T,
                fn_name,
                param_name,
                &domains_checked[checked],
            );

            if (positive_operator)
                matched_any_positive = matched_any_positive or numeric_constraint
            else
                matched_any_negative = matched_any_negative or numeric_constraint;

            // Numeric domain must be the only one or the last one
            finished = true;
        } else {
            // Need to validate domain match

            // Domain refinement
            switch (matched_domain) {
                .vector => {
                    // Vector type refinement, can only have one level
                    const one_level_check = oneLevelCheck(
                        type_signature,
                        &signature,
                        VT,
                        fn_name,
                        param_name,
                        .vector,
                        types.VectorType,
                        &domains_checked[checked],
                    );

                    if (positive_operator)
                        matched_any_positive = matched_any_positive or one_level_check
                    else
                        matched_any_negative = matched_any_negative or one_level_check;
                },
                .matrix => {
                    // Matrix type refinement
                    const two_level_check = twoLevelCheck(
                        type_signature,
                        &signature,
                        VT,
                        fn_name,
                        param_name,
                        .matrix,
                        types.MatrixKind,
                        types.MatrixStorage,
                        &domains_checked[checked],
                    );

                    if (positive_operator)
                        matched_any_positive = matched_any_positive or two_level_check
                    else
                        matched_any_negative = matched_any_negative or two_level_check;
                },
                .array => {
                    // Array type refinement, can only have one level
                    const one_level_check = oneLevelCheck(
                        type_signature,
                        &signature,
                        VT,
                        fn_name,
                        param_name,
                        .array,
                        types.ArrayType,
                        &domains_checked[checked],
                    );

                    if (positive_operator)
                        matched_any_positive = matched_any_positive or one_level_check
                    else
                        matched_any_negative = matched_any_negative or one_level_check;
                },
                else => unreachable,
            }
        }
    }

    if (mode == .silent) {
        return matched_any_positive and !matched_any_negative;
    } else {
        if (!matched_any_positive or matched_any_negative) {
            // Build error message, can be improved

            // For now, just print the type signature
            comptime var error_message: [1024]u8 = .{0} ** 1024;
            comptime var error_message_len = 0;

            const beginning = std.fmt.comptimePrint(
                "{s} requires parameter '{s}' to be a", // a matrix, a vector.Dense, a dense matrix(float) (if not both levels used)
                .{ fn_name, param_name },
            );

            @memcpy(
                error_message[0..beginning.len],
                beginning,
            );
            error_message_len += beginning.len;

            comptime var negated_message: [1024]u8 = .{0} ** 1024;
            comptime var negated_message_len = 0;

            comptime var i = 0;
            while (domains_checked[i]) |checked_domain| : (i += 1) {
                if (positive_operators[i]) {
                    if (i > 0) {
                        const or_str = " or ";
                        @memcpy(
                            error_message[error_message_len .. error_message_len + or_str.len],
                            or_str,
                        );
                        error_message_len += or_str.len;
                    }

                    switch (checked_domain.first_level) {
                        .numeric => {
                            const part =
                                std.fmt.comptimePrint(
                                    "'{s}'",
                                    .{checked_domain.first_level.numeric},
                                );

                            @memcpy(
                                error_message[error_message_len .. error_message_len + part.len],
                                part,
                            );
                            error_message_len += part.len;

                            continue;
                        },
                        else => {},
                    }

                    const part =
                        std.fmt.comptimePrint(
                            " '{s}",
                            .{checked_domain.first_level.toString()},
                        );

                    @memcpy(
                        error_message[error_message_len .. error_message_len + part.len],
                        part,
                    );
                    error_message_len += part.len;

                    if (checked_domain.element_type) |et| {
                        const elem_part =
                            std.fmt.comptimePrint(
                                "({s})'",
                                .{et},
                            );

                        @memcpy(
                            error_message[error_message_len .. error_message_len + elem_part.len],
                            elem_part,
                        );
                        error_message_len += elem_part.len;
                    } else {
                        // Put any
                        const any_str = "(any)'";
                        @memcpy(
                            error_message[error_message_len .. error_message_len + any_str.len],
                            any_str,
                        );
                        error_message_len += any_str.len;
                    }
                } else {
                    const part =
                        std.fmt.comptimePrint(
                            " '{s}",
                            .{checked_domain.first_level.toString()},
                        );

                    @memcpy(
                        negated_message[negated_message_len .. negated_message_len + part.len],
                        part,
                    );
                    negated_message_len += part.len;

                    if (checked_domain.element_type) |et| {
                        const elem_part =
                            std.fmt.comptimePrint(
                                "({s})'",
                                .{et},
                            );

                        @memcpy(
                            negated_message[negated_message_len .. negated_message_len + elem_part.len],
                            elem_part,
                        );
                        negated_message_len += elem_part.len;
                    } else {
                        // Put any
                        const any_str = "(any)'";
                        @memcpy(
                            negated_message[negated_message_len .. negated_message_len + any_str.len],
                            any_str,
                        );
                        negated_message_len += any_str.len;
                    }
                }
            }

            if (negated_message_len > 0) {
                const excluding_str = std.fmt.comptimePrint(
                    ", excluding {s}",
                    .{negated_message[0..negated_message_len]},
                );

                @memcpy(
                    error_message[error_message_len .. error_message_len + excluding_str.len],
                    excluding_str,
                );
                error_message_len += excluding_str.len;
            }

            const ending =
                std.fmt.comptimePrint(
                    ", but got '{s}'.",
                    .{@typeName(T)},
                );

            @memcpy(
                error_message[error_message_len .. error_message_len + ending.len],
                ending,
            );
            error_message_len += ending.len;

            @compileError(error_message[0..error_message_len]);
        }
    }
}

fn oneLevelCheck(
    comptime type_signature: []const u8,
    comptime signature: *[]const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
    comptime domain: types.Domain,
    comptime DomainType: type,
    comptime checked_domain: *?CheckedDomain,
) bool {
    comptime var type_base_name: [32:0]u8 = .{0} ** 32;
    @memcpy(type_base_name[0..domain.toString().len], domain.toString());
    comptime var filled_length = domain.toString().len;

    comptime var matched_subtype: ?DomainType = null;
    if (signature.*[0] == '.') {
        signature.* = signature.*[1..]; // Consume "."
        inline for (std.meta.fields(DomainType)) |subtype| {
            if (subtype.name.len <= signature.len and
                std.mem.eql(u8, signature.*[0..subtype.name.len], subtype.name))
            {
                // Consume matched part
                signature.* = signature.*[subtype.name.len..];
                matched_subtype = @enumFromInt(subtype.value);

                break;
            }
        }

        if (matched_subtype == null)
            @compileError(
                std.fmt.comptimePrint(
                    "Unrecognized {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
                ),
            );

        type_base_name[filled_length] = '.';
        const tag_name = @tagName(matched_subtype.?);
        type_base_name[filled_length + 1] = std.ascii.toUpper(tag_name[0]);
        filled_length += 2;
        @memcpy(type_base_name[filled_length .. filled_length + tag_name.len - 1], tag_name[1..]);
        filled_length += tag_name.len - 1;
    }

    // Check for element type specification
    comptime var element_type_specified: bool = false;
    comptime var correct_element_type: bool = true; // Default to true for no specification
    comptime var element_type_signature: []const u8 = "";
    if (signature.len > 0 and signature.*[0] == '(') {
        const closing_paren_index = std.mem.indexOfScalar(u8, signature.*, ')');
        if (closing_paren_index == null) {
            @compileError(
                std.fmt.comptimePrint(
                    "Unmatched '(' in type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ type_signature, param_name, fn_name },
                ),
            );
        }

        element_type_specified = true;
        element_type_signature = signature.*[1..closing_paren_index.?];

        // Check element type
        correct_element_type = checkNumericConstraints(
            type_signature,
            element_type_signature,
            types.Numeric(T),
            fn_name,
            param_name,
            checked_domain,
        );

        // Consume element type part, including parentheses
        signature.* = signature.*[closing_paren_index.? + 1 ..];

        checked_domain.*.?.element_type = element_type_signature;
    }

    // Validate type match
    comptime var correct_subdomain_type: bool = false;
    if (matched_subtype) |subtype| {
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(subtype);
        correct_subdomain_type = subtype.match()(T);
    } else {
        checked_domain.*.?.first_level = .{ .domain = domain };
        correct_subdomain_type = domain.match()(T);
    }

    return correct_subdomain_type and correct_element_type;
}

fn twoLevelCheck(
    comptime type_signature: []const u8,
    comptime signature: *[]const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
    comptime domain: types.Domain,
    comptime DomainType1: type, // Can be skipped if only one level is present
    comptime DomainType2: type,
    comptime checked_domain: *?CheckedDomain,
) bool {
    comptime var type_base_name: [32:0]u8 = .{0} ** 32;
    @memcpy(type_base_name[0..domain.toString().len], domain.toString());
    comptime var filled_length = domain.toString().len;

    comptime var matched_subtype1: ?DomainType1 = null;
    comptime var matched_subtype2: ?DomainType2 = null;

    // First subtype
    if (signature.*[0] == '.') {
        signature.* = signature.*[1..]; // Consume "."

        // Test first DomainType1
        inline for (std.meta.fields(DomainType1)) |subtype1| {
            if (subtype1.name.len <= signature.len and
                std.mem.eql(u8, signature.*[0..subtype1.name.len], subtype1.name))
            {
                // Consume matched part
                signature.* = signature.*[subtype1.name.len..];
                matched_subtype1 = @enumFromInt(subtype1.value);

                break;
            }
        }

        if (matched_subtype1 == null) {
            // Test second DomainType2 if first failed
            inline for (std.meta.fields(DomainType2)) |subtype2| {
                if (subtype2.name.len <= signature.len and
                    std.mem.eql(u8, signature.*[0..subtype2.name.len], subtype2.name))
                {
                    // Consume matched part
                    signature.* = signature.*[subtype2.name.len..];
                    matched_subtype2 = @enumFromInt(subtype2.value);

                    break;
                }
            }

            if (matched_subtype2 == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
                    ),
                );

            type_base_name[filled_length] = '.';
            const tag_name2 = @tagName(matched_subtype2.?);
            filled_length += 1;
            @memcpy(type_base_name[filled_length .. filled_length + tag_name2.len], tag_name2);
            filled_length += tag_name2.len;
        } else {
            type_base_name[filled_length] = '.';
            const tag_name1 = @tagName(matched_subtype1.?);
            filled_length += 1;
            @memcpy(type_base_name[filled_length .. filled_length + tag_name1.len], tag_name1);
            filled_length += tag_name1.len;
        }
    }

    // Second subtype
    if (matched_subtype2 == null) {
        if (signature.len > 0 and signature.*[0] == '.') {
            signature.* = signature.*[1..]; // Consume "."

            // Second is always DomainType2
            inline for (std.meta.fields(DomainType2)) |subtype2| {
                if (subtype2.name.len <= signature.len and
                    std.mem.eql(u8, signature.*[0..subtype2.name.len], subtype2.name))
                {
                    // Consume matched part
                    signature.* = signature.*[subtype2.name.len..];
                    matched_subtype2 = @enumFromInt(subtype2.value);

                    break;
                }
            }

            if (matched_subtype2 == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
                    ),
                );

            type_base_name[filled_length] = '.';
            const tag_name2 = @tagName(matched_subtype2.?);
            type_base_name[filled_length + 1] = std.ascii.toUpper(tag_name2[0]);
            filled_length += 2;
            @memcpy(type_base_name[filled_length .. filled_length + tag_name2.len - 1], tag_name2[1..]);
            filled_length += tag_name2.len - 1;
        }
    } else if (signature.len > 0 and signature.*[0] == '.') {
        @compileError(
            std.fmt.comptimePrint(
                "Found unexpected second {s} type in type signature '{s}' for parameter '{s}' of function '{s}'.",
                .{ std.mem.span(@as([*:0]const u8, @ptrCast(&type_base_name))), type_signature, param_name, fn_name },
            ),
        );
    }

    // Check for element type specification
    comptime var element_type_specified: bool = false;
    comptime var correct_element_type: bool = true; // Default to true for no specification
    comptime var element_type_signature: []const u8 = "";
    if (signature.len > 0 and signature.*[0] == '(') {
        const closing_paren_index = std.mem.indexOfScalar(u8, signature.*, ')');
        if (closing_paren_index == null) {
            @compileError(
                std.fmt.comptimePrint(
                    "Unmatched '(' in type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ type_signature, param_name, fn_name },
                ),
            );
        }

        element_type_specified = true;
        element_type_signature = signature.*[1..closing_paren_index.?];

        // Check element type
        correct_element_type = checkNumericConstraints(
            type_signature,
            element_type_signature,
            types.Numeric(T),
            fn_name,
            param_name,
            checked_domain,
        );

        // Consume element type part, including parentheses
        signature.* = signature.*[closing_paren_index.? + 1 ..];
    }

    // Validate type match
    comptime var correct_subdomain_type: bool = false;
    if (matched_subtype1 != null and matched_subtype2 != null) {
        // Both active
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(mixSubdomains(domain, matched_subtype1.?, matched_subtype2.?));
        correct_subdomain_type = matched_subtype1.?.match()(T) and matched_subtype2.?.match()(T);
    } else if (matched_subtype1) |subtype1| {
        // Only matched_subtype1 active
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(subtype1);
        correct_subdomain_type = subtype1.match()(T);
    } else if (matched_subtype2) |subtype2| {
        // Only matched_subtype2 active
        checked_domain.*.?.first_level = CheckedDomain.FirstLevel.fromAny(subtype2);
        correct_subdomain_type = subtype2.match()(T);
    } else {
        checked_domain.*.?.first_level = .{ .domain = domain };
        correct_subdomain_type = domain.match()(T);
    }

    return correct_subdomain_type and correct_element_type;
}

/// Given a numeric type signature (as a string) and a numeric type, checks
/// whether the type matches the signature, returning `true` if it does,
/// `false` otherwise. Does not trigger compile-time errors.
///
/// Parameters
/// ----------
/// comptime signature (`[]const u8`): The numeric type signature to check against.
/// See the documentation of `checkParameterType` for the format of type signatures.
///
/// comptime T (`type`): The numeric type to check against the signature.
/// Must be a supported numeric type.
///
/// Returns
/// -------
/// `bool`: `true` if the type matches the signature, `false` otherwise.
fn checkNumericConstraints(
    comptime type_signature: []const u8,
    comptime numeric_signature: []const u8,
    comptime T: type,
    comptime fn_name: []const u8,
    comptime param_name: []const u8,
    comptime checked_domain: *?CheckedDomain, // Is null if called from numeric domain, else only element_type is null
) bool {
    // First, fill in checked_domain
    if (checked_domain.* != null) {
        checked_domain.*.?.element_type = numeric_signature;
    } else {
        if (numeric_signature[numeric_signature.len - 1] == 0) {
            // Remove trailing 0 (implicit null terminator in string literals)
            checked_domain.* = .{ .first_level = .{ .numeric = numeric_signature[0 .. numeric_signature.len - 1] }, .element_type = null };
        } else {
            checked_domain.* = .{ .first_level = .{ .numeric = numeric_signature }, .element_type = null };
        }
    }

    // First check T is numeric
    const is_numeric = types.isNumeric(T);

    comptime var signature: []const u8 = std.mem.trimStart(u8, numeric_signature, " ");
    const actual_type: ?types.NumericType = if (is_numeric) types.numericType(T) else null;

    // Remove trailing 0 (implicit null terminator in string literals)
    signature = std.mem.trimEnd(u8, signature, &.{0});

    // If signature is only one word "numeric", accept any numeric type
    if (7 <= signature.len and
        std.mem.eql(u8, signature, "numeric"))
        return true and is_numeric;

    // Test principal operators
    // "=" means the same as no operator, but can only have one numeric type category, or all of them ("numeric")
    if (std.mem.startsWith(u8, signature, "=")) {
        signature = signature[1..];

        if (7 == signature.len and
            std.mem.eql(u8, signature, "numeric"))
            return true and is_numeric;

        // Exact match
        comptime var matched_type: ?types.NumericType = null;
        inline for (std.meta.fields(types.NumericType)) |num_type| {
            if (num_type.name.len <= signature.len and
                std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
            {
                // Consume matched part
                signature = signature[num_type.name.len..];
                matched_type = @enumFromInt(num_type.value);

                break;
            }
        }

        if (matched_type == null)
            return false and is_numeric;

        return is_numeric and matched_type.? == actual_type.?;
    }

    // Other principal operators: "<", "<=", ">", ">="
    comptime var lower_bound: ?types.NumericType = null;
    comptime var lower_inclusive: bool = false;
    comptime var upper_bound: ?types.NumericType = null;
    comptime var upper_inclusive: bool = false;
    comptime var finish_principal = false;
    while (!finish_principal) {
        signature = std.mem.trimStart(u8, signature, " ");

        if (std.mem.indexOfScalar(u8, signature, '<') == null and
            std.mem.indexOfScalar(u8, signature, '>') == null)
        {
            // Must be at the beginning, so no more principal operators
            finish_principal = true;
            break;
        }

        // If not <= or <, must be >= or >
        const is_lower = std.mem.startsWith(u8, signature, ">");
        signature = signature[1..]; // Consume "<" or ">"
        const is_inclusive = std.mem.startsWith(u8, signature, "=");
        signature = if (is_inclusive) signature[1..] else signature;

        // Match numeric type
        comptime var matched_type: ?types.NumericType = null;
        inline for (std.meta.fields(types.NumericType)) |num_type| {
            if (num_type.name.len <= signature.len and
                std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
            {
                // Consume matched part
                signature = signature[num_type.name.len..];
                matched_type = @enumFromInt(num_type.value);

                break;
            }
        }

        if (matched_type == null)
            return false and is_numeric;

        if (is_lower) {
            if (lower_bound != null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Multiple lower bound operators in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            lower_bound = matched_type;
            lower_inclusive = is_inclusive;
        } else {
            upper_bound = matched_type;
            upper_inclusive = is_inclusive;
        }
    }

    // Now check if actual_type matches the bounds, if they were specified, or the principal type
    comptime var matched_all_bounds = true;
    comptime var checked = 0;
    if (is_numeric) { // Ensure actual_type is not null
        if (lower_bound != null and upper_bound != null) {
            // Both bounds present
            if (lower_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.ge(lower_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.gt(lower_bound.?);

            if (upper_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.le(upper_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.lt(upper_bound.?);

            checked += 2;
        } else if (lower_bound != null) {
            // Only lower bound present
            if (lower_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.ge(lower_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.gt(lower_bound.?);

            checked += 1;
        } else if (upper_bound != null) {
            // Only upper bound present
            if (upper_inclusive)
                matched_all_bounds = matched_all_bounds and actual_type.?.le(upper_bound.?)
            else
                matched_all_bounds = matched_all_bounds and actual_type.?.lt(upper_bound.?);

            checked += 1;
        } // No bounds present -> either we have a standalone principal, or a category as the principal
    }

    // Now check remaining signature for positive/negative types and categories ("@")
    comptime var matched_any_positive = false; // Must match at least one positive type ("" or "|")
    comptime var matched_any_negative = false; // Must not match any negative type ("!")
    comptime var matched_all_categories = true; // All categories must match
    comptime var finished = false;
    while (!finished) : (checked += 1) {
        signature = std.mem.trimStart(u8, signature, " ");
        if (signature.len == 0 or signature[0] == 0) {
            // No more types to check
            finished = true;
            break;
        }

        // Check for "|" and "!" operators. Invalid if checked == 0
        if (signature[0] == '|' or signature[0] == '!') {
            comptime var positive_operator = true;
            if (signature[0] == '|') {
                if (checked == 0) {
                    @compileError(
                        std.fmt.comptimePrint(
                            "Operator '|' cannot be used before any principal numeric type or bound in type signature '{s}' for parameter '{s}' of function '{s}'.",
                            .{ type_signature, param_name, fn_name },
                        ),
                    );
                }

                positive_operator = true;
                signature = signature[1..];
            } else if (signature[0] == '!') {
                if (checked == 0) {
                    @compileError(
                        std.fmt.comptimePrint(
                            "Operator '!' cannot be used before any principal numeric type or bound in type signature '{s}' for parameter '{s}' of function '{s}'.",
                            .{ type_signature, param_name, fn_name },
                        ),
                    );
                }

                positive_operator = false;
                signature = signature[1..];
            }

            signature = std.mem.trimStart(u8, signature, " ");

            // Update checked count
            checked += 1;

            // Match numeric type
            comptime var matched_type: ?types.NumericType = null;
            inline for (std.meta.fields(types.NumericType)) |num_type| {
                if (num_type.name.len <= signature.len and
                    std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
                {
                    // Consume matched part
                    signature = signature[num_type.name.len..];
                    matched_type = @enumFromInt(num_type.value);
                    break;
                }
            }

            if (matched_type == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized numeric type in numeric type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            if (positive_operator) {
                if (is_numeric and matched_type.? == actual_type.?)
                    matched_any_positive = true;
            } else {
                if (is_numeric and matched_type.? == actual_type.?)
                    matched_any_negative = true;
            }

            continue;
        } else if (signature[0] == '@') {
            // Check for "@" category operator. Valid for any checked value
            signature = signature[1..]; // Consume "@"

            // Match numeric category
            comptime var matched_category: ?CategoryOperator = null;
            inline for (std.meta.fields(CategoryOperator)) |num_cat| {
                if (num_cat.name.len <= signature.len and
                    std.mem.eql(u8, signature[0..num_cat.name.len], num_cat.name))
                {
                    // Consume matched part
                    signature = signature[num_cat.name.len..];
                    matched_category = @enumFromInt(num_cat.value);

                    break;
                }
            }

            if (matched_category == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized numeric category in type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            // Check category match
            const category_match = matched_category.?.match()(T);
            matched_all_categories = matched_all_categories and category_match;
            matched_any_positive = matched_any_positive or category_match;

            continue;
        }

        // Check for no operator. Only valid if checked == 0
        if (checked == 0) {
            // Match numeric type
            comptime var matched_type: ?types.NumericType = null;
            inline for (std.meta.fields(types.NumericType)) |num_type| {
                if (num_type.name.len <= signature.len and
                    std.mem.eql(u8, signature[0..num_type.name.len], num_type.name))
                {
                    // Consume matched part
                    signature = signature[num_type.name.len..];
                    matched_type = @enumFromInt(num_type.value);

                    break;
                }
            }

            if (matched_type == null)
                @compileError(
                    std.fmt.comptimePrint(
                        "Unrecognized numeric type in numeric type signature '{s}' for parameter '{s}' of function '{s}'.",
                        .{ type_signature, param_name, fn_name },
                    ),
                );

            if (is_numeric and matched_type.? == actual_type.?)
                matched_any_positive = true;
        } else {
            @compileError(
                std.fmt.comptimePrint(
                    "Missing operator in numeric type signature '{s}' for parameter '{s}' of function '{s}'.",
                    .{ type_signature, param_name, fn_name },
                ),
            );
        }
    }

    return matched_all_bounds and matched_any_positive and !matched_any_negative and matched_all_categories and is_numeric;
}

// 3. DESIGN LANGUAGE AND IMPLEMENT checkParameterTypes FUNCTION

/// Checks if the input types match the provided type rules, triggering
/// compile-time errors if they do not.
///
/// Parameters
/// ----------
/// comptime type_rules (`[]const []const u8`): An array of type rules as strings.
/// Each rule specifies constraints on one or more parameters, using a simple
/// language. Each rule is of one of the following forms:
///
/// - `<selector> is {not} '<type_signature>'` (all parameters selected by `<selector>`
/// must match/not match the type signature)
/// - `if <selector1> is {not} '<type_signature1>' -> <selector2>|other is {not} '<type_signature2>'`
/// (conditional rule: if all parameters in the first selector match/do not match
/// the type signature then all parameters in the second selector/all parameters
/// not in the first selector must match/must not match the second type signature)
/// - `<selector1> <cmp> <selector2>` (all parameters selected by `<selector1>` must
/// compare to all parameters selected by `<selector2>` according to `<cmp>`)
/// - `count(is '<type_signature>') <cmp> N` (where N is an integer, and count
/// counts how many parameters match the type signature)
/// - `<selector> match {.numeric}` (checks if all parameters selected are equal,
/// or their numeric types are equal, if `.numeric` is specified)
///
/// where `<cmp>` is one of: `==`, `!=`, `<`, `<=`, `>`, `>=`, and follows the
/// order in the numeric type hierarchy. `<selector>` can be:
///
/// - `$<name>`, with `<name>` being the name of a parameter, such as `$x`,
/// - `$<name>.numeric`, to refer to the numeric type of a parameter,
/// - `$<name>.domain`, to refer to the domain type of a parameter (numeric, vector, etc.),
/// - `$<name>.child, to refer to the child type of a parameter, i.e., the
/// type a pointer points to (for non-pointer types this is the type itself),
/// - `[$name1, $name2, ...]` for multiple parameters, each also allowing the
/// `.` suffixes as above, or one after the `]`, if not used on individual parameters,
/// to apply to all parameters in the list,
/// - `any` to refer to all parameters, where at least one parameter must comply, also
/// allowing the `.` suffixes as above,
/// - `any of <selector>` to refer to all parameters in the given selector, where
/// at least one parameter must comply. Here <selector> must be a list of parameters, and
/// suffixes can also be applied, but must be in the selector, not on the `any`,
/// - `all` to refer to all parameters, which must all comply, also allowing the `.` suffixes as above,
/// - `all of <selector>` to refer to all parameters in the given selector, which must all comply. Here <selector> must be
/// a list of parameters, and suffixes can also be applied, but must be in the selector, not on the `all`,
/// - `none` to refer to no parameters, i.e., all parameters must not comply, also allowing the `.` suffixes as above, or
/// - `none of <selector>` to refer to no parameters in the given selector,
/// i.e., all parameters in the selector must not comply. Here <selector> must be
/// a list of parameters, and suffixes can also be applied, but must be in the selector, not on the `none`.
///
/// `{not}` is optional and negates the type signature match. `<type_signature>`
/// follows the same format as in `checkParameterType`. All rules must be satisfied
/// for the check to pass, i.e., they are combined with logical and. Some examples
/// of type rules:
///
/// - For a function with two float parameters where at least one must be float:
/// ```zig
/// &.{
///     "[$x, $y] is '<=float'",
///     "any is 'float'", // At least one parameter must be float
/// }
/// ```
///
/// - For a vector scaling function (`v = u * s`, with `u` and `v` vectors and `s` scalar):
/// ```zig
/// &.{
///     "[$u, $v] is 'vector'",
///     "$s is 'numeric'",
/// }
/// ```
///
/// - For a function that requires all parameters to be many-item pointers to numeric
/// types, some needing to be mutable:
/// ```zig
/// &.{
///     "[$a, $b, $c] is '[*]const numeric'",
///     "[$d, $e] is '[*]numeric'",
/// }
/// ```
///
/// - For an in-place function `add(o, a, b)` that adds `a` and `b` into `o`, and
/// we accept only numerics or matrices, but if any parameter is a matrix then all
/// must be matrices:
/// ```zig
/// &.{
///     "[$a, $b] is 'matrix | numeric'",
///     "$o is '* matrix | numeric'",
///     "if any.child is 'matrix' -> all.child is 'matrix'",
/// }
///
/// - For a function `add(x, y)` that adds two parameters and returns the result,
/// but most must be of the same domain:
/// ```zig
/// &.{
///    "all is 'any'", // All must be non-pointer types
///    "all.domain match", // All must be of the same domain
/// }
///
/// comptime args (`anytype`): A struct containing the types of the parameters
/// to check, with field names matching the parameter names. Must be of the form:
///```zig
/// .{
///     .param_name = param_type,
///     ...
/// }
/// ```
///
/// comptime fn_name (`[]const u8`): The name of the function for error messages.
///
/// Returns
/// -------
/// `void`: Triggers compile-time errors if the type rules are not satisfied.
pub fn checkParameterTypes(
    comptime rules: []const []const u8,
    comptime args: anytype,
    comptime fn_name: []const u8,
) void {
    // First, check args is correctly formed and extract parameter names and types
    const arginfo = @typeInfo(@TypeOf(args));
    if (arginfo != .@"struct")
        @compileError(
            std.fmt.comptimePrint(
                "Parameter 'args' to checkParameterTypes must be a struct with parameter names as fields in function '{s}'.",
                .{fn_name},
            ),
        );

    comptime var param_names: [arginfo.@"struct".fields.len][]const u8 = .{undefined} ** arginfo.@"struct".fields.len;
    comptime var param_types: [arginfo.@"struct".fields.len]type = .{undefined} ** arginfo.@"struct".fields.len;
    comptime var i: usize = 0;
    inline for (arginfo.@"struct".fields) |field| {
        param_names[i] = field.name;
        param_types[i] = @field(args, field.name);

        i += 1;
    }

    // Now process each rule
    // comptime var complies_all_rules = true;
    // inline for (rules) |_rule| {
    //     // To allow in-place slicing
    //     comptime var rule = std.mem.trimStart(u8, _rule, " ");

    //     // Parse first part of rule to determine its type, and call the appropriate function
    //     if (std.mem.startsWith(u8, rule, "if ")) {
    //         // Conditional rule
    //         complies_all_rules = complies_all_rules and ifRule(rule, param_names[0..i], param_types[0..i]);
    //     } else if (std.mem.startsWith(u8, rule, "count(")) {
    //         // Count rule
    //         complies_all_rules = complies_all_rules and countRule(rule, param_names[0..i], param_types[0..i]);
    //     }

    //     // Other rules must start with a selector
    // }

    @compileLog(selector(rules[0], &param_names));

    // After calling selector, skip selector.selector_end and then trim spaces
    // Then call checkParameterType with the extracted rule (always between two ' characters)
    // and applying the SelectorSubtype, if not null

    // - `<selector> is {not} '<type_signature>'` (all parameters selected by `<selector>`
    // must match/not match the type signature)
    // - `if <selector1> is {not} '<type_signature1>' -> <selector2>|other is {not} '<type_signature2>'`
    // (conditional rule: if all parameters in the first selector match/do not match
    // the type signature then all parameters in the second selector/all parameters
    // not in the first selector must match/must not match the second type signature)
    // - `<selector1> <cmp> <selector2>` (all parameters selected by `<selector1>` must
    // compare to all parameters selected by `<selector2>` according to `<cmp>`)
    // - `count(is '<type_signature>') <cmp> N` (where N is an integer, and count
    // counts how many parameters match the type signature)
    // - `<selector> match {.numeric}` (checks if all parameters selected are equal,
    // or their numeric types are equal, if `.numeric` is specified)

    @compileError("Debug");
}

const SelectorType = enum {
    inclusive,
    lean,
    exclusive,
};

const SelectorSubtype = enum {
    numeric,
    domain,
    child,

    pub fn match(self: SelectorSubtype) fn (type) type {
        return switch (self) {
            .numeric => types.Numeric,
            .domain => types.Domain, // problem: is the enum itself, no function that returns a type
            .child => types.Child, // problem: only works for pointers
        };
    }
};

/// Gets a selector and parses it to a usable form.
fn selector(
    comptime selector_: []const u8,
    comptime param_names: []const []const u8,
) struct {
    selected_indices: [max_parameters]usize,
    selected_subtypes: [max_parameters]?SelectorSubtype,
    selected_count: usize,
    selector_type: SelectorType,
    selector_end: usize,
} {
    comptime var sel = std.mem.trimStart(u8, selector_, " ");

    comptime var selected_indices: [max_parameters]usize = .{@as(u32, 0) -% 1} ** max_parameters;
    comptime var selected_subtypes: [max_parameters]?SelectorSubtype = .{null} ** max_parameters;
    comptime var selected_count: usize = 0;
    comptime var selector_type: SelectorType = .inclusive;
    comptime var selector_end: usize = 0;

    // Process selector
    switch (sel[0]) {
        '$' => {
            // Single parameter
            const first_dot_index = std.mem.indexOfScalar(u8, sel, '.');
            const first_space_index = std.mem.indexOfScalar(u8, sel, ' ');
            const param_name_end = @min(
                @min(
                    if (first_dot_index != null) first_dot_index.? else sel.len,
                    if (first_space_index != null) first_space_index.? else sel.len,
                ),
                sel.len,
            );

            const dot_suffix = first_dot_index != null and first_dot_index.? <= param_name_end;

            const param_name = sel[1..param_name_end];
            var found = false;
            inline for (param_names, 0..) |pname, index| {
                if (std.mem.eql(u8, pname, param_name)) {
                    selected_indices[selected_count] = index;
                    selected_count += 1;
                    selector_end = param_name_end;
                    found = true;

                    break;
                }
            }

            if (!found)
                @compileError(
                    std.fmt.comptimePrint(
                        "Parameter '{s}' in selector not found among function parameters.",
                        .{param_name},
                    ),
                );

            // Check for suffix
            if (dot_suffix) {
                comptime var suffix_end = std.mem.indexOfScalar(u8, sel[selector_end..], ' ');
                if (suffix_end == null)
                    suffix_end = sel.len - selector_end;

                const suffix = sel[selector_end + 1 .. selector_end + suffix_end.?];
                comptime var subtype: ?SelectorSubtype = null;
                inline for (std.meta.fields(SelectorSubtype)) |stype| {
                    if (stype.name.len == suffix.len and
                        std.mem.eql(u8, suffix, stype.name))
                    {
                        subtype = @enumFromInt(stype.value);

                        break;
                    }
                }

                if (subtype == null)
                    @compileError(
                        std.fmt.comptimePrint(
                            "Unrecognized selector subtype '.{s}' in selector for parameter '{s}'.",
                            .{ suffix, param_name },
                        ),
                    );

                selected_subtypes[selected_count - 1] = subtype;
                selector_end += suffix_end.?;
            }

            selector_type = .inclusive;
        },
        '[' => {
            // Multiple parameters
            const closing_bracket_index = std.mem.indexOfScalar(u8, sel, ']');
            if (closing_bracket_index == null)
                @compileError("Unmatched '[' in selector.");

            comptime var params_list: []const u8 = sel[1 .. closing_bracket_index.? + 1]; // Include closing ']'
            comptime var any_has_suffix = false;
            comptime var reached_end = false;
            while (params_list.len > 0 and !reached_end) {
                // Skip to next parameter (starts with '$')
                const dollar_index = std.mem.indexOfScalar(u8, params_list, '$');
                if (dollar_index == null)
                    break; // No more parameters

                params_list = params_list[dollar_index.? + 1 ..];

                // Find end of parameter name (space, comma, dot, or end of list)
                const first_dot_index = std.mem.indexOfScalar(u8, params_list, '.');
                comptime var first_space_index = std.mem.indexOfScalar(u8, params_list, ' ');
                comptime var first_comma_index = std.mem.indexOfScalar(u8, params_list, ',');
                comptime var closing_bracket_index_inner = std.mem.indexOfScalar(u8, params_list, ']').?; // Must exist

                const param_name_end: comptime_int = @min(
                    @min(
                        @min(
                            if (first_dot_index != null) first_dot_index.? else closing_bracket_index_inner,
                            if (first_space_index != null) first_space_index.? else closing_bracket_index_inner,
                        ),
                        if (first_comma_index != null) first_comma_index.? else closing_bracket_index_inner,
                    ),
                    closing_bracket_index_inner,
                );

                reached_end = param_name_end == closing_bracket_index_inner;
                const dot_suffix = first_dot_index != null and first_dot_index.? <= param_name_end;

                const param_name = params_list[0..param_name_end];
                var found = false;
                inline for (param_names, 0..) |pname, index| {
                    if (std.mem.eql(u8, pname, param_name)) {
                        if (std.mem.indexOfScalar(usize, selected_indices[0..selected_count], index) != null)
                            @compileError(
                                std.fmt.comptimePrint(
                                    "Parameter '{s}' in selector is duplicated.",
                                    .{param_name},
                                ),
                            );

                        selected_indices[selected_count] = index;
                        selected_count += 1;
                        found = true;

                        break;
                    }
                }

                if (!found)
                    @compileError(
                        std.fmt.comptimePrint(
                            "Parameter '{s}' in selector not found among function parameters.",
                            .{param_name},
                        ),
                    );

                if (reached_end)
                    break;

                // Check for suffix
                if (dot_suffix) {
                    any_has_suffix = true;

                    first_space_index = std.mem.indexOfScalar(u8, params_list[param_name_end..], ' ');
                    first_comma_index = std.mem.indexOfScalar(u8, params_list[param_name_end..], ',');
                    closing_bracket_index_inner = std.mem.indexOfScalar(u8, params_list[param_name_end..], ']').?; // Must exist
                    const suffix_end = @min(
                        @min(
                            if (first_space_index != null) first_space_index.? else closing_bracket_index_inner - param_name_end,
                            if (first_comma_index != null) first_comma_index.? else closing_bracket_index_inner - param_name_end,
                        ),
                        closing_bracket_index_inner - param_name_end,
                    );
                    reached_end = suffix_end == closing_bracket_index_inner - param_name_end;
                    const suffix = params_list[param_name_end + 1 .. param_name_end + suffix_end];
                    comptime var subtype: ?SelectorSubtype = null;
                    inline for (std.meta.fields(SelectorSubtype)) |stype| {
                        if (stype.name.len == suffix.len and
                            std.mem.eql(u8, suffix, stype.name))
                        {
                            subtype = @enumFromInt(stype.value);

                            break;
                        }
                    }

                    if (subtype == null)
                        @compileError(
                            std.fmt.comptimePrint(
                                "Unrecognized selector subtype '.{s}' in selector for parameter '{s}'.",
                                .{ suffix, param_name },
                            ),
                        );

                    selected_subtypes[selected_count - 1] = subtype;
                }
            }

            selector_end = closing_bracket_index.? + 1;

            // Check if list has a suffix
            comptime var suffix_start = std.mem.indexOfScalar(u8, sel[selector_end..], '.');
            if (suffix_start != null and suffix_start.? == 0) { // . must be immediately after ]
                if (any_has_suffix)
                    @compileError("Cannot have both per-parameter suffixes and a list-wide suffix in selector.");

                suffix_start = selector_end + suffix_start.?;
                comptime var suffix_end = std.mem.indexOfScalar(u8, sel[suffix_start.?..], ' ');
                if (suffix_end == null)
                    suffix_end = sel.len - suffix_start;

                const suffix = sel[suffix_start.? + 1 .. suffix_start.? + suffix_end.?];
                comptime var subtype: ?SelectorSubtype = null;
                inline for (std.meta.fields(SelectorSubtype)) |stype| {
                    if (stype.name.len == suffix.len and
                        std.mem.eql(u8, suffix, stype.name))
                    {
                        subtype = @enumFromInt(stype.value);

                        break;
                    }
                }

                if (subtype == null)
                    @compileError(
                        std.fmt.comptimePrint(
                            "Unrecognized selector subtype '.{s}' in selector for parameter list.",
                            .{suffix},
                        ),
                    );

                inline for (selected_subtypes[0..selected_count]) |*subtype_ptr| {
                    subtype_ptr.* = subtype;
                }

                selector_end += suffix_end.?;
            }

            selector_type = .inclusive;
        },
        'a' => switch (sel[1]) {
            'l' => {
                // "all" or "all of"
                if (std.mem.startsWith(u8, sel, "all of ")) {
                    // "all of <selector>"
                    comptime var inner_selector: []const u8 = sel[7..];
                    inner_selector = std.mem.trimStart(u8, inner_selector, " ");
                    const inner_result = selector(inner_selector, param_names);

                    inline for (inner_result.selected_indices[0..inner_result.selected_count], 0..) |idx, index| {
                        selected_indices[index] = idx;
                    }
                    selected_subtypes = inner_result.selected_subtypes;
                    selected_count = inner_result.selected_count;
                    selector_end = inner_result.selector_end + 7;

                    selector_type = .inclusive;
                } else if (std.mem.startsWith(u8, sel, "all")) {
                    // "all"
                    inline for (param_names, 0..) |_, index| {
                        selected_indices[index] = index;
                    }
                    // selected_subtypes remain null
                    selected_count = param_names.len;
                    selector_end = 3;

                    // Check for suffix
                    comptime var suffix_start = std.mem.indexOfScalar(u8, sel[selector_end..], '.');
                    if (suffix_start != null and suffix_start.? == 0) {
                        suffix_start = selector_end + suffix_start.?;
                        comptime var suffix_end = std.mem.indexOfScalar(u8, sel[suffix_start.?..], ' ');
                        if (suffix_end == null)
                            suffix_end = sel.len - suffix_start;

                        const suffix = sel[suffix_start.? + 1 .. suffix_start.? + suffix_end.?];
                        comptime var subtype: ?SelectorSubtype = null;
                        inline for (std.meta.fields(SelectorSubtype)) |stype| {
                            if (stype.name.len == suffix.len and
                                std.mem.eql(u8, suffix, stype.name))
                            {
                                subtype = @enumFromInt(stype.value);

                                break;
                            }
                        }

                        if (subtype == null)
                            @compileError(
                                std.fmt.comptimePrint(
                                    "Unrecognized selector subtype '.{s}' in selector for 'all'.",
                                    .{suffix},
                                ),
                            );

                        inline for (selected_subtypes[0..selected_count]) |*subtype_ptr| {
                            subtype_ptr.* = subtype;
                        }

                        selector_end += suffix_end.?;
                    }

                    selector_type = .inclusive;
                } else unreachable;
            },
            'n' => {
                // "any" or "any of"
                if (std.mem.startsWith(u8, sel, "any of ")) {
                    // "any of <selector>"
                    comptime var inner_selector: []const u8 = sel[7..];
                    inner_selector = std.mem.trimStart(u8, inner_selector, " ");
                    const inner_result = selector(inner_selector, param_names);

                    inline for (inner_result.selected_indices[0..inner_result.selected_count], 0..) |idx, index| {
                        selected_indices[index] = idx;
                    }
                    selected_subtypes = inner_result.selected_subtypes;
                    selected_count = inner_result.selected_count;
                    selector_end = inner_result.selector_end + 7;

                    selector_type = .lean;
                } else if (std.mem.startsWith(u8, sel, "any")) {
                    // "any"
                    inline for (param_names, 0..) |_, index| {
                        selected_indices[index] = index;
                    }
                    // selected_subtypes remain null
                    selected_count = param_names.len;
                    selector_end = 3;

                    // Check for suffix
                    comptime var suffix_start = std.mem.indexOfScalar(u8, sel[selector_end..], '.');
                    if (suffix_start != null and suffix_start.? == 0) {
                        suffix_start = selector_end + suffix_start.?;
                        comptime var suffix_end = std.mem.indexOfScalar(u8, sel[suffix_start.?..], ' ');
                        if (suffix_end == null)
                            suffix_end = sel.len - suffix_start;

                        const suffix = sel[suffix_start.? + 1 .. suffix_start.? + suffix_end.?];
                        comptime var subtype: ?SelectorSubtype = null;
                        inline for (std.meta.fields(SelectorSubtype)) |stype| {
                            if (stype.name.len == suffix.len and
                                std.mem.eql(u8, suffix, stype.name))
                            {
                                subtype = @enumFromInt(stype.value);

                                break;
                            }
                        }

                        if (subtype == null)
                            @compileError(
                                std.fmt.comptimePrint(
                                    "Unrecognized selector subtype '.{s}' in selector for 'any'.",
                                    .{suffix},
                                ),
                            );

                        inline for (selected_subtypes[0..selected_count]) |*subtype_ptr| {
                            subtype_ptr.* = subtype;
                        }

                        selector_end += suffix_end.?;
                    }

                    selector_type = .inclusive;
                } else unreachable;
            },
            else => unreachable,
        },
        'n' => {
            // "none" or "none of"
            if (std.mem.startsWith(u8, sel, "none of ")) {
                // "none of <selector>"
                comptime var inner_selector: []const u8 = sel[8..];
                inner_selector = std.mem.trimStart(u8, inner_selector, " ");
                const inner_result = selector(inner_selector, param_names);

                inline for (inner_result.selected_indices[0..inner_result.selected_count], 0..) |idx, index| {
                    selected_indices[index] = idx;
                }
                selected_subtypes = inner_result.selected_subtypes;
                selected_count = inner_result.selected_count;
                selector_end = inner_result.selector_end + 8;

                selector_type = .exclusive;
            } else if (std.mem.startsWith(u8, sel, "none")) {
                // "none"
                inline for (param_names, 0..) |_, index| {
                    selected_indices[index] = index;
                }
                // selected_subtypes remain null
                selected_count = param_names.len;
                selector_end = 4;

                // Check for suffix
                comptime var suffix_start = std.mem.indexOfScalar(u8, sel[selector_end..], '.');
                if (suffix_start != null and suffix_start.? == 0) {
                    suffix_start = selector_end + suffix_start.?;
                    comptime var suffix_end = std.mem.indexOfScalar(u8, sel[suffix_start.?..], ' ');
                    if (suffix_end == null)
                        suffix_end = sel.len - suffix_start;

                    const suffix = sel[suffix_start.? + 1 .. suffix_start.? + suffix_end.?];
                    comptime var subtype: ?SelectorSubtype = null;
                    inline for (std.meta.fields(SelectorSubtype)) |stype| {
                        if (stype.name.len == suffix.len and
                            std.mem.eql(u8, suffix, stype.name))
                        {
                            subtype = @enumFromInt(stype.value);

                            break;
                        }
                    }

                    if (subtype == null)
                        @compileError(
                            std.fmt.comptimePrint(
                                "Unrecognized selector subtype '.{s}' in selector for 'none'.",
                                .{suffix},
                            ),
                        );

                    inline for (selected_subtypes[0..selected_count]) |*subtype_ptr| {
                        subtype_ptr.* = subtype;
                    }

                    selector_end += suffix_end.?;
                }

                selector_type = .exclusive;
            } else unreachable;
        },
        else => unreachable,
    }

    return .{
        .selected_indices = selected_indices,
        .selected_subtypes = selected_subtypes,
        .selected_count = selected_count,
        .selector_type = selector_type,
        .selector_end = selector_end,
    };
}

// fn ifRule(
//     comptime rule: []const u8,
//     comptime param_names: []const []const u8,
//     comptime param_types: []const type,
// ) bool {
//     return true;
// }
