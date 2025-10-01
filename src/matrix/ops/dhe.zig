const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
const Uplo = types.Uplo;
const ops = @import("../ops.zig");
const constants = @import("../../constants.zig");
const int = @import("../../int.zig");

const matrix = @import("../../matrix.zig");

pub fn apply2(
    allocator: std.mem.Allocator,
    x: anytype,
    y: anytype,
    comptime op: anytype,
    ctx: anytype,
) !EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, Numeric(@TypeOf(x)), Numeric(@TypeOf(y)))) {
    const X: type = Numeric(@TypeOf(x));
    const Y: type = Numeric(@TypeOf(y));
    const R: type = EnsureMatrix(Coerce(@TypeOf(x), @TypeOf(y)), ReturnType2(op, X, Y));

    if (comptime !types.isHermitianDenseMatrix(@TypeOf(x))) {
        const ruplo: Uplo = comptime types.uploOf(@TypeOf(y));
        var result: R = if (comptime types.isHermitianDenseMatrix(R))
            try .init( // Hermitian
                allocator,
                y.size,
            )
        else
            try .init( // General
                allocator,
                y.size,
                y.size,
            );
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cu
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) { // Result is a general matrix
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }
                        }
                    } else { // cu cl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cu
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl cl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j + j * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j + j * y.ld], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu ru
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }
                        }
                    } else { // cu rl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl ru
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl rl
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x, y.data[j * y.ld + j]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x, y.data[j * y.ld + j], ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cu
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru cl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cu
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }
                        }
                    } else { // rl cl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i + j * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i + i * y.ld]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i + i * y.ld], ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru ru
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru rl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl ru
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        x,
                                        ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x, y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x, y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }
                        }
                    } else { // rl rl
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            x,
                                            ops.conjugate(y.data[i * y.ld + j], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x, y.data[i * y.ld + i]);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x, y.data[i * y.ld + i], ctx);
                            }
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isHermitianDenseMatrix(@TypeOf(y))) {
        const ruplo: Uplo = comptime types.uploOf(@TypeOf(x));
        var result: R = if (comptime types.isHermitianDenseMatrix(R))
            try .init( // Hermitian
                allocator,
                x.size,
            )
        else
            try .init( // General
                allocator,
                x.size,
                x.size,
            );
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu cu
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) { // Result is a general matrix
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }
                        }
                    } else { // cu cl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cl cu
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl cl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j + j * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j + j * x.ld], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cu ru
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }
                        }
                    } else { // cu rl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i < j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // cl ru
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }
                        }
                    } else { // cl rl
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[j + j * result.ld] = op(x.data[j * x.ld + j], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[j + j * result.ld] = try op(x.data[j * x.ld + j], y, ctx);
                            }

                            var i: u32 = j + 1;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru cu
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru cl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // rl cu
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j + i * x.ld], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j + i * x.ld], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }
                        }
                    } else { // rl cl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i + j * x.ld], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i + i * x.ld], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i + i * x.ld], y, ctx);
                            }
                        }
                    }
                }
            } else {
                if (comptime ruplo == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // ru ru
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else { // ru rl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }

                            var j: u32 = i + 1;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) { // rl ru
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                    );
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(
                                        ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                        y,
                                        ctx,
                                    );
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[j * x.ld + i], y);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[j * x.ld + i], y, ctx);
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }
                        }
                    } else { // rl rl
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j < i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }

                                if (comptime types.isComplex(X)) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(
                                            ops.conjugate(x.data[i * x.ld + j], .{}) catch unreachable,
                                            y,
                                            ctx,
                                        );
                                    }
                                }
                            }

                            if (comptime opinfo.@"fn".params.len == 2) {
                                result.data[i * result.ld + i] = op(x.data[i * x.ld + i], y);
                            } else if (comptime opinfo.@"fn".params.len == 3) {
                                result.data[i * result.ld + i] = try op(x.data[i * x.ld + i], y, ctx);
                            }
                        }
                    }
                }
            }
        }

        return result;
    }

    var result: R = try .init(allocator, x.size);
    errdefer result.deinit(allocator);

    const opinfo = @typeInfo(@TypeOf(op));
    if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cu rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cu rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // cl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // cl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // cl rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i + j * x.ld], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j + i * x.ld], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i + j * x.ld],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i + j * y.ld]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i + j * y.ld],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j + i * y.ld], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result)) == .upper) {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // ru rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // ru rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x)) == .upper) {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = ops.conjugate(op(x.data[i * x.ld + j], y.data[i * y.ld + j]), .{}) catch unreachable;
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = ops.conjugate(try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx), .{}) catch unreachable;
                                    }
                                }
                            }
                        } else { // rl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            ops.conjugate(x.data[j * x.ld + i], .{}) catch unreachable,
                                            y.data[i * y.ld + j],
                                            ctx,
                                        );
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y)) == .upper) { // rl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                        );
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(
                                            x.data[i * x.ld + j],
                                            ops.conjugate(y.data[j * y.ld + i], .{}) catch unreachable,
                                            ctx,
                                        );
                                    }
                                }
                            }
                        } else { // rl rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    return result;
}
