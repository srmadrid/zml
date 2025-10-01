const std = @import("std");

const types = @import("../../types.zig");
const ReturnType1 = types.ReturnType1;
const ReturnType2 = types.ReturnType2;
const EnsureMatrix = types.EnsureMatrix;
const Coerce = types.Coerce;
const Numeric = types.Numeric;
const Order = types.Order;
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

    if (comptime !types.isSymmetricDenseMatrix(@TypeOf(x))) {
        var result: R = try .init(allocator, y.size);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(y)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < y.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < y.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j + i * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j + i * y.ld], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i + j * y.ld]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i + j * y.ld], ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < y.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(y))) {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[j * y.ld + i]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[j * y.ld + i], ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < y.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x, y.data[i * y.ld + j]);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x, y.data[i * y.ld + j], ctx);
                                }
                            }
                        }
                    }
                }
            }
        }

        return result;
    } else if (comptime !types.isSymmetricDenseMatrix(@TypeOf(y))) {
        var result: R = try .init(allocator, x.size);
        errdefer result.deinit(allocator);

        const opinfo = @typeInfo(@TypeOf(op));
        if (comptime types.orderOf(@TypeOf(result)) == .col_major) {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = 0;
                            while (i <= j) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    } else {
                        var j: u32 = 0;
                        while (j < x.size) : (j += 1) {
                            var i: u32 = j;
                            while (i < x.size) : (i += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y, ctx);
                                }
                            }
                        }
                    }
                }
            }
        } else {
            if (comptime types.orderOf(@TypeOf(x)) == .col_major) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y, ctx);
                                }
                            }
                        }
                    }
                }
            } else {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = i;
                            while (j < x.size) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y, ctx);
                                }
                            }
                        }
                    } else {
                        var i: u32 = 0;
                        while (i < x.size) : (i += 1) {
                            var j: u32 = 0;
                            while (j <= i) : (j += 1) {
                                if (comptime opinfo.@"fn".params.len == 2) {
                                    result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y);
                                } else if (comptime opinfo.@"fn".params.len == 3) {
                                    result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y, ctx);
                                }
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
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cu cu
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
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
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
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cu ru
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
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu cl rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cu ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cl cu rl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl cl ru
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
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
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu ru cu
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
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cu rl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl ru cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // cl ru cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl rl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
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
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu ru ru
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
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cu rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cu rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j + i * result.ld] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j + i * result.ld] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // cl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // cl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i + j * result.ld] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i + j * result.ld] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
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
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cu cu
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
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru cl cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cu cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = 0;
                                while (i <= j) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // rl cu cl
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cl cu
                            var j: u32 = 0;
                            while (j < y.size) : (j += 1) {
                                var i: u32 = j;
                                while (i < y.size) : (i += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j + i * y.ld], ctx);
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
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cu ru
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
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru cl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cu ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i + j * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i + j * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // rl cu rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j + i * x.ld], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j + i * x.ld], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl cl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i + j * x.ld], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i + j * x.ld], y.data[j * y.ld + i], ctx);
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
            if (comptime types.orderOf(@TypeOf(y))) {
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru ru cu
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
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // ru rl cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl ru cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        } else { // rl ru cl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i + j * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i + j * y.ld], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl rl cu
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j + i * y.ld]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j + i * y.ld], ctx);
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
                if (comptime types.uploOf(@TypeOf(result))) {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru ru ru
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
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // ru rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // ru rl rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    if (comptime types.uploOf(@TypeOf(x))) {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl ru ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = i;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[j * result.ld + i] = op(x.data[i * x.ld + j], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[j * result.ld + i] = try op(x.data[i * x.ld + j], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        } else { // rl ru rl
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[j * x.ld + i], y.data[i * y.ld + j]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[j * x.ld + i], y.data[i * y.ld + j], ctx);
                                    }
                                }
                            }
                        }
                    } else {
                        if (comptime types.uploOf(@TypeOf(y))) { // rl rl ru
                            var i: u32 = 0;
                            while (i < x.size) : (i += 1) {
                                var j: u32 = 0;
                                while (j < x.size) : (j += 1) {
                                    if (comptime opinfo.@"fn".params.len == 2) {
                                        result.data[i * result.ld + j] = op(x.data[i * x.ld + j], y.data[j * y.ld + i]);
                                    } else if (comptime opinfo.@"fn".params.len == 3) {
                                        result.data[i * result.ld + j] = try op(x.data[i * x.ld + j], y.data[j * y.ld + i], ctx);
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
