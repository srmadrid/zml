const std = @import("std");
const core = @import("../core.zig");
const ci = @import("../c.zig");
const options = @import("options");
const Complex = std.math.Complex;
const ndarray = @import("../ndarray/ndarray.zig");

const Scalar = core.types.Scalar;

pub const Order = enum(c_uint) {
    RowMajor = 101,
    ColumnMajor = 102,
};

pub const Transpose = enum(c_uint) {
    NoTrans = 111,
    Trans = 112,
    ConjTrans = 113,
    ConjNoTrans = 114,
};

pub const Uplo = enum(c_uint) {
    Upper = 121,
    Lower = 122,
};

pub const Diag = enum(c_uint) {
    NonUnit = 131,
    Unit = 132,
};

pub const Side = enum(c_uint) {
    Left = 141,
    Right = 142,
};

// Level 1 BLAS
pub inline fn asum(comptime T: type, n: isize, x: [*]const T, incx: isize) Scalar(T) {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sasum(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dasum(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_scasum(@intCast(n), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_dzasum(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/asum.zig").asum(T, n, x, incx);
}
pub fn sasum(n: isize, x: [*]const f32, incx: isize) f32 {
    return asum(f32, n, x, incx);
}
pub fn dasum(n: isize, x: [*]const f64, incx: isize) f64 {
    return asum(f64, n, x, incx);
}
pub fn scasum(n: isize, x: [*]const Complex(f32), incx: isize) f32 {
    return asum(Complex(f32), n, x, incx);
}
pub fn dzasum(n: isize, x: [*]const Complex(f64), incx: isize) f64 {
    return asum(Complex(f64), n, x, incx);
}

pub inline fn axpy(comptime T: type, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_saxpy(@intCast(n), alpha, x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_daxpy(@intCast(n), alpha, x, @intCast(incx), y, @intCast(incy));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_caxpy(@intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zaxpy(@intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/axpy.zig").axpy(T, n, alpha, x, incx, y, incy);
}
pub fn saxpy(n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
    return axpy(f32, n, alpha, x, incx, y, incy);
}
pub fn daxpy(n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
    return axpy(f64, n, alpha, x, incx, y, incy);
}
pub fn caxpy(n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize) void {
    return axpy(Complex(f32), n, alpha, x, incx, y, incy);
}
pub fn zaxpy(n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize) void {
    return axpy(Complex(f64), n, alpha, x, incx, y, incy);
}

pub inline fn copy(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_scopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dcopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ccopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zcopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/copy.zig").copy(T, n, x, incx, y, incy);
}
pub fn scopy(n: isize, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
    return copy(f32, n, x, incx, y, incy);
}
pub fn dcopy(n: isize, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
    return copy(f64, n, x, incx, y, incy);
}
pub fn ccopy(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize) void {
    return copy(Complex(f32), n, x, incx, y, incy);
}
pub fn zcopy(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize) void {
    return copy(Complex(f64), n, x, incx, y, incy);
}

pub inline fn dot(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) Scalar(T) {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sdot(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_ddot(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/dot.zig").dot(T, n, x, incx, y, incy);
}
pub fn sdot(n: isize, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize) f32 {
    return dot(f32, n, x, incx, y, incy);
}
pub fn ddot(n: isize, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize) f64 {
    return dot(f64, n, x, incx, y, incy);
}

pub inline fn dotc(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    //const supported = core.types.numericType(T);

    //if (options.link_cblas != null) {
    //    switch (supported) {
    //        .Complex => {
    //            if (Scalar(T) == f32) {
    //                return ci.cblas_cdotc(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            } else if (Scalar(T) == f64) {
    //                return ci.cblas_zdotc(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            }
    //        },
    //        else => {},
    //    }
    //}

    return @import("blas/dotc.zig").dotc(T, n, x, incx, y, incy);
}
pub fn cdotc(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize) Complex(f32) {
    return dotc(Complex(f32), n, x, incx, y, incy);
}
pub fn zdotc(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize) Complex(f64) {
    return dotc(Complex(f64), n, x, incx, y, incy);
}

pub inline fn dotc_sub(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize, ret: *T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cdotc_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zdotc_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                }
            },
            else => {},
        }
    }

    return @import("blas/dotc_sub.zig").dotc_sub(T, n, x, incx, y, incy, ret);
}
pub fn cdotc_sub(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, ret: *Complex(f32)) void {
    return dotc_sub(Complex(f32), n, x, incx, y, incy, ret);
}
pub fn zdotc_sub(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, ret: *Complex(f64)) void {
    return dotc_sub(Complex(f64), n, x, incx, y, incy, ret);
}

pub inline fn dotu(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    //const supported = core.types.numericType(T);

    //if (options.link_cblas != null) {
    //    switch (supported) {
    //        .Complex => {
    //            if (Scalar(T) == f32) {
    //                return ci.cblas_cdotu(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            } else if (Scalar(T) == f64) {
    //                return ci.cblas_zdotu(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            }
    //        },
    //        else => {},
    //    }
    //}

    return @import("blas/dotu.zig").dotu(T, n, x, incx, y, incy);
}
pub fn cdotu(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize) Complex(f32) {
    return dotu(Complex(f32), n, x, incx, y, incy);
}
pub fn zdotu(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize) Complex(f64) {
    return dotu(Complex(f64), n, x, incx, y, incy);
}

pub inline fn dotu_sub(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize, ret: *T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cdotu_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zdotu_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                }
            },
            else => {},
        }
    }

    return @import("blas/dotu_sub.zig").dotu_sub(T, n, x, incx, y, incy, ret);
}
pub fn cdotu_sub(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, ret: *Complex(f32)) void {
    return dotu_sub(Complex(f32), n, x, incx, y, incy, ret);
}
pub fn zdotu_sub(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, ret: *Complex(f64)) void {
    return dotu_sub(Complex(f64), n, x, incx, y, incy, ret);
}

pub inline fn nrm2(comptime T: type, n: isize, x: [*]const T, incx: isize) Scalar(T) {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_snrm2(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dnrm2(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_scnrm2(@intCast(n), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_dznrm2(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/nrm2.zig").nrm2(T, n, x, incx);
}
pub fn snrm2(n: isize, x: [*]const f32, incx: isize) f32 {
    return nrm2(f32, n, x, incx);
}
pub fn dnrm2(n: isize, x: [*]const f64, incx: isize) f64 {
    return nrm2(f64, n, x, incx);
}
pub fn scnrm2(n: isize, x: [*]const Complex(f32), incx: isize) f32 {
    return nrm2(Complex(f32), n, x, incx);
}
pub fn dznrm2(n: isize, x: [*]const Complex(f64), incx: isize) f64 {
    return nrm2(Complex(f64), n, x, incx);
}

pub inline fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: Scalar(T), s: Scalar(T)) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_srot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                } else if (T == f64) {
                    return ci.cblas_drot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csrot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zdrot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                }
            },
            else => {},
        }
    }

    return @import("blas/rot.zig").rot(T, n, x, incx, y, incy, c, s);
}
pub fn srot(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, c: f32, s: f32) void {
    return rot(f32, n, x, incx, y, incy, c, s);
}
pub fn drot(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, c: f64, s: f64) void {
    return rot(f64, n, x, incx, y, incy, c, s);
}
pub fn csrot(n: isize, x: [*]Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize, c: f32, s: f32) void {
    return rot(Complex(f32), n, x, incx, y, incy, c, s);
}
pub fn zdrot(n: isize, x: [*]Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize, c: f64, s: f64) void {
    return rot(Complex(f64), n, x, incx, y, incy, c, s);
}

pub inline fn rotg(comptime T: type, a: *T, b: *T, c: *Scalar(T), s: *T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_srotg(a, b, c, s);
                } else if (T == f64) {
                    return ci.cblas_drotg(a, b, c, s);
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_crotg(a, b, c, s);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zrotg(a, b, c, s);
                }
            },
            else => {},
        }
    }

    return @import("blas/rotg.zig").rotg(T, a, b, c, s);
}
pub fn srotg(a: *f32, b: *f32, c: *f32, s: *f32) void {
    return rotg(f32, a, b, c, s);
}
pub fn drotg(a: *f64, b: *f64, c: *f64, s: *f64) void {
    return rotg(f64, a, b, c, s);
}
pub fn crotg(a: *Complex(f32), b: *Complex(f32), c: *f32, s: *Complex(f32)) void {
    return rotg(Complex(f32), a, b, c, s);
}
pub fn zrotg(a: *Complex(f64), b: *Complex(f64), c: *f64, s: *Complex(f64)) void {
    return rotg(Complex(f64), a, b, c, s);
}

pub inline fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: [*]const T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_srotm(@intCast(n), x, @intCast(incx), y, @intCast(incy), param);
                } else if (T == f64) {
                    return ci.cblas_drotm(@intCast(n), x, @intCast(incx), y, @intCast(incy), param);
                }
            },
            else => {},
        }
    }

    return @import("blas/rotm.zig").rotm(T, n, x, incx, y, incy, param);
}
pub fn srotm(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, param: [*]const f32) void {
    return rotm(f32, n, x, incx, y, incy, param);
}
pub fn drotm(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, param: [*]const f64) void {
    return rotm(f64, n, x, incx, y, incy, param);
}

pub inline fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: [*]T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_srotmg(d1, d2, x1, y1, param);
                } else if (T == f64) {
                    return ci.cblas_drotmg(d1, d2, x1, y1, param);
                }
            },
            else => {},
        }
    }

    return @import("blas/rotmg.zig").rotmg(T, d1, d2, x1, y1, param);
}
pub fn srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*]f32) void {
    return rotmg(f32, d1, d2, x1, y1, param);
}
pub fn drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*]f64) void {
    return rotmg(f64, d1, d2, x1, y1, param);
}

pub inline fn scal(comptime T: type, n: isize, a: T, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sscal(@intCast(n), a, x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dscal(@intCast(n), a, x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cscal(@intCast(n), &a, x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zscal(@intCast(n), &a, x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/scal.zig").scal(T, n, a, x, incx);
}
pub fn sscal(n: isize, alpha: f32, x: [*]f32, incx: isize) void {
    return scal(f32, n, alpha, x, incx);
}
pub fn dscal(n: isize, alpha: f64, x: [*]f64, incx: isize) void {
    return scal(f64, n, alpha, x, incx);
}
pub fn cscal(n: isize, alpha: Complex(f32), x: [*]Complex(f32), incx: isize) void {
    return scal(Complex(f32), n, alpha, x, incx);
}
pub fn zscal(n: isize, alpha: Complex(f64), x: [*]Complex(f64), incx: isize) void {
    return scal(Complex(f64), n, alpha, x, incx);
}
pub fn csscal(n: isize, alpha: f32, x: [*]Complex(f32), incx: isize) void {
    return scal(Complex(f32), n, Complex(f32){ .re = alpha, .im = 0.0 }, x, incx);
}
pub fn zdscal(n: isize, alpha: f64, x: [*]Complex(f64), incx: isize) void {
    return scal(Complex(f64), n, Complex(f64){ .re = alpha, .im = 0.0 }, x, incx);
}

pub inline fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/swap.zig").swap(T, n, x, incx, y, incy);
}
pub fn sswap(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize) void {
    return swap(f32, n, x, incx, y, incy);
}
pub fn dswap(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize) void {
    return swap(f64, n, x, incx, y, incy);
}
pub fn cswap(n: isize, x: [*]Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize) void {
    return swap(Complex(f32), n, x, incx, y, incy);
}
pub fn zswap(n: isize, x: [*]Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize) void {
    return swap(Complex(f64), n, x, incx, y, incy);
}

pub inline fn iamax(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_isamax(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_idamax(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_icamax(@intCast(n), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_izamax(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/iamax.zig").iamax(T, n, x, incx);
}
pub fn isamax(n: isize, x: [*]const f32, incx: isize) usize {
    return iamax(f32, n, x, incx);
}
pub fn idamax(n: isize, x: [*]const f64, incx: isize) usize {
    return iamax(f64, n, x, incx);
}
pub fn icamax(n: isize, x: [*]const Complex(f32), incx: isize) usize {
    return iamax(Complex(f32), n, x, incx);
}
pub fn izamax(n: isize, x: [*]const Complex(f64), incx: isize) usize {
    return iamax(Complex(f64), n, x, incx);
}

pub inline fn iamin(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_isamin(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_idamin(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_icamin(@intCast(n), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_izamin(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/iamin.zig").iamin(T, n, x, incx);
}
pub fn isamin(n: isize, x: [*]const f32, incx: isize) usize {
    return iamin(f32, n, x, incx);
}
pub fn idamin(n: isize, x: [*]const f64, incx: isize) usize {
    return iamin(f64, n, x, incx);
}
pub fn icamin(n: isize, x: [*]const Complex(f32), incx: isize) usize {
    return iamin(Complex(f32), n, x, incx);
}
pub fn izamin(n: isize, x: [*]const Complex(f64), incx: isize) usize {
    return iamin(Complex(f64), n, x, incx);
}

// Level 2 BLAS
pub inline fn gbmv(comptime T: type, order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sgbmv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dgbmv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgbmv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgbmv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/gbmv.zig").gbmv(T, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn sgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return gbmv(f32, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return gbmv(f64, order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn cgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return gbmv(Complex(f32), order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zgbmv(order: Order, transA: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return gbmv(Complex(f64), order, transA, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn gemv(comptime T: type, order: Order, transA: Transpose, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sgemv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dgemv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgemv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgemv(@intFromEnum(order), @intFromEnum(transA), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/gemv.zig").gemv(T, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn sgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return gemv(f32, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return gemv(f64, order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn cgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return gemv(Complex(f32), order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zgemv(order: Order, transA: Transpose, m: isize, n: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return gemv(Complex(f64), order, transA, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn ger(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sger(@intFromEnum(order), @intCast(m), @intCast(n), alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (T == f64) {
                    return ci.cblas_dger(@intFromEnum(order), @intCast(m), @intCast(n), alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/ger.zig").ger(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn sger(order: Order, m: isize, n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, A: [*]f32, lda: isize) void {
    return ger(f32, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn dger(order: Order, m: isize, n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, A: [*]f64, lda: isize) void {
    return ger(f64, order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn gerc(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgerc(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgerc(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/gerc.zig").gerc(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cgerc(order: Order, m: isize, n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, A: [*]Complex(f32), lda: isize) void {
    return gerc(Complex(f32), order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zgerc(order: Order, m: isize, n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, A: [*]Complex(f64), lda: isize) void {
    return gerc(Complex(f64), order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn geru(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgeru(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgeru(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/geru.zig").geru(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cgeru(order: Order, m: isize, n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, A: [*]Complex(f32), lda: isize) void {
    return geru(Complex(f32), order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zgeru(order: Order, m: isize, n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, A: [*]Complex(f64), lda: isize) void {
    return geru(Complex(f64), order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn hbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hbmv.zig").hbmv(T, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn chbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return hbmv(Complex(f32), order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zhbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return hbmv(Complex(f64), order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn hemv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chemv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhemv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hemv.zig").hemv(T, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn chemv(order: Order, uplo: Uplo, n: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return hemv(Complex(f32), order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zhemv(order: Order, uplo: Uplo, n: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return hemv(Complex(f64), order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn her(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cher(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), A, @intCast(lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zher(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/her.zig").her(T, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn cher(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const Complex(f32), incx: isize, A: [*]Complex(f32), lda: isize) void {
    return her(Complex(f32), order, uplo, n, alpha, x, incx, A, lda);
}
pub fn zher(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const Complex(f64), incx: isize, A: [*]Complex(f64), lda: isize) void {
    return her(Complex(f64), order, uplo, n, alpha, x, incx, A, lda);
}

pub inline fn her2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cher2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zher2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/her2.zig").her2(T, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cher2(order: Order, uplo: Uplo, n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, A: [*]Complex(f32), lda: isize) void {
    return her2(Complex(f32), order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zher2(order: Order, uplo: Uplo, n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, A: [*]Complex(f64), lda: isize) void {
    return her2(Complex(f64), order, uplo, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn hpmv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, Ap: [*]const T, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chpmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, Ap, x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhpmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, Ap, x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/hpmv.zig").hpmv(T, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn chpmv(order: Order, uplo: Uplo, n: isize, alpha: Complex(f32), Ap: [*]const Complex(f32), x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return hpmv(Complex(f32), order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn zhpmv(order: Order, uplo: Uplo, n: isize, alpha: Complex(f64), Ap: [*]const Complex(f64), x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return hpmv(Complex(f64), order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}

pub inline fn hpr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chpr(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), Ap);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhpr(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/hpr.zig").hpr(T, order, uplo, n, alpha, x, incx, Ap);
}
pub fn chpr(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const Complex(f32), incx: isize, Ap: [*]Complex(f32)) void {
    return hpr(Complex(f32), order, uplo, n, alpha, x, incx, Ap);
}
pub fn zhpr(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const Complex(f64), incx: isize, Ap: [*]Complex(f64)) void {
    return hpr(Complex(f64), order, uplo, n, alpha, x, incx, Ap);
}

pub inline fn hpr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, Ap: [*]T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chpr2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), Ap);
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhpr2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/hpr2.zig").hpr2(T, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn chpr2(order: Order, uplo: Uplo, n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, Ap: [*]Complex(f32)) void {
    return hpr2(Complex(f32), order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn zhpr2(order: Order, uplo: Uplo, n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, Ap: [*]Complex(f64)) void {
    return hpr2(Complex(f64), order, uplo, n, alpha, x, incx, y, incy, Ap);
}

pub inline fn sbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dsbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/sbmv.zig").sbmv(T, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn ssbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return sbmv(f32, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dsbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return sbmv(f64, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn spmv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, Ap: [*]const T, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sspmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, Ap, x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dspmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, Ap, x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/spmv.zig").spmv(T, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn sspmv(order: Order, uplo: Uplo, n: isize, alpha: f32, Ap: [*]const f32, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return spmv(f32, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}
pub fn dspmv(order: Order, uplo: Uplo, n: isize, alpha: f64, Ap: [*]const f64, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return spmv(f64, order, uplo, n, alpha, Ap, x, incx, beta, y, incy);
}

pub inline fn spr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: Scalar(T), x: [*]const T, incx: isize, Ap: [*]T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sspr(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), Ap);
                } else if (T == f64) {
                    return ci.cblas_dspr(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/spr.zig").spr(T, order, uplo, n, alpha, x, incx, Ap);
}
pub fn sspr(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, Ap: [*]f32) void {
    return spr(f32, order, uplo, n, alpha, x, incx, Ap);
}
pub fn dspr(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, Ap: [*]f64) void {
    return spr(f64, order, uplo, n, alpha, x, incx, Ap);
}

pub inline fn spr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, Ap: [*]T) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sspr2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), y, @intCast(incy), Ap);
                } else if (T == f64) {
                    return ci.cblas_dspr2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), y, @intCast(incy), Ap);
                }
            },
            else => {},
        }
    }

    return @import("blas/spr2.zig").spr2(T, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn sspr2(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, Ap: [*]f32) void {
    return spr2(f32, order, uplo, n, alpha, x, incx, y, incy, Ap);
}
pub fn dspr2(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, Ap: [*]f64) void {
    return spr2(f64, order, uplo, n, alpha, x, incx, y, incy, Ap);
}

pub inline fn symv(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssymv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dsymv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("blas/symv.zig").symv(T, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn ssymv(order: Order, uplo: Uplo, n: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return symv(f32, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dsymv(order: Order, uplo: Uplo, n: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return symv(f64, order, uplo, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub inline fn syr(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssyr(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), A, @intCast(lda));
                } else if (T == f64) {
                    return ci.cblas_dsyr(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr.zig").syr(T, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn ssyr(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, A: [*]f32, lda: isize) void {
    return syr(f32, order, uplo, n, alpha, x, incx, A, lda);
}
pub fn dsyr(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, A: [*]f64, lda: isize) void {
    return syr(f64, order, uplo, n, alpha, x, incx, A, lda);
}

pub inline fn syr2(comptime T: type, order: Order, uplo: Uplo, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssyr2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (T == f64) {
                    return ci.cblas_dsyr2(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr2.zig").syr2(T, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn ssyr2(order: Order, uplo: Uplo, n: isize, alpha: f32, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize, A: [*]f32, lda: isize) void {
    return syr2(f32, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}
pub fn dsyr2(order: Order, uplo: Uplo, n: isize, alpha: f64, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize, A: [*]f64, lda: isize) void {
    return syr2(f64, order, uplo, n, alpha, x, incx, y, incy, A, lda);
}

pub inline fn tbmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_stbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dtbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztbmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tbmv.zig").tbmv(T, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn stbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return tbmv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn dtbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return tbmv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ctbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const Complex(f32), lda: isize, x: [*]Complex(f32), incx: isize) void {
    return tbmv(Complex(f32), order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ztbmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const Complex(f64), lda: isize, x: [*]Complex(f64), incx: isize) void {
    return tbmv(Complex(f64), order, uplo, transA, diag, n, k, A, lda, x, incx);
}

pub inline fn tbsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_stbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dtbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztbsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), @intCast(k), A, @intCast(lda), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tbsv.zig").tbsv(T, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn stbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return tbsv(f32, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn dtbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return tbsv(f64, order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ctbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const Complex(f32), lda: isize, x: [*]Complex(f32), incx: isize) void {
    return tbsv(Complex(f32), order, uplo, transA, diag, n, k, A, lda, x, incx);
}
pub fn ztbsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, k: isize, A: [*]const Complex(f64), lda: isize, x: [*]Complex(f64), incx: isize) void {
    return tbsv(Complex(f64), order, uplo, transA, diag, n, k, A, lda, x, incx);
}

pub inline fn tpmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const T, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_stpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dtpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztpmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tpmv.zig").tpmv(T, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn stpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f32, x: [*]f32, incx: isize) void {
    return tpmv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn dtpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f64, x: [*]f64, incx: isize) void {
    return tpmv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ctpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const Complex(f32), x: [*]Complex(f32), incx: isize) void {
    return tpmv(Complex(f32), order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ztpmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const Complex(f64), x: [*]Complex(f64), incx: isize) void {
    return tpmv(Complex(f64), order, uplo, transA, diag, n, Ap, x, incx);
}

pub inline fn tpsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const T, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_stpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dtpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztpsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), Ap, x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/tpsv.zig").tpsv(T, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn stpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f32, x: [*]f32, incx: isize) void {
    return tpsv(f32, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn dtpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const f64, x: [*]f64, incx: isize) void {
    return tpsv(f64, order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ctpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const Complex(f32), x: [*]Complex(f32), incx: isize) void {
    return tpsv(Complex(f32), order, uplo, transA, diag, n, Ap, x, incx);
}
pub fn ztpsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, Ap: [*]const Complex(f64), x: [*]Complex(f64), incx: isize) void {
    return tpsv(Complex(f64), order, uplo, transA, diag, n, Ap, x, incx);
}

pub inline fn trmv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_strmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dtrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrmv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/trmv.zig").trmv(T, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn strmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return trmv(f32, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn dtrmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return trmv(f64, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ctrmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const Complex(f32), lda: isize, x: [*]Complex(f32), incx: isize) void {
    return trmv(Complex(f32), order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ztrmv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const Complex(f64), lda: isize, x: [*]Complex(f64), incx: isize) void {
    return trmv(Complex(f64), order, uplo, transA, diag, n, A, lda, x, incx);
}

pub inline fn trsv(comptime T: type, order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const T, lda: isize, x: [*]T, incx: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_strsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dtrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrsv(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(n), A, @intCast(lda), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("blas/trsv.zig").trsv(T, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn strsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f32, lda: isize, x: [*]f32, incx: isize) void {
    return trsv(f32, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn dtrsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const f64, lda: isize, x: [*]f64, incx: isize) void {
    return trsv(f64, order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ctrsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const Complex(f32), lda: isize, x: [*]Complex(f32), incx: isize) void {
    return trsv(Complex(f32), order, uplo, transA, diag, n, A, lda, x, incx);
}
pub fn ztrsv(order: Order, uplo: Uplo, transA: Transpose, diag: Diag, n: isize, A: [*]const Complex(f64), lda: isize, x: [*]Complex(f64), incx: isize) void {
    return trsv(Complex(f64), order, uplo, transA, diag, n, A, lda, x, incx);
}

// Level 3 BLAS
pub inline fn gemm(comptime T: type, order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), @intCast(m), @intCast(n), @intCast(k), alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                } else if (T == f64) {
                    return ci.cblas_dgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), @intCast(m), @intCast(n), @intCast(k), alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), @intCast(m), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zgemm(@intFromEnum(order), @intFromEnum(transA), @intFromEnum(transB), @intCast(m), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/gemm.zig").gemm(T, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn sgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]const f32, ldb: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return gemm(f32, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn dgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]const f64, ldb: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return gemm(f64, order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn cgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]const Complex(f32), ldb: isize, beta: Complex(f32), C: [*]Complex(f32), ldc: isize) void {
    return gemm(Complex(f32), order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zgemm(order: Order, transA: Transpose, transB: Transpose, m: isize, n: isize, k: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]const Complex(f64), ldb: isize, beta: Complex(f64), C: [*]Complex(f64), ldc: isize) void {
    return gemm(Complex(f64), order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn hemm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_chemm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zhemm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/hemm.zig").hemm(T, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn chemm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]const Complex(f32), ldb: isize, beta: Complex(f32), C: [*]Complex(f32), ldc: isize) void {
    return hemm(Complex(f32), order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zhemm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]const Complex(f64), ldb: isize, beta: Complex(f64), C: [*]Complex(f64), ldc: isize) void {
    return hemm(Complex(f64), order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn herk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Scalar(T), A: [*]const T, lda: isize, beta: Scalar(T), C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cherk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), alpha, A, @intCast(lda), beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zherk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), alpha, A, @intCast(lda), beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/herk.zig").herk(T, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn cherk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f32, A: [*]const Complex(f32), lda: isize, beta: f32, C: [*]Complex(f32), ldc: isize) void {
    return herk(Complex(f32), order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn zherk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f64, A: [*]const Complex(f64), lda: isize, beta: f64, C: [*]Complex(f64), ldc: isize) void {
    return herk(Complex(f64), order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

pub inline fn her2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: Scalar(T), C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_cher2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zher2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/her2k.zig").her2k(T, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn cher2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]const Complex(f32), ldb: isize, beta: f32, C: [*]Complex(f32), ldc: isize) void {
    return her2k(Complex(f32), order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zher2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]const Complex(f64), ldb: isize, beta: f64, C: [*]Complex(f64), ldc: isize) void {
    return her2k(Complex(f64), order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn symm(comptime T: type, order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intCast(m), @intCast(n), alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                } else if (T == f64) {
                    return ci.cblas_dsymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intCast(m), @intCast(n), alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zsymm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/symm.zig").symm(T, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn ssymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]const f32, ldb: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return symm(f32, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn dsymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]const f64, ldb: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return symm(f64, order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn csymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]const Complex(f32), ldb: isize, beta: Complex(f32), C: [*]Complex(f32), ldc: isize) void {
    return symm(Complex(f32), order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zsymm(order: Order, side: Side, uplo: Uplo, m: isize, n: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]const Complex(f64), ldb: isize, beta: Complex(f64), C: [*]Complex(f64), ldc: isize) void {
    return symm(Complex(f64), order, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn syrk(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), alpha, A, @intCast(lda), beta, C, @intCast(ldc));
                } else if (T == f64) {
                    return ci.cblas_dsyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), alpha, A, @intCast(lda), beta, C, @intCast(ldc));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), &beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zsyrk(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), &beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/syrk.zig").syrk(T, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn ssyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return syrk(f32, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn dsyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return syrk(f64, order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn csyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, beta: Complex(f32), C: [*]Complex(f32), ldc: isize) void {
    return syrk(Complex(f32), order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}
pub fn zsyrk(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, beta: Complex(f64), C: [*]Complex(f64), ldc: isize) void {
    return syrk(Complex(f64), order, uplo, trans, n, k, alpha, A, lda, beta, C, ldc);
}

pub inline fn syr2k(comptime T: type, order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, B: [*]const T, ldb: isize, beta: T, C: [*]T, ldc: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_ssyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                } else if (T == f64) {
                    return ci.cblas_dsyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), alpha, A, @intCast(lda), B, @intCast(ldb), beta, C, @intCast(ldc));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_csyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_zsyr2k(@intFromEnum(order), @intFromEnum(uplo), @intFromEnum(trans), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), B, @intCast(ldb), &beta, C, @intCast(ldc));
                }
            },
            else => {},
        }
    }

    return @import("blas/syr2k.zig").syr2k(T, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn ssyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]const f32, ldb: isize, beta: f32, C: [*]f32, ldc: isize) void {
    return syr2k(f32, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn dsyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]const f64, ldb: isize, beta: f64, C: [*]f64, ldc: isize) void {
    return syr2k(f64, order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn csyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]const Complex(f32), ldb: isize, beta: Complex(f32), C: [*]Complex(f32), ldc: isize) void {
    return syr2k(Complex(f32), order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}
pub fn zsyr2k(order: Order, uplo: Uplo, trans: Transpose, n: isize, k: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]const Complex(f64), ldb: isize, beta: Complex(f64), C: [*]Complex(f64), ldc: isize) void {
    return syr2k(Complex(f64), order, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

pub inline fn trmm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_strmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), alpha, A, @intCast(lda), B, @intCast(ldb));
                } else if (T == f64) {
                    return ci.cblas_dtrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), alpha, A, @intCast(lda), B, @intCast(ldb));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrmm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb));
                }
            },
            else => {},
        }
    }

    return @import("blas/trmm.zig").trmm(T, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn strmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]f32, ldb: isize) void {
    return trmm(f32, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn dtrmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]f64, ldb: isize) void {
    return trmm(f64, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ctrmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]Complex(f32), ldb: isize) void {
    return trmm(Complex(f32), order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ztrmm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]Complex(f64), ldb: isize) void {
    return trmm(Complex(f64), order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

pub inline fn trsm(comptime T: type, order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, B: [*]T, ldb: isize) void {
    const supported = core.types.numericType(T);

    if (options.link_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_strsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), alpha, A, @intCast(lda), B, @intCast(ldb));
                } else if (T == f64) {
                    return ci.cblas_dtrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), alpha, A, @intCast(lda), B, @intCast(ldb));
                }
            },
            .Complex => {
                if (Scalar(T) == f32) {
                    return ci.cblas_ctrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb));
                } else if (Scalar(T) == f64) {
                    return ci.cblas_ztrsm(@intFromEnum(order), @intFromEnum(side), @intFromEnum(uplo), @intFromEnum(transA), @intFromEnum(diag), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), B, @intCast(ldb));
                }
            },
            else => {},
        }
    }

    return @import("blas/trsm.zig").trsm(T, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn strsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f32, A: [*]const f32, lda: isize, B: [*]f32, ldb: isize) void {
    return trsm(f32, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn dtrsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: f64, A: [*]const f64, lda: isize, B: [*]f64, ldb: isize) void {
    return trsm(f64, order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ctrsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, B: [*]Complex(f32), ldb: isize) void {
    return trsm(Complex(f32), order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}
pub fn ztrsm(order: Order, side: Side, uplo: Uplo, transA: Transpose, diag: Diag, m: isize, n: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, B: [*]Complex(f64), ldb: isize) void {
    return trsm(Complex(f64), order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb);
}

test {
    const test_level_1 = true;
    const test_level_2 = true;
    const test_level_3 = true;

    if (test_level_1) {
        _ = @import("blas/asum.zig");
        _ = @import("blas/axpy.zig");
        _ = @import("blas/copy.zig");
        _ = @import("blas/dot.zig");
        _ = @import("blas/dotc.zig");
        _ = @import("blas/dotc_sub.zig");
        _ = @import("blas/dotu.zig");
        _ = @import("blas/dotu_sub.zig");
        _ = @import("blas/nrm2.zig");
        _ = @import("blas/rot.zig");
        _ = @import("blas/rotg.zig");
        _ = @import("blas/rotm.zig");
        _ = @import("blas/rotmg.zig");
        _ = @import("blas/scal.zig");
        _ = @import("blas/swap.zig");
        _ = @import("blas/iamax.zig");
        _ = @import("blas/iamin.zig");
    }

    if (test_level_2) {
        _ = @import("blas/gbmv.zig");
        _ = @import("blas/gemv.zig");
        _ = @import("blas/ger.zig");
        _ = @import("blas/gerc.zig");
        _ = @import("blas/geru.zig");
        _ = @import("blas/hbmv.zig");
        _ = @import("blas/hemv.zig");
        _ = @import("blas/her.zig");
        _ = @import("blas/her2.zig");
        _ = @import("blas/hpmv.zig");
        _ = @import("blas/hpr.zig");
        _ = @import("blas/hpr2.zig");
        _ = @import("blas/sbmv.zig");
        _ = @import("blas/spmv.zig");
        _ = @import("blas/spr.zig");
        _ = @import("blas/spr2.zig");
        _ = @import("blas/symv.zig");
        _ = @import("blas/syr.zig");
        _ = @import("blas/syr2.zig");
        _ = @import("blas/tbmv.zig");
        _ = @import("blas/tbsv.zig");
        _ = @import("blas/tpmv.zig");
        _ = @import("blas/tpsv.zig");
        _ = @import("blas/trmv.zig");
        _ = @import("blas/trsv.zig");
    }

    if (test_level_3) {
        _ = @import("blas/gemm.zig");
        _ = @import("blas/hemm.zig");
        _ = @import("blas/herk.zig");
        _ = @import("blas/her2k.zig");
        _ = @import("blas/symm.zig");
        _ = @import("blas/syrk.zig");
        _ = @import("blas/syr2k.zig");
        _ = @import("blas/trmm.zig");
        _ = @import("blas/trsm.zig");
    }
}
