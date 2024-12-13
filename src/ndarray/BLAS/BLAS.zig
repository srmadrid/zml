const std = @import("std");
const core = @import("../../core/core.zig");
const ci = @import("../../c.zig");
const options = @import("options");
const Complex = std.math.Complex;
const ndarray = @import("../ndarray.zig");
const Order = ndarray.Order;
const Transpose = ndarray.Transpose;
const Uplo = ndarray.Uplo;
const Diag = ndarray.Diag;
const Side = ndarray.Side;

const scalar = core.supported.scalar;

// Level 1 BLAS
pub inline fn asum(comptime T: type, n: isize, x: [*]const T, incx: isize) scalar(T) {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sasum(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dasum(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_scasum(@intCast(n), x, @intCast(incx));
                } else if (scalar(T) == f64) {
                    return ci.cblas_dzasum(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("asum.zig").asum(T, n, x, incx);
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

pub inline fn axpy(comptime T: type, n: isize, a: T, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_saxpy(@intCast(n), a, x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_daxpy(@intCast(n), a, x, @intCast(incx), y, @intCast(incy));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_caxpy(@intCast(n), &a, x, @intCast(incx), y, @intCast(incy));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zaxpy(@intCast(n), &a, x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("axpy.zig").axpy(T, n, a, x, incx, y, incy);
}
pub fn saxpy(n: isize, a: f32, x: [*]const f32, incx: isize, y: [*]f32, incy: isize) void {
    return axpy(f32, n, a, x, incx, y, incy);
}
pub fn daxpy(n: isize, a: f64, x: [*]const f64, incx: isize, y: [*]f64, incy: isize) void {
    return axpy(f64, n, a, x, incx, y, incy);
}
pub fn caxpy(n: isize, a: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]Complex(f32), incy: isize) void {
    return axpy(Complex(f32), n, a, x, incx, y, incy);
}
pub fn zaxpy(n: isize, a: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]Complex(f64), incy: isize) void {
    return axpy(Complex(f64), n, a, x, incx, y, incy);
}

pub fn copy(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_scopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dcopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_ccopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zcopy(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("copy.zig").copy(T, n, x, incx, y, incy);
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

pub fn dot(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) scalar(T) {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
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

    return @import("dot.zig").dot(T, n, x, incx, y, incy);
}
pub fn sdot(n: isize, x: [*]const f32, incx: isize, y: [*]const f32, incy: isize) f32 {
    return dot(f32, n, x, incx, y, incy);
}
pub fn ddot(n: isize, x: [*]const f64, incx: isize, y: [*]const f64, incy: isize) f64 {
    return dot(f64, n, x, incx, y, incy);
}

pub fn dotc(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    //const supported = core.supported.whatSupportedNumericType(T);

    //if (options.use_cblas) {
    //    switch (supported) {
    //        .Complex => {
    //            if (scalar(T) == f32) {
    //                return ci.cblas_cdotc(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            } else if (scalar(T) == f64) {
    //                return ci.cblas_zdotc(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            }
    //        },
    //        else => {},
    //    }
    //}

    return @import("dotc.zig").dotc(T, n, x, incx, y, incy);
}
pub fn cdotc(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize) Complex(f32) {
    return dotc(Complex(f32), n, x, incx, y, incy);
}
pub fn zdotc(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize) Complex(f64) {
    return dotc(Complex(f64), n, x, incx, y, incy);
}

pub fn dotc_sub(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize, ret: *T) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cdotc_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                } else if (scalar(T) == f64) {
                    return ci.cblas_zdotc_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                }
            },
            else => {},
        }
    }

    return @import("dotc_sub.zig").dotc_sub(T, n, x, incx, y, incy, ret);
}
pub fn cdotc_sub(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, ret: *Complex(f32)) void {
    return dotc_sub(Complex(f32), n, x, incx, y, incy, ret);
}
pub fn zdotc_sub(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, ret: *Complex(f64)) void {
    return dotc_sub(Complex(f64), n, x, incx, y, incy, ret);
}

pub fn dotu(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize) T {
    //const supported = core.supported.whatSupportedNumericType(T);

    //if (options.use_cblas) {
    //    switch (supported) {
    //        .Complex => {
    //            if (scalar(T) == f32) {
    //                return ci.cblas_cdotu(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            } else if (scalar(T) == f64) {
    //                return ci.cblas_zdotu(@intCast(n), x, @intCast(incx), y, @intCast(incy));
    //            }
    //        },
    //        else => {},
    //    }
    //}

    return @import("dotu.zig").dotu(T, n, x, incx, y, incy);
}
pub fn cdotu(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize) Complex(f32) {
    return dotu(Complex(f32), n, x, incx, y, incy);
}
pub fn zdotu(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize) Complex(f64) {
    return dotu(Complex(f64), n, x, incx, y, incy);
}

pub fn dotu_sub(comptime T: type, n: isize, x: [*]const T, incx: isize, y: [*]const T, incy: isize, ret: *T) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cdotu_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                } else if (scalar(T) == f64) {
                    return ci.cblas_zdotu_sub(@intCast(n), x, @intCast(incx), y, @intCast(incy), ret);
                }
            },
            else => {},
        }
    }

    return @import("dotu_sub.zig").dotu_sub(T, n, x, incx, y, incy, ret);
}
pub fn cdotu_sub(n: isize, x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, ret: *Complex(f32)) void {
    return dotu_sub(Complex(f32), n, x, incx, y, incy, ret);
}
pub fn zdotu_sub(n: isize, x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, ret: *Complex(f64)) void {
    return dotu_sub(Complex(f64), n, x, incx, y, incy, ret);
}

pub fn nrm2(comptime T: type, n: isize, x: [*]const T, incx: isize) scalar(T) {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_snrm2(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dnrm2(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_scnrm2(@intCast(n), x, @intCast(incx));
                } else if (scalar(T) == f64) {
                    return ci.cblas_dznrm2(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("nrm2.zig").nrm2(T, n, x, incx);
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

pub fn rot(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, c: scalar(T), s: scalar(T)) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_srot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                } else if (T == f64) {
                    return ci.cblas_drot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_csrot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                } else if (scalar(T) == f64) {
                    return ci.cblas_zdrot(@intCast(n), x, @intCast(incx), y, @intCast(incy), c, s);
                }
            },
            else => {},
        }
    }

    return @import("rot.zig").rot(T, n, x, incx, y, incy, c, s);
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

pub fn rotg(comptime T: type, a: *T, b: *T, c: *scalar(T), s: *T) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_srotg(a, b, c, s);
                } else if (T == f64) {
                    return ci.cblas_drotg(a, b, c, s);
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_crotg(a, b, c, s);
                } else if (scalar(T) == f64) {
                    return ci.cblas_zrotg(a, b, c, s);
                }
            },
            else => {},
        }
    }

    return @import("rotg.zig").rotg(T, a, b, c, s);
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

pub fn rotm(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize, param: [*]const T) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
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

    return @import("rotm.zig").rotm(T, n, x, incx, y, incy, param);
}
pub fn srotm(n: isize, x: [*]f32, incx: isize, y: [*]f32, incy: isize, param: [*]const f32) void {
    return rotm(f32, n, x, incx, y, incy, param);
}
pub fn drotm(n: isize, x: [*]f64, incx: isize, y: [*]f64, incy: isize, param: [*]const f64) void {
    return rotm(f64, n, x, incx, y, incy, param);
}

pub fn rotmg(comptime T: type, d1: *T, d2: *T, x1: *T, y1: T, param: [*]T) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
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

    return @import("rotmg.zig").rotmg(T, d1, d2, x1, y1, param);
}
pub fn srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*]f32) void {
    return rotmg(f32, d1, d2, x1, y1, param);
}
pub fn drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*]f64) void {
    return rotmg(f64, d1, d2, x1, y1, param);
}

pub fn scal(comptime T: type, n: isize, a: T, x: [*]T, incx: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sscal(@intCast(n), a, x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_dscal(@intCast(n), a, x, @intCast(incx));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cscal(@intCast(n), &a, x, @intCast(incx));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zscal(@intCast(n), &a, x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("scal.zig").scal(T, n, a, x, incx);
}
pub fn sscal(n: isize, a: f32, x: [*]f32, incx: isize) void {
    return scal(f32, n, a, x, incx);
}
pub fn dscal(n: isize, a: f64, x: [*]f64, incx: isize) void {
    return scal(f64, n, a, x, incx);
}
pub fn cscal(n: isize, a: Complex(f32), x: [*]Complex(f32), incx: isize) void {
    return scal(Complex(f32), n, a, x, incx);
}
pub fn zscal(n: isize, a: Complex(f64), x: [*]Complex(f64), incx: isize) void {
    return scal(Complex(f64), n, a, x, incx);
}

pub fn swap(comptime T: type, n: isize, x: [*]T, incx: isize, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zswap(@intCast(n), x, @intCast(incx), y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("swap.zig").swap(T, n, x, incx, y, incy);
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

pub fn iamax(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_isamax(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_idamax(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_icamax(@intCast(n), x, @intCast(incx));
                } else if (scalar(T) == f64) {
                    return ci.cblas_izamax(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("iamax.zig").iamax(T, n, x, incx);
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

pub fn iamin(comptime T: type, n: isize, x: [*]const T, incx: isize) usize {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_isamin(@intCast(n), x, @intCast(incx));
                } else if (T == f64) {
                    return ci.cblas_idamin(@intCast(n), x, @intCast(incx));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_icamin(@intCast(n), x, @intCast(incx));
                } else if (scalar(T) == f64) {
                    return ci.cblas_izamin(@intCast(n), x, @intCast(incx));
                }
            },
            else => {},
        }
    }

    return @import("iamin.zig").iamin(T, n, x, incx);
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
pub fn gbmv(comptime T: type, order: Order, trans: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sgbmv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dgbmv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cgbmv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zgbmv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), @intCast(kl), @intCast(ku), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("gbmv.zig").gbmv(T, order, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn sgbmv(order: Order, trans: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: f32, A: [*]const f32, lda: isize, x: [*]const f32, incx: isize, beta: f32, y: [*]f32, incy: isize) void {
    return gbmv(f32, order, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn dgbmv(order: Order, trans: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: f64, A: [*]const f64, lda: isize, x: [*]const f64, incx: isize, beta: f64, y: [*]f64, incy: isize) void {
    return gbmv(f64, order, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn cgbmv(order: Order, trans: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return gbmv(Complex(f32), order, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zgbmv(order: Order, trans: Transpose, m: isize, n: isize, kl: isize, ku: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return gbmv(Complex(f64), order, trans, m, n, kl, ku, alpha, A, lda, x, incx, beta, y, incy);
}

pub fn gemv(comptime T: type, order: Order, trans: Transpose, m: isize, n: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_sgemv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_dgemv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cgemv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zgemv(@intFromEnum(order), @intFromEnum(trans), @intCast(m), @intCast(n), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("gemv.zig").gemv(T, order, trans, m, n, alpha, A, lda, x, incx, beta, y, incy);
}

pub fn ger(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
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

    return @import("ger.zig").ger(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub fn gerc(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cgerc(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zgerc(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("gerc.zig").gerc(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cgerc(order: Order, m: isize, n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, A: [*]Complex(f32), lda: isize) void {
    return gerc(Complex(f32), order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zgerc(order: Order, m: isize, n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, A: [*]Complex(f64), lda: isize) void {
    return gerc(Complex(f64), order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub fn geru(comptime T: type, order: Order, m: isize, n: isize, alpha: T, x: [*]const T, incx: isize, y: [*]const T, incy: isize, A: [*]T, lda: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_cgeru(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zgeru(@intFromEnum(order), @intCast(m), @intCast(n), &alpha, x, @intCast(incx), y, @intCast(incy), A, @intCast(lda));
                }
            },
            else => {},
        }
    }

    return @import("geru.zig").geru(T, order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn cgeru(order: Order, m: isize, n: isize, alpha: Complex(f32), x: [*]const Complex(f32), incx: isize, y: [*]const Complex(f32), incy: isize, A: [*]Complex(f32), lda: isize) void {
    return geru(Complex(f32), order, m, n, alpha, x, incx, y, incy, A, lda);
}
pub fn zgeru(order: Order, m: isize, n: isize, alpha: Complex(f64), x: [*]const Complex(f64), incx: isize, y: [*]const Complex(f64), incy: isize, A: [*]Complex(f64), lda: isize) void {
    return geru(Complex(f64), order, m, n, alpha, x, incx, y, incy, A, lda);
}

pub fn hbmv(comptime T: type, order: Order, uplo: Uplo, n: isize, k: isize, alpha: T, A: [*]const T, lda: isize, x: [*]const T, incx: isize, beta: T, y: [*]T, incy: isize) void {
    const supported = core.supported.whatSupportedNumericType(T);

    if (options.use_cblas != null) {
        switch (supported) {
            .BuiltinFloat => {
                if (T == f32) {
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                } else if (T == f64) {
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), alpha, A, @intCast(lda), x, @intCast(incx), beta, y, @intCast(incy));
                }
            },
            .Complex => {
                if (scalar(T) == f32) {
                    return ci.cblas_chbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                } else if (scalar(T) == f64) {
                    return ci.cblas_zhbmv(@intFromEnum(order), @intFromEnum(uplo), @intCast(n), @intCast(k), &alpha, A, @intCast(lda), x, @intCast(incx), &beta, y, @intCast(incy));
                }
            },
            else => {},
        }
    }

    return @import("hbmv.zig").hbmv(T, order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn chbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: Complex(f32), A: [*]const Complex(f32), lda: isize, x: [*]const Complex(f32), incx: isize, beta: Complex(f32), y: [*]Complex(f32), incy: isize) void {
    return hbmv(Complex(f32), order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}
pub fn zhbmv(order: Order, uplo: Uplo, n: isize, k: isize, alpha: Complex(f64), A: [*]const Complex(f64), lda: isize, x: [*]const Complex(f64), incx: isize, beta: Complex(f64), y: [*]Complex(f64), incy: isize) void {
    return hbmv(Complex(f64), order, uplo, n, k, alpha, A, lda, x, incx, beta, y, incy);
}

test {
    std.testing.refAllDeclsRecursive(@This());

    // Level 1 BLAS
    _ = @import("asum.zig");
    _ = @import("axpy.zig");
    _ = @import("copy.zig");
    _ = @import("dot.zig");
    _ = @import("dotc.zig");
    _ = @import("dotc_sub.zig");
    _ = @import("dotu.zig");
    _ = @import("dotu_sub.zig");
    _ = @import("nrm2.zig");
    _ = @import("rot.zig");
    _ = @import("rotg.zig");
    _ = @import("rotm.zig");
    _ = @import("rotmg.zig");
    _ = @import("scal.zig");
    _ = @import("swap.zig");
    _ = @import("iamax.zig");
    _ = @import("iamin.zig");

    // Level 2 BLAS
    _ = @import("gbmv.zig");
    _ = @import("gemv.zig");
    _ = @import("ger.zig");
    _ = @import("gerc.zig");
    _ = @import("geru.zig");
    _ = @import("hbmv.zig");
}
