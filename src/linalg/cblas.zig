const complex = @import("../complex.zig");
const cf32 = complex.cf32;
const cf64 = complex.cf64;

pub const row_major: c_int = 101;
pub const col_major: c_int = 102;

pub const no_trans: c_int = 111;
pub const trans: c_int = 112;
pub const conj_trans: c_int = 113;
pub const conj_no_trans: c_int = 114;

pub const upper: c_int = 121;
pub const lower: c_int = 122;

pub const non_unit: c_int = 131;
pub const unit: c_int = 132;

pub const left: c_int = 141;
pub const right: c_int = 142;

// Level 1
extern fn cblas_sasum(n: isize, x: [*c]const f32, incx: isize) f32;
extern fn cblas_dasum(n: isize, x: [*c]const f64, incx: isize) f64;
extern fn cblas_scasum(n: isize, x: [*c]const cf32, incx: isize) f32;
extern fn cblas_dzasum(n: isize, x: [*c]const cf64, incx: isize) f64;
pub const sasum = cblas_sasum;
pub const dasum = cblas_dasum;
pub const scasum = cblas_scasum;
pub const dzasum = cblas_dzasum;

extern fn cblas_saxpy(n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]f32, incy: isize) void;
extern fn cblas_daxpy(n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]f64, incy: isize) void;
extern fn cblas_caxpy(n: isize, alpha: [*c]const cf32, x: [*c]const cf32, incx: isize, y: [*c]cf32, incy: isize) void;
extern fn cblas_zaxpy(n: isize, alpha: [*c]const cf64, x: [*c]const cf64, incx: isize, y: [*c]cf64, incy: isize) void;
pub const saxpy = cblas_saxpy;
pub const daxpy = cblas_daxpy;
pub const caxpy = cblas_caxpy;
pub const zaxpy = cblas_zaxpy;

extern fn cblas_scopy(n: isize, x: [*c]const f32, incx: isize, y: [*c]f32, incy: isize) void;
extern fn cblas_dcopy(n: isize, x: [*c]const f64, incx: isize, y: [*c]f64, incy: isize) void;
extern fn cblas_ccopy(n: isize, x: [*c]const cf32, incx: isize, y: [*c]cf32, incy: isize) void;
extern fn cblas_zcopy(n: isize, x: [*c]const cf64, incx: isize, y: [*c]cf64, incy: isize) void;
pub const scopy = cblas_scopy;
pub const dcopy = cblas_dcopy;
pub const ccopy = cblas_ccopy;
pub const zcopy = cblas_zcopy;

extern fn cblas_sdot(n: isize, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize) f32;
extern fn cblas_ddot(n: isize, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize) f64;
pub const sdot = cblas_sdot;
pub const ddot = cblas_ddot;

extern fn cblas_cdotc_sub(n: isize, x: [*c]const cf32, incx: isize, y: [*c]const cf32, incy: isize, dotc: [*c]cf32) void;
extern fn cblas_zdotc_sub(n: isize, x: [*c]const cf64, incx: isize, y: [*c]const cf64, incy: isize, dotc: [*c]cf64) void;
extern fn cblas_cdotu_sub(n: isize, x: [*c]const cf32, incx: isize, y: [*c]const cf32, incy: isize, dotu: [*c]cf32) void;
extern fn cblas_zdotu_sub(n: isize, x: [*c]const cf64, incx: isize, y: [*c]const cf64, incy: isize, dotu: [*c]cf64) void;
pub const cdotc_sub = cblas_cdotc_sub;
pub const zdotc_sub = cblas_zdotc_sub;
pub const cdotu_sub = cblas_cdotu_sub;
pub const zdotu_sub = cblas_zdotu_sub;

extern fn cblas_snrm2(n: isize, x: [*c]const f32, incx: isize) f32;
extern fn cblas_dnrm2(n: isize, x: [*c]const f64, incx: isize) f64;
extern fn cblas_scnrm2(n: isize, x: [*c]const cf32, incx: isize) f32;
extern fn cblas_dznrm2(n: isize, x: [*c]const cf64, incx: isize) f64;
pub const snrm2 = cblas_snrm2;
pub const dnrm2 = cblas_dnrm2;
pub const scnrm2 = cblas_scnrm2;
pub const dznrm2 = cblas_dznrm2;

extern fn cblas_srot(n: isize, x: [*c]f32, incx: isize, y: [*c]f32, incy: isize, c: f32, s: f32) void;
extern fn cblas_drot(n: isize, x: [*c]f64, incx: isize, y: [*c]f64, incy: isize, c: f64, s: f64) void;
extern fn cblas_csrot(n: isize, x: [*c]cf32, incx: isize, y: [*c]cf32, incy: isize, c: f32, s: f32) void;
extern fn cblas_zdrot(n: isize, x: [*c]cf64, incx: isize, y: [*c]cf64, incy: isize, c: f64, s: f64) void;
extern fn cblas_crot(n: isize, x: [*c]cf32, incx: isize, y: [*c]cf32, incy: isize, c: f32, s: [*c]const cf32) void;
extern fn cblas_zrot(n: isize, x: [*c]cf64, incx: isize, y: [*c]cf64, incy: isize, c: f64, s: [*c]const cf64) void;
pub const srot = cblas_srot;
pub const drot = cblas_drot;
pub const csrot = cblas_csrot;
pub const zdrot = cblas_zdrot;
pub const crot = cblas_crot;
pub const zrot = cblas_zrot;

extern fn cblas_srotg(a: *f32, b: *f32, c: *f32, s: *f32) void;
extern fn cblas_drotg(a: *f64, b: *f64, c: *f64, s: *f64) void;
extern fn cblas_crotg(a: *cf32, b: *cf32, c: *f32, s: *cf32) void;
extern fn cblas_zrotg(a: *cf64, b: *cf64, c: *f64, s: *cf64) void;
pub const srotg = cblas_srotg;
pub const drotg = cblas_drotg;
pub const crotg = cblas_crotg;
pub const zrotg = cblas_zrotg;

extern fn cblas_srotm(n: isize, x: [*c]f32, incx: isize, y: [*c]f32, incy: isize, param: [*c]const f32) void;
extern fn cblas_drotm(n: isize, x: [*c]f64, incx: isize, y: [*c]f64, incy: isize, param: [*c]const f64) void;
pub const srotm = cblas_srotm;
pub const drotm = cblas_drotm;

extern fn cblas_srotmg(d1: *f32, d2: *f32, x1: *f32, y1: f32, param: [*c]f32) void;
extern fn cblas_drotmg(d1: *f64, d2: *f64, x1: *f64, y1: f64, param: [*c]f64) void;
pub const srotmg = cblas_srotmg;
pub const drotmg = cblas_drotmg;

extern fn cblas_sscal(n: isize, alpha: f32, x: [*c]f32, incx: isize) void;
extern fn cblas_dscal(n: isize, alpha: f64, x: [*c]f64, incx: isize) void;
extern fn cblas_cscal(n: isize, alpha: [*c]const cf32, x: [*c]cf32, incx: isize) void;
extern fn cblas_zscal(n: isize, alpha: [*c]const cf64, x: [*c]cf64, incx: isize) void;
extern fn cblas_csscal(n: isize, alpha: f32, x: [*c]cf32, incx: isize) void;
extern fn cblas_zdscal(n: isize, alpha: f64, x: [*c]cf64, incx: isize) void;
pub const sscal = cblas_sscal;
pub const dscal = cblas_dscal;
pub const cscal = cblas_cscal;
pub const zscal = cblas_zscal;
pub const csscal = cblas_csscal;
pub const zdscal = cblas_zdscal;

extern fn cblas_sswap(n: isize, x: [*c]f32, incx: isize, y: [*c]f32, incy: isize) void;
extern fn cblas_dswap(n: isize, x: [*c]f64, incx: isize, y: [*c]f64, incy: isize) void;
extern fn cblas_cswap(n: isize, x: [*c]cf32, incx: isize, y: [*c]cf32, incy: isize) void;
extern fn cblas_zswap(n: isize, x: [*c]cf64, incx: isize, y: [*c]cf64, incy: isize) void;
pub const sswap = cblas_sswap;
pub const dswap = cblas_dswap;
pub const cswap = cblas_cswap;
pub const zswap = cblas_zswap;

extern fn cblas_isamax(n: isize, x: [*c]const f32, incx: isize) usize;
extern fn cblas_idamax(n: isize, x: [*c]const f64, incx: isize) usize;
extern fn cblas_icamax(n: isize, x: [*c]const cf32, incx: isize) usize;
extern fn cblas_izamax(n: isize, x: [*c]const cf64, incx: isize) usize;
pub const isamax = cblas_isamax;
pub const idamax = cblas_idamax;
pub const icamax = cblas_icamax;
pub const izamax = cblas_izamax;

extern fn cblas_isamin(n: isize, x: [*c]const f32, incx: isize) usize;
extern fn cblas_idamin(n: isize, x: [*c]const f64, incx: isize) usize;
extern fn cblas_icamin(n: isize, x: [*c]const cf32, incx: isize) usize;
extern fn cblas_izamin(n: isize, x: [*c]const cf64, incx: isize) usize;
pub const isamin = cblas_isamin;
pub const idamin = cblas_idamin;
pub const icamin = cblas_icamin;
pub const izamin = cblas_izamin;

// Level 2
extern fn cblas_sgbmv(layout: c_int, transa: c_int, m: isize, n: isize, kl: isize, ku: isize, alpha: f32, a: [*c]const f32, lda: isize, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void;
extern fn cblas_dgbmv(layout: c_int, transa: c_int, m: isize, n: isize, kl: isize, ku: isize, alpha: f64, a: [*c]const f64, lda: isize, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void;
extern fn cblas_cgbmv(layout: c_int, transa: c_int, m: isize, n: isize, kl: isize, ku: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, x: [*c]const cf32, incx: isize, beta: [*c]const cf32, y: [*c]cf32, incy: isize) void;
extern fn cblas_zgbmv(layout: c_int, transa: c_int, m: isize, n: isize, kl: isize, ku: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, x: [*c]const cf64, incx: isize, beta: [*c]const cf64, y: [*c]cf64, incy: isize) void;
pub const sgbmv = cblas_sgbmv;
pub const dgbmv = cblas_dgbmv;
pub const cgbmv = cblas_cgbmv;
pub const zgbmv = cblas_zgbmv;

extern fn cblas_sgemv(layout: c_int, transa: c_int, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void;
extern fn cblas_dgemv(layout: c_int, transa: c_int, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void;
extern fn cblas_cgemv(layout: c_int, transa: c_int, m: isize, n: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, x: [*c]const cf32, incx: isize, beta: [*c]const cf32, y: [*c]cf32, incy: isize) void;
extern fn cblas_zgemv(layout: c_int, transa: c_int, m: isize, n: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, x: [*c]const cf64, incx: isize, beta: [*c]const cf64, y: [*c]cf64, incy: isize) void;
pub const sgemv = cblas_sgemv;
pub const dgemv = cblas_dgemv;
pub const cgemv = cblas_cgemv;
pub const zgemv = cblas_zgemv;

extern fn cblas_sger(layout: c_int, m: isize, n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize, a: [*c]f32, lda: isize) void;
extern fn cblas_dger(layout: c_int, m: isize, n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize, a: [*c]f64, lda: isize) void;
extern fn cblas_cgerc(layout: c_int, m: isize, n: isize, alpha: [*c]const cf32, x: [*c]const cf32, incx: isize, y: [*c]const cf32, incy: isize, a: [*c]cf32, lda: isize) void;
extern fn cblas_zgerc(layout: c_int, m: isize, n: isize, alpha: [*c]const cf64, x: [*c]const cf64, incx: isize, y: [*c]const cf64, incy: isize, a: [*c]cf64, lda: isize) void;
extern fn cblas_cgeru(layout: c_int, m: isize, n: isize, alpha: [*c]const cf32, x: [*c]const cf32, incx: isize, y: [*c]const cf32, incy: isize, a: [*c]cf32, lda: isize) void;
extern fn cblas_zgeru(layout: c_int, m: isize, n: isize, alpha: [*c]const cf64, x: [*c]const cf64, incx: isize, y: [*c]const cf64, incy: isize, a: [*c]cf64, lda: isize) void;
pub const sger = cblas_sger;
pub const dger = cblas_dger;
pub const cgerc = cblas_cgerc;
pub const zgerc = cblas_zgerc;
pub const cgeru = cblas_cgeru;
pub const zgeru = cblas_zgeru;

extern fn cblas_chbmv(layout: c_int, uplo: c_int, n: isize, k: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, x: [*c]const cf32, incx: isize, beta: [*c]const cf32, y: [*c]cf32, incy: isize) void;
extern fn cblas_zhbmv(layout: c_int, uplo: c_int, n: isize, k: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, x: [*c]const cf64, incx: isize, beta: [*c]const cf64, y: [*c]cf64, incy: isize) void;
pub const chbmv = cblas_chbmv;
pub const zhbmv = cblas_zhbmv;

extern fn cblas_chemv(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, x: [*c]const cf32, incx: isize, beta: [*c]const cf32, y: [*c]cf32, incy: isize) void;
extern fn cblas_zhemv(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, x: [*c]const cf64, incx: isize, beta: [*c]const cf64, y: [*c]cf64, incy: isize) void;
pub const chemv = cblas_chemv;
pub const zhemv = cblas_zhemv;

extern fn cblas_cher(layout: c_int, uplo: c_int, n: isize, alpha: f32, x: [*c]const cf32, incx: isize, a: [*c]cf32, lda: isize) void;
extern fn cblas_zher(layout: c_int, uplo: c_int, n: isize, alpha: f64, x: [*c]const cf64, incx: isize, a: [*c]cf64, lda: isize) void;
pub const cher = cblas_cher;
pub const zher = cblas_zher;

extern fn cblas_cher2(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf32, x: [*c]const cf32, incx: isize, y: [*c]const cf32, incy: isize, a: [*c]cf32, lda: isize) void;
extern fn cblas_zher2(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf64, x: [*c]const cf64, incx: isize, y: [*c]const cf64, incy: isize, a: [*c]cf64, lda: isize) void;
pub const cher2 = cblas_cher2;
pub const zher2 = cblas_zher2;

extern fn cblas_chpmv(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf32, ap: [*c]const cf32, x: [*c]const cf32, incx: isize, beta: [*c]const cf32, y: [*c]cf32, incy: isize) void;
extern fn cblas_zhpmv(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf64, ap: [*c]const cf64, x: [*c]const cf64, incx: isize, beta: [*c]const cf64, y: [*c]cf64, incy: isize) void;
pub const chpmv = cblas_chpmv;
pub const zhpmv = cblas_zhpmv;

extern fn cblas_chpr(layout: c_int, uplo: c_int, n: isize, alpha: f32, x: [*c]const cf32, incx: isize, ap: [*c]cf32) void;
extern fn cblas_zhpr(layout: c_int, uplo: c_int, n: isize, alpha: f64, x: [*c]const cf64, incx: isize, ap: [*c]cf64) void;
pub const chpr = cblas_chpr;
pub const zhpr = cblas_zhpr;

extern fn cblas_chpr2(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf32, x: [*c]const cf32, incx: isize, y: [*c]const cf32, incy: isize, ap: [*c]cf32) void;
extern fn cblas_zhpr2(layout: c_int, uplo: c_int, n: isize, alpha: [*c]const cf64, x: [*c]const cf64, incx: isize, y: [*c]const cf64, incy: isize, ap: [*c]cf64) void;
pub const chpr2 = cblas_chpr2;
pub const zhpr2 = cblas_zhpr2;

extern fn cblas_ssbmv(layout: c_int, uplo: c_int, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void;
extern fn cblas_dsbmv(layout: c_int, uplo: c_int, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void;
pub const ssbmv = cblas_ssbmv;
pub const dsbmv = cblas_dsbmv;

extern fn cblas_sspmv(layout: c_int, uplo: c_int, n: isize, alpha: f32, ap: [*c]const f32, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void;
extern fn cblas_dspmv(layout: c_int, uplo: c_int, n: isize, alpha: f64, ap: [*c]const f64, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void;
pub const sspmv = cblas_sspmv;
pub const dspmv = cblas_dspmv;

extern fn cblas_sspr(layout: c_int, uplo: c_int, n: isize, alpha: f32, x: [*c]const f32, incx: isize, ap: [*c]f32) void;
extern fn cblas_dspr(layout: c_int, uplo: c_int, n: isize, alpha: f64, x: [*c]const f64, incx: isize, ap: [*c]f64) void;
pub const sspr = cblas_sspr;
pub const dspr = cblas_dspr;

extern fn cblas_sspr2(layout: c_int, uplo: c_int, n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize, ap: [*c]f32) void;
extern fn cblas_dspr2(layout: c_int, uplo: c_int, n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize, ap: [*c]f64) void;
pub const sspr2 = cblas_sspr2;
pub const dspr2 = cblas_dspr2;

extern fn cblas_ssymv(layout: c_int, uplo: c_int, n: isize, alpha: f32, a: [*c]const f32, lda: isize, x: [*c]const f32, incx: isize, beta: f32, y: [*c]f32, incy: isize) void;
extern fn cblas_dsymv(layout: c_int, uplo: c_int, n: isize, alpha: f64, a: [*c]const f64, lda: isize, x: [*c]const f64, incx: isize, beta: f64, y: [*c]f64, incy: isize) void;
pub const ssymv = cblas_ssymv;
pub const dsymv = cblas_dsymv;

extern fn cblas_ssyr(layout: c_int, uplo: c_int, n: isize, alpha: f32, x: [*c]const f32, incx: isize, a: [*c]f32, lda: isize) void;
extern fn cblas_dsyr(layout: c_int, uplo: c_int, n: isize, alpha: f64, x: [*c]const f64, incx: isize, a: [*c]f64, lda: isize) void;
pub const ssyr = cblas_ssyr;
pub const dsyr = cblas_dsyr;

extern fn cblas_ssyr2(layout: c_int, uplo: c_int, n: isize, alpha: f32, x: [*c]const f32, incx: isize, y: [*c]const f32, incy: isize, a: [*c]f32, lda: isize) void;
extern fn cblas_dsyr2(layout: c_int, uplo: c_int, n: isize, alpha: f64, x: [*c]const f64, incx: isize, y: [*c]const f64, incy: isize, a: [*c]f64, lda: isize) void;
pub const ssyr2 = cblas_ssyr2;
pub const dsyr2 = cblas_dsyr2;

extern fn cblas_stbmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const f32, lda: isize, x: [*c]f32, incx: isize) void;
extern fn cblas_dtbmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const f64, lda: isize, x: [*c]f64, incx: isize) void;
extern fn cblas_ctbmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const cf32, lda: isize, x: [*c]cf32, incx: isize) void;
extern fn cblas_ztbmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const cf64, lda: isize, x: [*c]cf64, incx: isize) void;
pub const stbmv = cblas_stbmv;
pub const dtbmv = cblas_dtbmv;
pub const ctbmv = cblas_ctbmv;
pub const ztbmv = cblas_ztbmv;

extern fn cblas_stbsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const f32, lda: isize, x: [*c]f32, incx: isize) void;
extern fn cblas_dtbsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const f64, lda: isize, x: [*c]f64, incx: isize) void;
extern fn cblas_ctbsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const cf32, lda: isize, x: [*c]cf32, incx: isize) void;
extern fn cblas_ztbsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, k: isize, a: [*c]const cf64, lda: isize, x: [*c]cf64, incx: isize) void;
pub const stbsv = cblas_stbsv;
pub const dtbsv = cblas_dtbsv;
pub const ctbsv = cblas_ctbsv;
pub const ztbsv = cblas_ztbsv;

extern fn cblas_stpmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const f32, x: [*c]f32, incx: isize) void;
extern fn cblas_dtpmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const f64, x: [*c]f64, incx: isize) void;
extern fn cblas_ctpmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const cf32, x: [*c]cf32, incx: isize) void;
extern fn cblas_ztpmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const cf64, x: [*c]cf64, incx: isize) void;
pub const stpmv = cblas_stpmv;
pub const dtpmv = cblas_dtpmv;
pub const ctpmv = cblas_ctpmv;
pub const ztpmv = cblas_ztpmv;

extern fn cblas_stpsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const f32, x: [*c]f32, incx: isize) void;
extern fn cblas_dtpsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const f64, x: [*c]f64, incx: isize) void;
extern fn cblas_ctpsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const cf32, x: [*c]cf32, incx: isize) void;
extern fn cblas_ztpsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, ap: [*c]const cf64, x: [*c]cf64, incx: isize) void;
pub const stpsv = cblas_stpsv;
pub const dtpsv = cblas_dtpsv;
pub const ctpsv = cblas_ctpsv;
pub const ztpsv = cblas_ztpsv;

extern fn cblas_strmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const f32, lda: isize, x: [*c]f32, incx: isize) void;
extern fn cblas_dtrmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const f64, lda: isize, x: [*c]f64, incx: isize) void;
extern fn cblas_ctrmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const cf32, lda: isize, x: [*c]cf32, incx: isize) void;
extern fn cblas_ztrmv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const cf64, lda: isize, x: [*c]cf64, incx: isize) void;
pub const strmv = cblas_strmv;
pub const dtrmv = cblas_dtrmv;
pub const ctrmv = cblas_ctrmv;
pub const ztrmv = cblas_ztrmv;

extern fn cblas_strsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const f32, lda: isize, x: [*c]f32, incx: isize) void;
extern fn cblas_dtrsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const f64, lda: isize, x: [*c]f64, incx: isize) void;
extern fn cblas_ctrsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const cf32, lda: isize, x: [*c]cf32, incx: isize) void;
extern fn cblas_ztrsv(layout: c_int, uplo: c_int, transa: c_int, diag: c_int, n: isize, a: [*c]const cf64, lda: isize, x: [*c]cf64, incx: isize) void;
pub const strsv = cblas_strsv;
pub const dtrsv = cblas_dtrsv;
pub const ctrsv = cblas_ctrsv;
pub const ztrsv = cblas_ztrsv;

// Level 3
extern fn cblas_sgemm(layout: c_int, transa: c_int, transb: c_int, m: isize, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, beta: f32, c: [*c]f32, ldc: isize) void;
extern fn cblas_dgemm(layout: c_int, transa: c_int, transb: c_int, m: isize, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, beta: f64, c: [*c]f64, ldc: isize) void;
extern fn cblas_cgemm(layout: c_int, transa: c_int, transb: c_int, m: isize, n: isize, k: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, beta: [*c]const cf32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zgemm(layout: c_int, transa: c_int, transb: c_int, m: isize, n: isize, k: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, beta: [*c]const cf64, c: [*c]cf64, ldc: isize) void;
pub const sgemm = cblas_sgemm;
pub const dgemm = cblas_dgemm;
pub const cgemm = cblas_cgemm;
pub const zgemm = cblas_zgemm;

extern fn cblas_chemm(layout: c_int, side: c_int, uplo: c_int, m: isize, n: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, beta: [*c]const cf32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zhemm(layout: c_int, side: c_int, uplo: c_int, m: isize, n: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, beta: [*c]const cf64, c: [*c]cf64, ldc: isize) void;
pub const chemm = cblas_chemm;
pub const zhemm = cblas_zhemm;

extern fn cblas_cherk(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: f32, a: [*c]const cf32, lda: isize, beta: f32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zherk(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: f64, a: [*c]const cf64, lda: isize, beta: f64, c: [*c]cf64, ldc: isize) void;
pub const cherk = cblas_cherk;
pub const zherk = cblas_zherk;

extern fn cblas_cher2k(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, beta: f32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zher2k(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, beta: f64, c: [*c]cf64, ldc: isize) void;
pub const cher2k = cblas_cher2k;
pub const zher2k = cblas_zher2k;

extern fn cblas_ssymm(layout: c_int, side: c_int, uplo: c_int, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, beta: f32, c: [*c]f32, ldc: isize) void;
extern fn cblas_dsymm(layout: c_int, side: c_int, uplo: c_int, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, beta: f64, c: [*c]f64, ldc: isize) void;
extern fn cblas_csymm(layout: c_int, side: c_int, uplo: c_int, m: isize, n: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, beta: [*c]const cf32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zsymm(layout: c_int, side: c_int, uplo: c_int, m: isize, n: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, beta: [*c]const cf64, c: [*c]cf64, ldc: isize) void;
pub const ssymm = cblas_ssymm;
pub const dsymm = cblas_dsymm;
pub const csymm = cblas_csymm;
pub const zsymm = cblas_zsymm;

extern fn cblas_ssyrk(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, beta: f32, c: [*c]f32, ldc: isize) void;
extern fn cblas_dsyrk(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, beta: f64, c: [*c]f64, ldc: isize) void;
extern fn cblas_csyrk(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, beta: [*c]const cf32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zsyrk(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, beta: [*c]const cf64, c: [*c]cf64, ldc: isize) void;
pub const ssyrk = cblas_ssyrk;
pub const dsyrk = cblas_dsyrk;
pub const csyrk = cblas_csyrk;
pub const zsyrk = cblas_zsyrk;

extern fn cblas_ssyr2k(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, beta: f32, c: [*c]f32, ldc: isize) void;
extern fn cblas_dsyr2k(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, beta: f64, c: [*c]f64, ldc: isize) void;
extern fn cblas_csyr2k(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, beta: [*c]const cf32, c: [*c]cf32, ldc: isize) void;
extern fn cblas_zsyr2k(layout: c_int, uplo: c_int, trans: c_int, n: isize, k: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, beta: [*c]const cf64, c: [*c]cf64, ldc: isize) void;
pub const ssyr2k = cblas_ssyr2k;
pub const dsyr2k = cblas_dsyr2k;
pub const csyr2k = cblas_csyr2k;
pub const zsyr2k = cblas_zsyr2k;

extern fn cblas_strmm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) void;
extern fn cblas_dtrmm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) void;
extern fn cblas_ctrmm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) void;
extern fn cblas_ztrmm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) void;
pub const strmm = cblas_strmm;
pub const dtrmm = cblas_dtrmm;
pub const ctrmm = cblas_ctrmm;
pub const ztrmm = cblas_ztrmm;

extern fn cblas_strsm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: f32, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) void;
extern fn cblas_dtrsm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: f64, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) void;
extern fn cblas_ctrsm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: [*c]const cf32, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) void;
extern fn cblas_ztrsm(layout: c_int, side: c_int, uplo: c_int, transa: c_int, diag: c_int, m: isize, n: isize, alpha: [*c]const cf64, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) void;
pub const strsm = cblas_strsm;
pub const dtrsm = cblas_dtrsm;
pub const ctrsm = cblas_ctrsm;
pub const ztrsm = cblas_ztrsm;
