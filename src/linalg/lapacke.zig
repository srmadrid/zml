const complex = @import("../complex.zig");
const cf32 = complex.cf32;
const cf64 = complex.cf64;

extern fn LAPACKE_sbdsdc(layout: c_int, uplo: u8, compq: u8, n: isize, d: [*c]f32, e: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, q: [*c]f32, iq: [*c]isize) isize;
extern fn LAPACKE_dbdsdc(layout: c_int, uplo: u8, compq: u8, n: isize, d: [*c]f64, e: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, q: [*c]f64, iq: [*c]isize) isize;
pub const sbdsdc = LAPACKE_sbdsdc;
pub const dbdsdc = LAPACKE_dbdsdc;

extern fn LAPACKE_sbdsqr(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f32, e: [*c]f32, vt: [*c]f32, ldvt: isize, u: [*c]f32, ldu: isize, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dbdsqr(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f64, e: [*c]f64, vt: [*c]f64, ldvt: isize, u: [*c]f64, ldu: isize, c: [*c]f64, ldc: isize) isize;
extern fn LAPACKE_cbdsqr(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f32, e: [*c]f32, vt: [*c]cf32, ldvt: isize, u: [*c]cf32, ldu: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zbdsqr(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f64, e: [*c]f64, vt: [*c]cf64, ldvt: isize, u: [*c]cf64, ldu: isize, c: [*c]cf64, ldc: isize) isize;
pub const sbdsqr = LAPACKE_sbdsqr;
pub const dbdsqr = LAPACKE_dbdsqr;
pub const cbdsqr = LAPACKE_cbdsqr;
pub const zbdsqr = LAPACKE_zbdsqr;

extern fn LAPACKE_sbdsvdx(layout: c_int, uplo: u8, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, ns: [*c]isize, s: [*c]f32, z: [*c]f32, ldz: isize, superb: [*c]isize) isize;
extern fn LAPACKE_dbdsvdx(layout: c_int, uplo: u8, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, ns: [*c]isize, s: [*c]f64, z: [*c]f64, ldz: isize, superb: [*c]isize) isize;
pub const sbdsvdx = LAPACKE_sbdsvdx;
pub const dbdsvdx = LAPACKE_dbdsvdx;

extern fn LAPACKE_sdisna(job: u8, m: isize, n: isize, d: [*c]const f32, sep: [*c]f32) isize;
extern fn LAPACKE_ddisna(job: u8, m: isize, n: isize, d: [*c]const f64, sep: [*c]f64) isize;
pub const sdisna = LAPACKE_sdisna;
pub const ddisna = LAPACKE_ddisna;

extern fn LAPACKE_sgbbrd(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]f32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]f32, ldq: isize, pt: [*c]f32, ldpt: isize, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dgbbrd(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]f64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]f64, ldq: isize, pt: [*c]f64, ldpt: isize, c: [*c]f64, ldc: isize) isize;
extern fn LAPACKE_cgbbrd(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]cf32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]cf32, ldq: isize, pt: [*c]cf32, ldpt: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zgbbrd(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]cf64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]cf64, ldq: isize, pt: [*c]cf64, ldpt: isize, c: [*c]cf64, ldc: isize) isize;
pub const sgbbrd = LAPACKE_sgbbrd;
pub const dgbbrd = LAPACKE_dgbbrd;
pub const cgbbrd = LAPACKE_cgbbrd;
pub const zgbbrd = LAPACKE_zgbbrd;

extern fn LAPACKE_sgbcon(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dgbcon(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cgbcon(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zgbcon(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const sgbcon = LAPACKE_sgbcon;
pub const dgbcon = LAPACKE_dgbcon;
pub const cgbcon = LAPACKE_cgbcon;
pub const zgbcon = LAPACKE_zgbcon;

extern fn LAPACKE_sgbequ(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgbequ(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgbequ(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgbequ(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgbequ = LAPACKE_sgbequ;
pub const dgbequ = LAPACKE_dgbequ;
pub const cgbequ = LAPACKE_cgbequ;
pub const zgbequ = LAPACKE_zgbequ;

extern fn LAPACKE_sgbequb(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgbequb(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgbequb(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgbequb(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgbequb = LAPACKE_sgbequb;
pub const dgbequb = LAPACKE_dgbequb;
pub const cgbequb = LAPACKE_cgbequb;
pub const zgbequb = LAPACKE_zgbequb;

extern fn LAPACKE_sgbrfs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, afb: [*c]const f32, ldafb: isize, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dgbrfs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, afb: [*c]const f64, ldafb: isize, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cgbrfs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, afb: [*c]const cf32, ldafb: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zgbrfs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, afb: [*c]const cf64, ldafb: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sgbrfs = LAPACKE_sgbrfs;
pub const dgbrfs = LAPACKE_dgbrfs;
pub const cgbrfs = LAPACKE_cgbrfs;
pub const zgbrfs = LAPACKE_zgbrfs;

extern fn LAPACKE_sgbrfsx(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, afb: [*c]const f32, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dgbrfsx(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, afb: [*c]const f64, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_cgbrfsx(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, afb: [*c]const cf32, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zgbrfsx(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, afb: [*c]const cf64, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const sgbrfsx = LAPACKE_sgbrfsx;
pub const dgbrfsx = LAPACKE_dgbrfsx;
pub const cgbrfsx = LAPACKE_cgbrfsx;
pub const zgbrfsx = LAPACKE_zgbrfsx;

extern fn LAPACKE_sgbsv(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f32, ldab: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgbsv(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f64, ldab: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgbsv(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgbsv(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const sgbsv = LAPACKE_sgbsv;
pub const dgbsv = LAPACKE_dgbsv;
pub const cgbsv = LAPACKE_cgbsv;
pub const zgbsv = LAPACKE_zgbsv;

extern fn LAPACKE_sgbsvx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f32, ldab: isize, afb: [*c]f32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, rpivot: [*c]f32) isize;
extern fn LAPACKE_dgbsvx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f64, ldab: isize, afb: [*c]f64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, rpivot: [*c]f64) isize;
extern fn LAPACKE_cgbsvx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, afb: [*c]cf32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, rpivot: [*c]f32) isize;
extern fn LAPACKE_zgbsvx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, afb: [*c]cf64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, rpivot: [*c]f64) isize;
pub const sgbsvx = LAPACKE_sgbsvx;
pub const dgbsvx = LAPACKE_dgbsvx;
pub const cgbsvx = LAPACKE_cgbsvx;
pub const zgbsvx = LAPACKE_zgbsvx;

extern fn LAPACKE_sgbsvxx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f32, ldab: isize, afb: [*c]f32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dgbsvxx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f64, ldab: isize, afb: [*c]f64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_cgbsvxx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, afb: [*c]cf32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zgbsvxx(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, afb: [*c]cf64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const sgbsvxx = LAPACKE_sgbsvxx;
pub const dgbsvxx = LAPACKE_dgbsvxx;
pub const cgbsvxx = LAPACKE_cgbsvxx;
pub const zgbsvxx = LAPACKE_zgbsvxx;

extern fn LAPACKE_sgbtrf(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]f32, ldab: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgbtrf(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]f64, ldab: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgbtrf(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]cf32, ldab: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgbtrf(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]cf64, ldab: isize, ipiv: [*c]isize) isize;
pub const sgbtrf = LAPACKE_sgbtrf;
pub const dgbtrf = LAPACKE_dgbtrf;
pub const cgbtrf = LAPACKE_cgbtrf;
pub const zgbtrf = LAPACKE_zgbtrf;

extern fn LAPACKE_sgbtrs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgbtrs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgbtrs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgbtrs(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const sgbtrs = LAPACKE_sgbtrs;
pub const dgbtrs = LAPACKE_dgbtrs;
pub const cgbtrs = LAPACKE_cgbtrs;
pub const zgbtrs = LAPACKE_zgbtrs;

extern fn LAPACKE_sgebak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f32, m: isize, v: [*c]f32, ldv: isize) isize;
extern fn LAPACKE_dgebak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f64, m: isize, v: [*c]f64, ldv: isize) isize;
extern fn LAPACKE_cgebak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f32, m: isize, v: [*c]cf32, ldv: isize) isize;
extern fn LAPACKE_zgebak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f64, m: isize, v: [*c]cf64, ldv: isize) isize;
pub const sgebak = LAPACKE_sgebak;
pub const dgebak = LAPACKE_dgebak;
pub const cgebak = LAPACKE_cgebak;
pub const zgebak = LAPACKE_zgebak;

extern fn LAPACKE_sgebal(layout: c_int, job: u8, n: isize, a: [*c]f32, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32) isize;
extern fn LAPACKE_dgebal(layout: c_int, job: u8, n: isize, a: [*c]f64, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64) isize;
extern fn LAPACKE_cgebal(layout: c_int, job: u8, n: isize, a: [*c]cf32, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32) isize;
extern fn LAPACKE_zgebal(layout: c_int, job: u8, n: isize, a: [*c]cf64, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64) isize;
pub const sgebal = LAPACKE_sgebal;
pub const dgebal = LAPACKE_dgebal;
pub const cgebal = LAPACKE_cgebal;
pub const zgebal = LAPACKE_zgebal;

extern fn LAPACKE_sgebrd(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, d: [*c]f32, e: [*c]f32, tauq: [*c]f32, taup: [*c]f32) isize;
extern fn LAPACKE_dgebrd(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, d: [*c]f64, e: [*c]f64, tauq: [*c]f64, taup: [*c]f64) isize;
extern fn LAPACKE_cgebrd(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, d: [*c]f32, e: [*c]f32, tauq: [*c]cf32, taup: [*c]cf32) isize;
extern fn LAPACKE_zgebrd(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, d: [*c]f64, e: [*c]f64, tauq: [*c]cf64, taup: [*c]cf64) isize;
pub const sgebrd = LAPACKE_sgebrd;
pub const dgebrd = LAPACKE_dgebrd;
pub const cgebrd = LAPACKE_cgebrd;
pub const zgebrd = LAPACKE_zgebrd;

extern fn LAPACKE_sgecon(layout: c_int, norm: u8, n: isize, a: [*c]const f32, lda: isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dgecon(layout: c_int, norm: u8, n: isize, a: [*c]const f64, lda: isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cgecon(layout: c_int, norm: u8, n: isize, a: [*c]const cf32, lda: isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zgecon(layout: c_int, norm: u8, n: isize, a: [*c]const cf64, lda: isize, anorm: f64, rcond: [*c]f64) isize;
pub const sgecon = LAPACKE_sgecon;
pub const dgecon = LAPACKE_dgecon;
pub const cgecon = LAPACKE_cgecon;
pub const zgecon = LAPACKE_zgecon;

extern fn LAPACKE_sgeequ(layout: c_int, m: isize, n: isize, a: [*c]const f32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgeequ(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgeequ(layout: c_int, m: isize, n: isize, a: [*c]const cf32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgeequ(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgeequ = LAPACKE_sgeequ;
pub const dgeequ = LAPACKE_dgeequ;
pub const cgeequ = LAPACKE_cgeequ;
pub const zgeequ = LAPACKE_zgeequ;

extern fn LAPACKE_sgeequb(layout: c_int, m: isize, n: isize, a: [*c]const f32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgeequb(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgeequb(layout: c_int, m: isize, n: isize, a: [*c]const cf32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgeequb(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgeequb = LAPACKE_sgeequb;
pub const dgeequb = LAPACKE_dgeequb;
pub const cgeequb = LAPACKE_cgeequb;
pub const zgeequb = LAPACKE_zgeequb;

extern fn LAPACKE_sgees(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f32, [*c]const f32) isize, n: isize, a: [*c]f32, lda: isize, sdim: [*c]isize, wr: [*c]f32, wi: [*c]f32, vs: [*c]f32, ldvs: isize) isize;
extern fn LAPACKE_dgees(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f64, [*c]const f64) isize, n: isize, a: [*c]f64, lda: isize, sdim: [*c]isize, wr: [*c]f64, wi: [*c]f64, vs: [*c]f64, ldvs: isize) isize;
extern fn LAPACKE_cgees(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf32) isize, n: isize, a: [*c]cf32, lda: isize, sdim: [*c]isize, w: [*c]cf32, vs: [*c]cf32, ldvs: isize) isize;
extern fn LAPACKE_zgees(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf64) isize, n: isize, a: [*c]cf64, lda: isize, sdim: [*c]isize, w: [*c]cf64, vs: [*c]cf64, ldvs: isize) isize;
pub const sgees = LAPACKE_sgees;
pub const dgees = LAPACKE_dgees;
pub const cgees = LAPACKE_cgees;
pub const zgees = LAPACKE_zgees;

extern fn LAPACKE_sgeesx(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f32, [*c]const f32) isize, sense: u8, n: isize, a: [*c]f32, lda: isize, sdim: [*c]isize, wr: [*c]f32, wi: [*c]f32, vs: [*c]f32, ldvs: isize, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_dgeesx(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f64, [*c]const f64) isize, sense: u8, n: isize, a: [*c]f64, lda: isize, sdim: [*c]isize, wr: [*c]f64, wi: [*c]f64, vs: [*c]f64, ldvs: isize, rconde: [*c]f64, rcondv: [*c]f64) isize;
extern fn LAPACKE_cgeesx(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf32) isize, sense: u8, n: isize, a: [*c]cf32, lda: isize, sdim: [*c]isize, w: [*c]cf32, vs: [*c]cf32, ldvs: isize, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_zgeesx(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf64) isize, sense: u8, n: isize, a: [*c]cf64, lda: isize, sdim: [*c]isize, w: [*c]cf64, vs: [*c]cf64, ldvs: isize, rconde: [*c]f64, rcondv: [*c]f64) isize;
pub const sgeesx = LAPACKE_sgeesx;
pub const dgeesx = LAPACKE_dgeesx;
pub const cgeesx = LAPACKE_cgeesx;
pub const zgeesx = LAPACKE_zgeesx;

extern fn LAPACKE_sgeev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f32, lda: isize, wr: [*c]f32, wi: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize) isize;
extern fn LAPACKE_dgeev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f64, lda: isize, wr: [*c]f64, wi: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize) isize;
extern fn LAPACKE_cgeev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize) isize;
extern fn LAPACKE_zgeev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize) isize;
pub const sgeev = LAPACKE_sgeev;
pub const dgeev = LAPACKE_dgeev;
pub const cgeev = LAPACKE_cgeev;
pub const zgeev = LAPACKE_zgeev;

extern fn LAPACKE_sgeevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f32, lda: isize, wr: [*c]f32, wi: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32, abnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_dgeevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f64, lda: isize, wr: [*c]f64, wi: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64, abnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64) isize;
extern fn LAPACKE_cgeevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32, abnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_zgeevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64, abnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64) isize;
pub const sgeevx = LAPACKE_sgeevx;
pub const dgeevx = LAPACKE_dgeevx;
pub const cgeevx = LAPACKE_cgeevx;
pub const zgeevx = LAPACKE_zgeevx;

extern fn LAPACKE_sgehrd(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgehrd(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgehrd(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgehrd(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgehrd = LAPACKE_sgehrd;
pub const dgehrd = LAPACKE_dgehrd;
pub const cgehrd = LAPACKE_cgehrd;
pub const zgehrd = LAPACKE_zgehrd;

extern fn LAPACKE_sgejsv(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]f32, lda: isize, sva: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, stat: [*c]f32, istat: [*c]isize) isize;
extern fn LAPACKE_dgejsv(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]f64, lda: isize, sva: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, stat: [*c]f64, istat: [*c]isize) isize;
extern fn LAPACKE_cgejsv(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, sva: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, stat: [*c]f32, istat: [*c]isize) isize;
extern fn LAPACKE_zgejsv(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, sva: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, stat: [*c]f64, istat: [*c]isize) isize;
pub const sgejsv = LAPACKE_sgejsv;
pub const dgejsv = LAPACKE_dgejsv;
pub const cgejsv = LAPACKE_cgejsv;
pub const zgejsv = LAPACKE_zgejsv;

extern fn LAPACKE_sgelq2(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgelq2(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgelq2(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgelq2(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgelq2 = LAPACKE_sgelq2;
pub const dgelq2 = LAPACKE_dgelq2;
pub const cgelq2 = LAPACKE_cgelq2;
pub const zgelq2 = LAPACKE_zgelq2;

extern fn LAPACKE_sgelqf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgelqf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgelqf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgelqf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgelqf = LAPACKE_sgelqf;
pub const dgelqf = LAPACKE_dgelqf;
pub const cgelqf = LAPACKE_cgelqf;
pub const zgelqf = LAPACKE_zgelqf;

extern fn LAPACKE_sgels(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgels(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgels(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgels(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const sgels = LAPACKE_sgels;
pub const dgels = LAPACKE_dgels;
pub const cgels = LAPACKE_cgels;
pub const zgels = LAPACKE_zgels;

extern fn LAPACKE_sgelsd(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize) isize;
extern fn LAPACKE_dgelsd(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize) isize;
extern fn LAPACKE_cgelsd(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize) isize;
extern fn LAPACKE_zgelsd(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize) isize;
pub const sgelsd = LAPACKE_sgelsd;
pub const dgelsd = LAPACKE_dgelsd;
pub const cgelsd = LAPACKE_cgelsd;
pub const zgelsd = LAPACKE_zgelsd;

extern fn LAPACKE_sgelss(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize) isize;
extern fn LAPACKE_dgelss(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize) isize;
extern fn LAPACKE_cgelss(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize) isize;
extern fn LAPACKE_zgelss(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize) isize;
pub const sgelss = LAPACKE_sgelss;
pub const dgelss = LAPACKE_dgelss;
pub const cgelss = LAPACKE_cgelss;
pub const zgelss = LAPACKE_zgelss;

extern fn LAPACKE_sgelsy(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, jpvt: [*c]isize, rcond: f32, rank: [*c]isize) isize;
extern fn LAPACKE_dgelsy(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, jpvt: [*c]isize, rcond: f64, rank: [*c]isize) isize;
extern fn LAPACKE_cgelsy(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, jpvt: [*c]isize, rcond: f32, rank: [*c]isize) isize;
extern fn LAPACKE_zgelsy(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, jpvt: [*c]isize, rcond: f64, rank: [*c]isize) isize;
pub const sgelsy = LAPACKE_sgelsy;
pub const dgelsy = LAPACKE_dgelsy;
pub const cgelsy = LAPACKE_cgelsy;
pub const zgelsy = LAPACKE_zgelsy;

extern fn LAPACKE_sgeqlf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgeqlf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgeqlf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgeqlf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgeqlf = LAPACKE_sgeqlf;
pub const dgeqlf = LAPACKE_dgeqlf;
pub const cgeqlf = LAPACKE_cgeqlf;
pub const zgeqlf = LAPACKE_zgeqlf;

extern fn LAPACKE_sgeqp3(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, jpvt: [*c]isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgeqp3(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, jpvt: [*c]isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgeqp3(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, jpvt: [*c]isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgeqp3(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, jpvt: [*c]isize, tau: [*c]cf64) isize;
pub const sgeqp3 = LAPACKE_sgeqp3;
pub const dgeqp3 = LAPACKE_dgeqp3;
pub const cgeqp3 = LAPACKE_cgeqp3;
pub const zgeqp3 = LAPACKE_zgeqp3;

extern fn LAPACKE_sgeqpf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, jpvt: [*c]isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgeqpf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, jpvt: [*c]isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgeqpf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, jpvt: [*c]isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgeqpf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, jpvt: [*c]isize, tau: [*c]cf64) isize;
pub const sgeqpf = LAPACKE_sgeqpf;
pub const dgeqpf = LAPACKE_dgeqpf;
pub const cgeqpf = LAPACKE_cgeqpf;
pub const zgeqpf = LAPACKE_zgeqpf;

extern fn LAPACKE_sgeqr2(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgeqr2(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgeqr2(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgeqr2(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgeqr2 = LAPACKE_sgeqr2;
pub const dgeqr2 = LAPACKE_dgeqr2;
pub const cgeqr2 = LAPACKE_cgeqr2;
pub const zgeqr2 = LAPACKE_zgeqr2;

extern fn LAPACKE_sgeqrf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgeqrf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgeqrf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgeqrf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgeqrf = LAPACKE_sgeqrf;
pub const dgeqrf = LAPACKE_dgeqrf;
pub const cgeqrf = LAPACKE_cgeqrf;
pub const zgeqrf = LAPACKE_zgeqrf;

extern fn LAPACKE_sgeqrfp(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgeqrfp(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgeqrfp(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgeqrfp(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgeqrfp = LAPACKE_sgeqrfp;
pub const dgeqrfp = LAPACKE_dgeqrfp;
pub const cgeqrfp = LAPACKE_cgeqrfp;
pub const zgeqrfp = LAPACKE_zgeqrfp;

extern fn LAPACKE_sgerfs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dgerfs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cgerfs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zgerfs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sgerfs = LAPACKE_sgerfs;
pub const dgerfs = LAPACKE_dgerfs;
pub const cgerfs = LAPACKE_cgerfs;
pub const zgerfs = LAPACKE_zgerfs;

extern fn LAPACKE_sgerfsx(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dgerfsx(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_cgerfsx(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zgerfsx(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const sgerfsx = LAPACKE_sgerfsx;
pub const dgerfsx = LAPACKE_dgerfsx;
pub const cgerfsx = LAPACKE_cgerfsx;
pub const zgerfsx = LAPACKE_zgerfsx;

extern fn LAPACKE_sgerqf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dgerqf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_cgerqf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zgerqf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const sgerqf = LAPACKE_sgerqf;
pub const dgerqf = LAPACKE_dgerqf;
pub const cgerqf = LAPACKE_cgerqf;
pub const zgerqf = LAPACKE_zgerqf;

extern fn LAPACKE_sgesdd(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]f32, lda: isize, s: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize) isize;
extern fn LAPACKE_dgesdd(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]f64, lda: isize, s: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize) isize;
extern fn LAPACKE_cgesdd(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, s: [*c]f32, u: [*c]cf32, ldu: isize, vt: [*c]cf32, ldvt: isize) isize;
extern fn LAPACKE_zgesdd(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, s: [*c]f64, u: [*c]cf64, ldu: isize, vt: [*c]cf64, ldvt: isize) isize;
pub const sgesdd = LAPACKE_sgesdd;
pub const dgesdd = LAPACKE_dgesdd;
pub const cgesdd = LAPACKE_cgesdd;
pub const zgesdd = LAPACKE_zgesdd;

extern fn LAPACKE_sgesv(layout: c_int, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgesv(layout: c_int, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgesv(layout: c_int, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgesv(layout: c_int, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
extern fn LAPACKE_dsgesv(layout: c_int, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, iter: [*c]isize) isize;
extern fn LAPACKE_zcgesv(layout: c_int, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, iter: [*c]isize) isize;
pub const sgesv = LAPACKE_sgesv;
pub const dgesv = LAPACKE_dgesv;
pub const cgesv = LAPACKE_cgesv;
pub const zgesv = LAPACKE_zgesv;
pub const dsgesv = LAPACKE_dsgesv;
pub const zcgesv = LAPACKE_zcgesv;

extern fn LAPACKE_sgesvd(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]f32, lda: isize, s: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, superb: [*c]f32) isize;
extern fn LAPACKE_dgesvd(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]f64, lda: isize, s: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, superb: [*c]f64) isize;
extern fn LAPACKE_cgesvd(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, s: [*c]f32, u: [*c]cf32, ldu: isize, vt: [*c]cf32, ldvt: isize, superb: [*c]f32) isize;
extern fn LAPACKE_zgesvd(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, s: [*c]f64, u: [*c]cf64, ldu: isize, vt: [*c]cf64, ldvt: isize, superb: [*c]f64) isize;
pub const sgesvd = LAPACKE_sgesvd;
pub const dgesvd = LAPACKE_dgesvd;
pub const cgesvd = LAPACKE_cgesvd;
pub const zgesvd = LAPACKE_zgesvd;

extern fn LAPACKE_sgesvdx(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, ns: [*c]isize, s: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, superb: [*c]isize) isize;
extern fn LAPACKE_dgesvdx(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, ns: [*c]isize, s: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, superb: [*c]isize) isize;
extern fn LAPACKE_cgesvdx(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, ns: [*c]isize, s: [*c]f32, u: [*c]cf32, ldu: isize, vt: [*c]cf32, ldvt: isize, superb: [*c]isize) isize;
extern fn LAPACKE_zgesvdx(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, ns: [*c]isize, s: [*c]f64, u: [*c]cf64, ldu: isize, vt: [*c]cf64, ldvt: isize, superb: [*c]isize) isize;
pub const sgesvdx = LAPACKE_sgesvdx;
pub const dgesvdx = LAPACKE_dgesvdx;
pub const cgesvdx = LAPACKE_cgesvdx;
pub const zgesvdx = LAPACKE_zgesvdx;

extern fn LAPACKE_sgesvdq(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f32, lda: isize, s: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, numrank: [*c]isize) isize;
extern fn LAPACKE_dgesvdq(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f64, lda: isize, s: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, numrank: [*c]isize) isize;
extern fn LAPACKE_cgesvdq(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, s: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, numrank: [*c]isize) isize;
extern fn LAPACKE_zgesvdq(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, s: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, numrank: [*c]isize) isize;
pub const sgesvdq = LAPACKE_sgesvdq;
pub const dgesvdq = LAPACKE_dgesvdq;
pub const cgesvdq = LAPACKE_cgesvdq;
pub const zgesvdq = LAPACKE_zgesvdq;

extern fn LAPACKE_sgesvj(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f32, lda: isize, sva: [*c]f32, mv: isize, v: [*c]f32, ldv: isize, stat: [*c]f32) isize;
extern fn LAPACKE_dgesvj(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f64, lda: isize, sva: [*c]f64, mv: isize, v: [*c]f64, ldv: isize, stat: [*c]f64) isize;
extern fn LAPACKE_cgesvj(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, sva: [*c]f32, mv: isize, v: [*c]cf32, ldv: isize, stat: [*c]f32) isize;
extern fn LAPACKE_zgesvj(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, sva: [*c]f64, mv: isize, v: [*c]cf64, ldv: isize, stat: [*c]f64) isize;
pub const sgesvj = LAPACKE_sgesvj;
pub const dgesvj = LAPACKE_dgesvj;
pub const cgesvj = LAPACKE_cgesvj;
pub const zgesvj = LAPACKE_zgesvj;

extern fn LAPACKE_sgesvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, rpivot: [*c]f32) isize;
extern fn LAPACKE_dgesvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, rpivot: [*c]f64) isize;
extern fn LAPACKE_cgesvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, rpivot: [*c]f32) isize;
extern fn LAPACKE_zgesvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, rpivot: [*c]f64) isize;
pub const sgesvx = LAPACKE_sgesvx;
pub const dgesvx = LAPACKE_dgesvx;
pub const cgesvx = LAPACKE_cgesvx;
pub const zgesvx = LAPACKE_zgesvx;

extern fn LAPACKE_sgesvxx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dgesvxx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_cgesvxx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zgesvxx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const sgesvxx = LAPACKE_sgesvxx;
pub const dgesvxx = LAPACKE_dgesvxx;
pub const cgesvxx = LAPACKE_cgesvxx;
pub const zgesvxx = LAPACKE_zgesvxx;

extern fn LAPACKE_sgetf2(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgetf2(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgetf2(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgetf2(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const sgetf2 = LAPACKE_sgetf2;
pub const dgetf2 = LAPACKE_dgetf2;
pub const cgetf2 = LAPACKE_cgetf2;
pub const zgetf2 = LAPACKE_zgetf2;

extern fn LAPACKE_sgetrf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgetrf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgetrf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgetrf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const sgetrf = LAPACKE_sgetrf;
pub const dgetrf = LAPACKE_dgetrf;
pub const cgetrf = LAPACKE_cgetrf;
pub const zgetrf = LAPACKE_zgetrf;

extern fn LAPACKE_sgetrf2(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgetrf2(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgetrf2(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgetrf2(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const sgetrf2 = LAPACKE_sgetrf2;
pub const dgetrf2 = LAPACKE_dgetrf2;
pub const cgetrf2 = LAPACKE_cgetrf2;
pub const zgetrf2 = LAPACKE_zgetrf2;

extern fn LAPACKE_sgetri(layout: c_int, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_dgetri(layout: c_int, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_cgetri(layout: c_int, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zgetri(layout: c_int, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize) isize;
pub const sgetri = LAPACKE_sgetri;
pub const dgetri = LAPACKE_dgetri;
pub const cgetri = LAPACKE_cgetri;
pub const zgetri = LAPACKE_zgetri;

extern fn LAPACKE_sgetrs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgetrs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgetrs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgetrs(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const sgetrs = LAPACKE_sgetrs;
pub const dgetrs = LAPACKE_dgetrs;
pub const cgetrs = LAPACKE_cgetrs;
pub const zgetrs = LAPACKE_zgetrs;

extern fn LAPACKE_sggbak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f32, rscale: [*c]const f32, m: isize, v: [*c]f32, ldv: isize) isize;
extern fn LAPACKE_dggbak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f64, rscale: [*c]const f64, m: isize, v: [*c]f64, ldv: isize) isize;
extern fn LAPACKE_cggbak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f32, rscale: [*c]const f32, m: isize, v: [*c]cf32, ldv: isize) isize;
extern fn LAPACKE_zggbak(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f64, rscale: [*c]const f64, m: isize, v: [*c]cf64, ldv: isize) isize;
pub const sggbak = LAPACKE_sggbak;
pub const dggbak = LAPACKE_dggbak;
pub const cggbak = LAPACKE_cggbak;
pub const zggbak = LAPACKE_zggbak;

extern fn LAPACKE_sggbal(layout: c_int, job: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32) isize;
extern fn LAPACKE_dggbal(layout: c_int, job: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64) isize;
extern fn LAPACKE_cggbal(layout: c_int, job: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32) isize;
extern fn LAPACKE_zggbal(layout: c_int, job: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64) isize;
pub const sggbal = LAPACKE_sggbal;
pub const dggbal = LAPACKE_dggbal;
pub const cggbal = LAPACKE_cggbal;
pub const zggbal = LAPACKE_zggbal;

extern fn LAPACKE_sgges(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f32, [*c]const f32, [*c]const f32) isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, sdim: [*c]isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vsl: [*c]f32, ldvsl: isize, vsr: [*c]f32, ldvsr: isize) isize;
extern fn LAPACKE_dgges(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f64, [*c]const f64, [*c]const f64) isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, sdim: [*c]isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vsl: [*c]f64, ldvsl: isize, vsr: [*c]f64, ldvsr: isize) isize;
extern fn LAPACKE_cgges(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf32, [*c]const cf32) isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, sdim: [*c]isize, alpha: [*c]cf32, beta: [*c]cf32, vsl: [*c]cf32, ldvsl: isize, vsr: [*c]cf32, ldvsr: isize) isize;
extern fn LAPACKE_zgges(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf64, [*c]const cf64) isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, sdim: [*c]isize, alpha: [*c]cf64, beta: [*c]cf64, vsl: [*c]cf64, ldvsl: isize, vsr: [*c]cf64, ldvsr: isize) isize;
pub const sgges = LAPACKE_sgges;
pub const dgges = LAPACKE_dgges;
pub const cgges = LAPACKE_cgges;
pub const zgges = LAPACKE_zgges;

extern fn LAPACKE_sgges3(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f32, [*c]const f32, [*c]const f32) isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, sdim: [*c]isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vsl: [*c]f32, ldvsl: isize, vsr: [*c]f32, ldvsr: isize) isize;
extern fn LAPACKE_dgges3(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f64, [*c]const f64, [*c]const f64) isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, sdim: [*c]isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vsl: [*c]f64, ldvsl: isize, vsr: [*c]f64, ldvsr: isize) isize;
extern fn LAPACKE_cgges3(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf32, [*c]const cf32) isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, sdim: [*c]isize, alpha: [*c]cf32, beta: [*c]cf32, vsl: [*c]cf32, ldvsl: isize, vsr: [*c]cf32, ldvsr: isize) isize;
extern fn LAPACKE_zgges3(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf64, [*c]const cf64) isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, sdim: [*c]isize, alpha: [*c]cf64, beta: [*c]cf64, vsl: [*c]cf64, ldvsl: isize, vsr: [*c]cf64, ldvsr: isize) isize;
pub const sgges3 = LAPACKE_sgges3;
pub const dgges3 = LAPACKE_dgges3;
pub const cgges3 = LAPACKE_cgges3;
pub const zgges3 = LAPACKE_zgges3;

extern fn LAPACKE_sggesx(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f32, [*c]const f32, [*c]const f32) isize, sense: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, sdim: [*c]isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vsl: [*c]f32, ldvsl: isize, vsr: [*c]f32, ldvsr: isize, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_dggesx(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f64, [*c]const f64, [*c]const f64) isize, sense: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, sdim: [*c]isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vsl: [*c]f64, ldvsl: isize, vsr: [*c]f64, ldvsr: isize, rconde: [*c]f64, rcondv: [*c]f64) isize;
extern fn LAPACKE_cggesx(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf32, [*c]const cf32) isize, sense: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, sdim: [*c]isize, alpha: [*c]cf32, beta: [*c]cf32, vsl: [*c]cf32, ldvsl: isize, vsr: [*c]cf32, ldvsr: isize, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_zggesx(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf64, [*c]const cf64) isize, sense: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, sdim: [*c]isize, alpha: [*c]cf64, beta: [*c]cf64, vsl: [*c]cf64, ldvsl: isize, vsr: [*c]cf64, ldvsr: isize, rconde: [*c]f64, rcondv: [*c]f64) isize;
pub const sggesx = LAPACKE_sggesx;
pub const dggesx = LAPACKE_dggesx;
pub const cggesx = LAPACKE_cggesx;
pub const zggesx = LAPACKE_zggesx;

extern fn LAPACKE_sggev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize) isize;
extern fn LAPACKE_dggev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize) isize;
extern fn LAPACKE_cggev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize) isize;
extern fn LAPACKE_zggev(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize) isize;
pub const sggev = LAPACKE_sggev;
pub const dggev = LAPACKE_dggev;
pub const cggev = LAPACKE_cggev;
pub const zggev = LAPACKE_zggev;

extern fn LAPACKE_sggev3(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize) isize;
extern fn LAPACKE_dggev3(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize) isize;
extern fn LAPACKE_cggev3(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize) isize;
extern fn LAPACKE_zggev3(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize) isize;
pub const sggev3 = LAPACKE_sggev3;
pub const dggev3 = LAPACKE_dggev3;
pub const cggev3 = LAPACKE_cggev3;
pub const zggev3 = LAPACKE_zggev3;

extern fn LAPACKE_sggevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32, abnrm: [*c]f32, bbnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_dggevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64, abnrm: [*c]f64, bbnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64) isize;
extern fn LAPACKE_cggevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32, abnrm: [*c]f32, bbnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32) isize;
extern fn LAPACKE_zggevx(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64, abnrm: [*c]f64, bbnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64) isize;
pub const sggevx = LAPACKE_sggevx;
pub const dggevx = LAPACKE_dggevx;
pub const cggevx = LAPACKE_cggevx;
pub const zggevx = LAPACKE_zggevx;

extern fn LAPACKE_sggglm(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, d: [*c]f32, x: [*c]f32, y: [*c]f32) isize;
extern fn LAPACKE_dggglm(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, d: [*c]f64, x: [*c]f64, y: [*c]f64) isize;
extern fn LAPACKE_cggglm(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, d: [*c]cf32, x: [*c]cf32, y: [*c]cf32) isize;
extern fn LAPACKE_zggglm(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, d: [*c]cf64, x: [*c]cf64, y: [*c]cf64) isize;
pub const sggglm = LAPACKE_sggglm;
pub const dggglm = LAPACKE_dggglm;
pub const cggglm = LAPACKE_cggglm;
pub const zggglm = LAPACKE_zggglm;

extern fn LAPACKE_sgghrd(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dgghrd(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_cgghrd(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zgghrd(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize) isize;
pub const sgghrd = LAPACKE_sgghrd;
pub const dgghrd = LAPACKE_dgghrd;
pub const cgghrd = LAPACKE_cgghrd;
pub const zgghrd = LAPACKE_zgghrd;

extern fn LAPACKE_sgghd3(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dgghd3(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_cgghd3(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zgghd3(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize) isize;
pub const sgghd3 = LAPACKE_sgghd3;
pub const dgghd3 = LAPACKE_dgghd3;
pub const cgghd3 = LAPACKE_cgghd3;
pub const zgghd3 = LAPACKE_zgghd3;

extern fn LAPACKE_sgglse(layout: c_int, m: isize, n: isize, p: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, c: [*c]f32, d: [*c]f32, x: [*c]f32) isize;
extern fn LAPACKE_dgglse(layout: c_int, m: isize, n: isize, p: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, c: [*c]f64, d: [*c]f64, x: [*c]f64) isize;
extern fn LAPACKE_cgglse(layout: c_int, m: isize, n: isize, p: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, c: [*c]cf32, d: [*c]cf32, x: [*c]cf32) isize;
extern fn LAPACKE_zgglse(layout: c_int, m: isize, n: isize, p: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, c: [*c]cf64, d: [*c]cf64, x: [*c]cf64) isize;
pub const sgglse = LAPACKE_sgglse;
pub const dgglse = LAPACKE_dgglse;
pub const cgglse = LAPACKE_cgglse;
pub const zgglse = LAPACKE_zgglse;

extern fn LAPACKE_sggqrf(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f32, lda: isize, taua: [*c]f32, b: [*c]f32, ldb: isize, taub: [*c]f32) isize;
extern fn LAPACKE_dggqrf(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f64, lda: isize, taua: [*c]f64, b: [*c]f64, ldb: isize, taub: [*c]f64) isize;
extern fn LAPACKE_cggqrf(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf32, lda: isize, taua: [*c]cf32, b: [*c]cf32, ldb: isize, taub: [*c]cf32) isize;
extern fn LAPACKE_zggqrf(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf64, lda: isize, taua: [*c]cf64, b: [*c]cf64, ldb: isize, taub: [*c]cf64) isize;
pub const sggqrf = LAPACKE_sggqrf;
pub const dggqrf = LAPACKE_dggqrf;
pub const cggqrf = LAPACKE_cggqrf;
pub const zggqrf = LAPACKE_zggqrf;

extern fn LAPACKE_sggrqf(layout: c_int, m: isize, p: isize, n: isize, a: [*c]f32, lda: isize, taua: [*c]f32, b: [*c]f32, ldb: isize, taub: [*c]f32) isize;
extern fn LAPACKE_dggrqf(layout: c_int, m: isize, p: isize, n: isize, a: [*c]f64, lda: isize, taua: [*c]f64, b: [*c]f64, ldb: isize, taub: [*c]f64) isize;
extern fn LAPACKE_cggrqf(layout: c_int, m: isize, p: isize, n: isize, a: [*c]cf32, lda: isize, taua: [*c]cf32, b: [*c]cf32, ldb: isize, taub: [*c]cf32) isize;
extern fn LAPACKE_zggrqf(layout: c_int, m: isize, p: isize, n: isize, a: [*c]cf64, lda: isize, taua: [*c]cf64, b: [*c]cf64, ldb: isize, taub: [*c]cf64) isize;
pub const sggrqf = LAPACKE_sggrqf;
pub const dggrqf = LAPACKE_dggrqf;
pub const cggrqf = LAPACKE_cggrqf;
pub const zggrqf = LAPACKE_zggrqf;

extern fn LAPACKE_sggsvd(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dggsvd(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cggsvd(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_zggsvd(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, iwork: [*c]isize) isize;
pub const sggsvd = LAPACKE_sggsvd;
pub const dggsvd = LAPACKE_dggsvd;
pub const cggsvd = LAPACKE_cggsvd;
pub const zggsvd = LAPACKE_zggsvd;

extern fn LAPACKE_sggsvd3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dggsvd3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cggsvd3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_zggsvd3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, iwork: [*c]isize) isize;
pub const sggsvd3 = LAPACKE_sggsvd3;
pub const dggsvd3 = LAPACKE_dggsvd3;
pub const cggsvd3 = LAPACKE_cggsvd3;
pub const zggsvd3 = LAPACKE_zggsvd3;

extern fn LAPACKE_sggsvp(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize) isize;
extern fn LAPACKE_dggsvp(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize) isize;
extern fn LAPACKE_cggsvp(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize) isize;
extern fn LAPACKE_zggsvp(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize) isize;
pub const sggsvp = LAPACKE_sggsvp;
pub const dggsvp = LAPACKE_dggsvp;
pub const cggsvp = LAPACKE_cggsvp;
pub const zggsvp = LAPACKE_zggsvp;

extern fn LAPACKE_sggsvp3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize) isize;
extern fn LAPACKE_dggsvp3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize) isize;
extern fn LAPACKE_cggsvp3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize) isize;
extern fn LAPACKE_zggsvp3(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize) isize;
pub const sggsvp3 = LAPACKE_sggsvp3;
pub const dggsvp3 = LAPACKE_dggsvp3;
pub const cggsvp3 = LAPACKE_cggsvp3;
pub const zggsvp3 = LAPACKE_zggsvp3;

extern fn LAPACKE_sgtcon(norm: u8, n: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, du2: [*c]const f32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dgtcon(norm: u8, n: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, du2: [*c]const f64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cgtcon(norm: u8, n: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, du2: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zgtcon(norm: u8, n: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, du2: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const sgtcon = LAPACKE_sgtcon;
pub const dgtcon = LAPACKE_dgtcon;
pub const cgtcon = LAPACKE_cgtcon;
pub const zgtcon = LAPACKE_zgtcon;

extern fn LAPACKE_sgtrfs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, dlf: [*c]const f32, df: [*c]const f32, duf: [*c]const f32, du2: [*c]const f32, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dgtrfs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, dlf: [*c]const f64, df: [*c]const f64, duf: [*c]const f64, du2: [*c]const f64, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cgtrfs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, dlf: [*c]const cf32, df: [*c]const cf32, duf: [*c]const cf32, du2: [*c]const cf32, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zgtrfs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, dlf: [*c]const cf64, df: [*c]const cf64, duf: [*c]const cf64, du2: [*c]const cf64, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sgtrfs = LAPACKE_sgtrfs;
pub const dgtrfs = LAPACKE_dgtrfs;
pub const cgtrfs = LAPACKE_cgtrfs;
pub const zgtrfs = LAPACKE_zgtrfs;

extern fn LAPACKE_sgtsv(layout: c_int, n: isize, nrhs: isize, dl: [*c]f32, d: [*c]f32, du: [*c]f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgtsv(layout: c_int, n: isize, nrhs: isize, dl: [*c]f64, d: [*c]f64, du: [*c]f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgtsv(layout: c_int, n: isize, nrhs: isize, dl: [*c]cf32, d: [*c]cf32, du: [*c]cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgtsv(layout: c_int, n: isize, nrhs: isize, dl: [*c]cf64, d: [*c]cf64, du: [*c]cf64, b: [*c]cf64, ldb: isize) isize;
pub const sgtsv = LAPACKE_sgtsv;
pub const dgtsv = LAPACKE_dgtsv;
pub const cgtsv = LAPACKE_cgtsv;
pub const zgtsv = LAPACKE_zgtsv;

extern fn LAPACKE_sgtsvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, dlf: [*c]f32, df: [*c]f32, duf: [*c]f32, du2: [*c]f32, ipiv: [*c]isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dgtsvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, dlf: [*c]f64, df: [*c]f64, duf: [*c]f64, du2: [*c]f64, ipiv: [*c]isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cgtsvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, dlf: [*c]cf32, df: [*c]cf32, duf: [*c]cf32, du2: [*c]cf32, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zgtsvx(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, dlf: [*c]cf64, df: [*c]cf64, duf: [*c]cf64, du2: [*c]cf64, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sgtsvx = LAPACKE_sgtsvx;
pub const dgtsvx = LAPACKE_dgtsvx;
pub const cgtsvx = LAPACKE_cgtsvx;
pub const zgtsvx = LAPACKE_zgtsvx;

extern fn LAPACKE_sgttrf(n: isize, dl: [*c]f32, d: [*c]f32, du: [*c]f32, du2: [*c]f32, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgttrf(n: isize, dl: [*c]f64, d: [*c]f64, du: [*c]f64, du2: [*c]f64, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgttrf(n: isize, dl: [*c]cf32, d: [*c]cf32, du: [*c]cf32, du2: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgttrf(n: isize, dl: [*c]cf64, d: [*c]cf64, du: [*c]cf64, du2: [*c]cf64, ipiv: [*c]isize) isize;
pub const sgttrf = LAPACKE_sgttrf;
pub const dgttrf = LAPACKE_dgttrf;
pub const cgttrf = LAPACKE_cgttrf;
pub const zgttrf = LAPACKE_zgttrf;

extern fn LAPACKE_sgttrs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, du2: [*c]const f32, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgttrs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, du2: [*c]const f64, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgttrs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, du2: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgttrs(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, du2: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const sgttrs = LAPACKE_sgttrs;
pub const dgttrs = LAPACKE_dgttrs;
pub const cgttrs = LAPACKE_cgttrs;
pub const zgttrs = LAPACKE_zgttrs;

extern fn LAPACKE_chbev(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhbev(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chbev = LAPACKE_chbev;
pub const zhbev = LAPACKE_zhbev;

extern fn LAPACKE_chbevd(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhbevd(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chbevd = LAPACKE_chbevd;
pub const zhbevd = LAPACKE_zhbevd;

extern fn LAPACKE_chbevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, q: [*c]cf32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhbevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, q: [*c]cf64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const chbevx = LAPACKE_chbevx;
pub const zhbevx = LAPACKE_zhbevx;

extern fn LAPACKE_chbgst(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]const cf32, ldbb: isize, x: [*c]cf32, ldx: isize) isize;
extern fn LAPACKE_zhbgst(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]const cf64, ldbb: isize, x: [*c]cf64, ldx: isize) isize;
pub const chbgst = LAPACKE_chbgst;
pub const zhbgst = LAPACKE_zhbgst;

extern fn LAPACKE_chbgv(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]cf32, ldbb: isize, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhbgv(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]cf64, ldbb: isize, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chbgv = LAPACKE_chbgv;
pub const zhbgv = LAPACKE_zhbgv;

extern fn LAPACKE_chbgvd(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]cf32, ldbb: isize, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhbgvd(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]cf64, ldbb: isize, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chbgvd = LAPACKE_chbgvd;
pub const zhbgvd = LAPACKE_zhbgvd;

extern fn LAPACKE_chbgvx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]cf32, ldbb: isize, q: [*c]cf32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhbgvx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]cf64, ldbb: isize, q: [*c]cf64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const chbgvx = LAPACKE_chbgvx;
pub const zhbgvx = LAPACKE_zhbgvx;

extern fn LAPACKE_chbtrd(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]cf32, ldq: isize) isize;
extern fn LAPACKE_zhbtrd(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]cf64, ldq: isize) isize;
pub const chbtrd = LAPACKE_chbtrd;
pub const zhbtrd = LAPACKE_zhbtrd;

extern fn LAPACKE_checon(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zhecon(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const checon = LAPACKE_checon;
pub const zhecon = LAPACKE_zhecon;

extern fn LAPACKE_cheequb(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zheequb(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const cheequb = LAPACKE_cheequb;
pub const zheequb = LAPACKE_zheequb;

extern fn LAPACKE_cheev(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_zheev(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64) isize;
pub const cheev = LAPACKE_cheev;
pub const zheev = LAPACKE_zheev;

extern fn LAPACKE_cheevd(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_zheevd(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64) isize;
pub const cheevd = LAPACKE_cheevd;
pub const zheevd = LAPACKE_zheevd;

extern fn LAPACKE_cheevr(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_zheevr(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, isuppz: [*c]isize) isize;
pub const cheevr = LAPACKE_cheevr;
pub const zheevr = LAPACKE_zheevr;

extern fn LAPACKE_cheevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zheevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const cheevx = LAPACKE_cheevx;
pub const zheevx = LAPACKE_zheevx;

extern fn LAPACKE_chegst(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]const cf32, ldb: isize) isize;
extern fn LAPACKE_zhegst(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]const cf64, ldb: isize) isize;
pub const chegst = LAPACKE_chegst;
pub const zhegst = LAPACKE_zhegst;

extern fn LAPACKE_chegv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, w: [*c]f32) isize;
extern fn LAPACKE_zhegv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, w: [*c]f64) isize;
pub const chegv = LAPACKE_chegv;
pub const zhegv = LAPACKE_zhegv;

extern fn LAPACKE_chegvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, w: [*c]f32) isize;
extern fn LAPACKE_zhegvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, w: [*c]f64) isize;
pub const chegvd = LAPACKE_chegvd;
pub const zhegvd = LAPACKE_zhegvd;

extern fn LAPACKE_chegvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhegvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const chegvx = LAPACKE_chegvx;
pub const zhegvx = LAPACKE_zhegvx;

extern fn LAPACKE_cherfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zherfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const cherfs = LAPACKE_cherfs;
pub const zherfs = LAPACKE_zherfs;

extern fn LAPACKE_cherfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zherfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const cherfsx = LAPACKE_cherfsx;
pub const zherfsx = LAPACKE_zherfsx;

extern fn LAPACKE_chesv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhesv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const chesv = LAPACKE_chesv;
pub const zhesv = LAPACKE_zhesv;

extern fn LAPACKE_chesvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zhesvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const chesvx = LAPACKE_chesvx;
pub const zhesvx = LAPACKE_zhesvx;

extern fn LAPACKE_chesvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zhesvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const chesvxx = LAPACKE_chesvxx;
pub const zhesvxx = LAPACKE_zhesvxx;

extern fn LAPACKE_chetrd(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, d: [*c]f32, e: [*c]f32, tau: [*c]cf32) isize;
extern fn LAPACKE_zhetrd(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, d: [*c]f64, e: [*c]f64, tau: [*c]cf64) isize;
pub const chetrd = LAPACKE_chetrd;
pub const zhetrd = LAPACKE_zhetrd;

extern fn LAPACKE_chetrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zhetrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const chetrf = LAPACKE_chetrf;
pub const zhetrf = LAPACKE_zhetrf;

extern fn LAPACKE_chetri(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zhetri(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize) isize;
pub const chetri = LAPACKE_chetri;
pub const zhetri = LAPACKE_zhetri;

extern fn LAPACKE_chetrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhetrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const chetrs = LAPACKE_chetrs;
pub const zhetrs = LAPACKE_zhetrs;

extern fn LAPACKE_chfrk(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f32, a: [*c]const cf32, lda: isize, beta: f32, c: [*c]cf32) isize;
extern fn LAPACKE_zhfrk(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f64, a: [*c]const cf64, lda: isize, beta: f64, c: [*c]cf64) isize;
pub const chfrk = LAPACKE_chfrk;
pub const zhfrk = LAPACKE_zhfrk;

extern fn LAPACKE_shgeqz(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f32, ldh: isize, t: [*c]f32, ldt: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dhgeqz(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f64, ldh: isize, t: [*c]f64, ldt: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_chgeqz(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf32, ldh: isize, t: [*c]cf32, ldt: isize, alpha: [*c]cf32, beta: [*c]cf32, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhgeqz(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf64, ldh: isize, t: [*c]cf64, ldt: isize, alpha: [*c]cf64, beta: [*c]cf64, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize) isize;
pub const shgeqz = LAPACKE_shgeqz;
pub const dhgeqz = LAPACKE_dhgeqz;
pub const chgeqz = LAPACKE_chgeqz;
pub const zhgeqz = LAPACKE_zhgeqz;

extern fn LAPACKE_chpcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zhpcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const chpcon = LAPACKE_chpcon;
pub const zhpcon = LAPACKE_zhpcon;

extern fn LAPACKE_chpev(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhpev(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chpev = LAPACKE_chpev;
pub const zhpev = LAPACKE_zhpev;

extern fn LAPACKE_chpevd(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhpevd(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chpevd = LAPACKE_chpevd;
pub const zhpevd = LAPACKE_zhpevd;

extern fn LAPACKE_chpevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhpevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const chpevx = LAPACKE_chpevx;
pub const zhpevx = LAPACKE_zhpevx;

extern fn LAPACKE_chpgst(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]const cf32) isize;
extern fn LAPACKE_zhpgst(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]const cf64) isize;
pub const chpgst = LAPACKE_chpgst;
pub const zhpgst = LAPACKE_zhpgst;

extern fn LAPACKE_chpgv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhpgv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chpgv = LAPACKE_chpgv;
pub const zhpgv = LAPACKE_zhpgv;

extern fn LAPACKE_chpgvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhpgvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chpgvd = LAPACKE_chpgvd;
pub const zhpgvd = LAPACKE_zhpgvd;

extern fn LAPACKE_chpgvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]cf32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhpgvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]cf64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const chpgvx = LAPACKE_chpgvx;
pub const zhpgvx = LAPACKE_zhpgvx;

extern fn LAPACKE_chprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]const cf32, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zhprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]const cf64, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const chprfs = LAPACKE_chprfs;
pub const zhprfs = LAPACKE_zhprfs;

extern fn LAPACKE_chpsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhpsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const chpsv = LAPACKE_chpsv;
pub const zhpsv = LAPACKE_zhpsv;

extern fn LAPACKE_chpsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]cf32, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zhpsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]cf64, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const chpsvx = LAPACKE_chpsvx;
pub const zhpsvx = LAPACKE_zhpsvx;

extern fn LAPACKE_chptrd(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, d: [*c]f32, e: [*c]f32, tau: [*c]cf32) isize;
extern fn LAPACKE_zhptrd(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, d: [*c]f64, e: [*c]f64, tau: [*c]cf64) isize;
pub const chptrd = LAPACKE_chptrd;
pub const zhptrd = LAPACKE_zhptrd;

extern fn LAPACKE_chptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zhptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]isize) isize;
pub const chptrf = LAPACKE_chptrf;
pub const zhptrf = LAPACKE_zhptrf;

extern fn LAPACKE_chptri(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zhptri(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]const isize) isize;
pub const chptri = LAPACKE_chptri;
pub const zhptri = LAPACKE_zhptri;

extern fn LAPACKE_chptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const chptrs = LAPACKE_chptrs;
pub const zhptrs = LAPACKE_zhptrs;

extern fn LAPACKE_shsein(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]isize, n: isize, h: [*c]const f32, ldh: isize, wr: [*c]f32, wi: [*c]const f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, mm: isize, m: [*c]isize, ifaill: [*c]isize, ifailr: [*c]isize) isize;
extern fn LAPACKE_dhsein(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]isize, n: isize, h: [*c]const f64, ldh: isize, wr: [*c]f64, wi: [*c]const f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, mm: isize, m: [*c]isize, ifaill: [*c]isize, ifailr: [*c]isize) isize;
extern fn LAPACKE_chsein(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]const isize, n: isize, h: [*c]const cf32, ldh: isize, w: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, mm: isize, m: [*c]isize, ifaill: [*c]isize, ifailr: [*c]isize) isize;
extern fn LAPACKE_zhsein(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]const isize, n: isize, h: [*c]const cf64, ldh: isize, w: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, mm: isize, m: [*c]isize, ifaill: [*c]isize, ifailr: [*c]isize) isize;
pub const shsein = LAPACKE_shsein;
pub const dhsein = LAPACKE_dhsein;
pub const chsein = LAPACKE_chsein;
pub const zhsein = LAPACKE_zhsein;

extern fn LAPACKE_shseqr(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f32, ldh: isize, wr: [*c]f32, wi: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dhseqr(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f64, ldh: isize, wr: [*c]f64, wi: [*c]f64, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_chseqr(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf32, ldh: isize, w: [*c]cf32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhseqr(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf64, ldh: isize, w: [*c]cf64, z: [*c]cf64, ldz: isize) isize;
pub const shseqr = LAPACKE_shseqr;
pub const dhseqr = LAPACKE_dhseqr;
pub const chseqr = LAPACKE_chseqr;
pub const zhseqr = LAPACKE_zhseqr;

extern fn LAPACKE_clacgv(n: isize, x: [*c]cf32, incx: isize) isize;
extern fn LAPACKE_zlacgv(n: isize, x: [*c]cf64, incx: isize) isize;
pub const clacgv = LAPACKE_clacgv;
pub const zlacgv = LAPACKE_zlacgv;

extern fn LAPACKE_slacn2(n: isize, v: [*c]f32, x: [*c]f32, isgn: [*c]isize, est: [*c]f32, kase: [*c]isize, isave: [*c]isize) isize;
extern fn LAPACKE_dlacn2(n: isize, v: [*c]f64, x: [*c]f64, isgn: [*c]isize, est: [*c]f64, kase: [*c]isize, isave: [*c]isize) isize;
extern fn LAPACKE_clacn2(n: isize, v: [*c]cf32, x: [*c]cf32, est: [*c]f32, kase: [*c]isize, isave: [*c]isize) isize;
extern fn LAPACKE_zlacn2(n: isize, v: [*c]cf64, x: [*c]cf64, est: [*c]f64, kase: [*c]isize, isave: [*c]isize) isize;
pub const slacn2 = LAPACKE_slacn2;
pub const dlacn2 = LAPACKE_dlacn2;
pub const clacn2 = LAPACKE_clacn2;
pub const zlacn2 = LAPACKE_zlacn2;

extern fn LAPACKE_slacpy(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dlacpy(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_clacpy(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zlacpy(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const slacpy = LAPACKE_slacpy;
pub const dlacpy = LAPACKE_dlacpy;
pub const clacpy = LAPACKE_clacpy;
pub const zlacpy = LAPACKE_zlacpy;

extern fn LAPACKE_clacp2(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zlacp2(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const clacp2 = LAPACKE_clacp2;
pub const zlacp2 = LAPACKE_zlacp2;

extern fn LAPACKE_slag2d(layout: c_int, m: isize, n: isize, sa: [*c]const f32, ldsa: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_dlag2s(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, sa: [*c]f32, ldsa: isize) isize;
extern fn LAPACKE_clag2z(layout: c_int, m: isize, n: isize, sa: [*c]const cf32, ldsa: isize, a: [*c]cf64, lda: isize) isize;
extern fn LAPACKE_zlag2c(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, sa: [*c]cf32, ldsa: isize) isize;
pub const slag2d = LAPACKE_slag2d;
pub const dlag2s = LAPACKE_dlag2s;
pub const clag2z = LAPACKE_clag2z;
pub const zlag2c = LAPACKE_zlag2c;

extern fn LAPACKE_slagge(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f32, a: [*c]f32, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_dlagge(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f64, a: [*c]f64, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_clagge(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f32, a: [*c]cf32, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_zlagge(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f64, a: [*c]cf64, lda: isize, iseed: [*c]isize) isize;
pub const slagge = LAPACKE_slagge;
pub const dlagge = LAPACKE_dlagge;
pub const clagge = LAPACKE_clagge;
pub const zlagge = LAPACKE_zlagge;

extern fn LAPACKE_slamch(cmach: u8) f32;
extern fn LAPACKE_dlamch(cmach: u8) f64;
pub const slamch = LAPACKE_slamch;
pub const dlamch = LAPACKE_dlamch;

extern fn LAPACKE_slangb(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize) f32;
extern fn LAPACKE_dlangb(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize) f64;
extern fn LAPACKE_clangb(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize) f32;
extern fn LAPACKE_zlangb(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize) f64;
pub const slangb = LAPACKE_slangb;
pub const dlangb = LAPACKE_dlangb;
pub const clangb = LAPACKE_clangb;
pub const zlangb = LAPACKE_zlangb;

extern fn LAPACKE_slange(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const f32, lda: isize) f32;
extern fn LAPACKE_dlange(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const f64, lda: isize) f64;
extern fn LAPACKE_clange(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize) f32;
extern fn LAPACKE_zlange(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize) f64;
pub const slange = LAPACKE_slange;
pub const dlange = LAPACKE_dlange;
pub const clange = LAPACKE_clange;
pub const zlange = LAPACKE_zlange;

extern fn LAPACKE_clanhe(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf32, lda: isize) f32;
extern fn LAPACKE_zlanhe(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf64, lda: isize) f64;
pub const clanhe = LAPACKE_clanhe;
pub const zlanhe = LAPACKE_zlanhe;

extern fn LAPACKE_clacrm(layout: c_int, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zlacrm(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]cf64, ldc: isize) isize;
pub const clacrm = LAPACKE_clacrm;
pub const zlacrm = LAPACKE_zlacrm;

extern fn LAPACKE_clarcm(layout: c_int, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zlarcm(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize) isize;
pub const clarcm = LAPACKE_clarcm;
pub const zlarcm = LAPACKE_zlarcm;

extern fn LAPACKE_slansy(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const f32, lda: isize) f32;
extern fn LAPACKE_dlansy(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const f64, lda: isize) f64;
extern fn LAPACKE_clansy(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf32, lda: isize) f32;
extern fn LAPACKE_zlansy(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf64, lda: isize) f64;
pub const slansy = LAPACKE_slansy;
pub const dlansy = LAPACKE_dlansy;
pub const clansy = LAPACKE_clansy;
pub const zlansy = LAPACKE_zlansy;

extern fn LAPACKE_slantr(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const f32, lda: isize) f32;
extern fn LAPACKE_dlantr(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const f64, lda: isize) f64;
extern fn LAPACKE_clantr(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize) f32;
extern fn LAPACKE_zlantr(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize) f64;
pub const slantr = LAPACKE_slantr;
pub const dlantr = LAPACKE_dlantr;
pub const clantr = LAPACKE_clantr;
pub const zlantr = LAPACKE_zlantr;

extern fn LAPACKE_slarfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dlarfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, c: [*c]f64, ldc: isize) isize;
extern fn LAPACKE_clarfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zlarfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, c: [*c]cf64, ldc: isize) isize;
pub const slarfb = LAPACKE_slarfb;
pub const dlarfb = LAPACKE_dlarfb;
pub const clarfb = LAPACKE_clarfb;
pub const zlarfb = LAPACKE_zlarfb;

extern fn LAPACKE_slarfg(n: isize, alpha: [*c]f32, x: [*c]f32, incx: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dlarfg(n: isize, alpha: [*c]f64, x: [*c]f64, incx: isize, tau: [*c]f64) isize;
extern fn LAPACKE_clarfg(n: isize, alpha: [*c]cf32, x: [*c]cf32, incx: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zlarfg(n: isize, alpha: [*c]cf64, x: [*c]cf64, incx: isize, tau: [*c]cf64) isize;
pub const slarfg = LAPACKE_slarfg;
pub const dlarfg = LAPACKE_dlarfg;
pub const clarfg = LAPACKE_clarfg;
pub const zlarfg = LAPACKE_zlarfg;

extern fn LAPACKE_slarft(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const f32, ldv: isize, tau: [*c]const f32, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dlarft(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const f64, ldv: isize, tau: [*c]const f64, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_clarft(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const cf32, ldv: isize, tau: [*c]const cf32, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zlarft(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const cf64, ldv: isize, tau: [*c]const cf64, t: [*c]cf64, ldt: isize) isize;
pub const slarft = LAPACKE_slarft;
pub const dlarft = LAPACKE_dlarft;
pub const clarft = LAPACKE_clarft;
pub const zlarft = LAPACKE_zlarft;

extern fn LAPACKE_slarfx(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const f32, tau: f32, c: [*c]f32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_dlarfx(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const f64, tau: f64, c: [*c]f64, ldc: isize, work: [*c]f64) isize;
extern fn LAPACKE_clarfx(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const cf32, tau: cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zlarfx(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const cf64, tau: cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64) isize;
pub const slarfx = LAPACKE_slarfx;
pub const dlarfx = LAPACKE_dlarfx;
pub const clarfx = LAPACKE_clarfx;
pub const zlarfx = LAPACKE_zlarfx;

extern fn LAPACKE_slarnv(idist: isize, iseed: [*c]isize, n: isize, x: [*c]f32) isize;
extern fn LAPACKE_dlarnv(idist: isize, iseed: [*c]isize, n: isize, x: [*c]f64) isize;
extern fn LAPACKE_clarnv(idist: isize, iseed: [*c]isize, n: isize, x: [*c]cf32) isize;
extern fn LAPACKE_zlarnv(idist: isize, iseed: [*c]isize, n: isize, x: [*c]cf64) isize;
pub const slarnv = LAPACKE_slarnv;
pub const dlarnv = LAPACKE_dlarnv;
pub const clarnv = LAPACKE_clarnv;
pub const zlarnv = LAPACKE_zlarnv;

extern fn LAPACKE_slascl(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f32, cto: f32, m: isize, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlascl(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f64, cto: f64, m: isize, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_clascl(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f32, cto: f32, m: isize, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlascl(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f64, cto: f64, m: isize, n: isize, a: [*c]cf64, lda: isize) isize;
pub const slascl = LAPACKE_slascl;
pub const dlascl = LAPACKE_dlascl;
pub const clascl = LAPACKE_clascl;
pub const zlascl = LAPACKE_zlascl;

extern fn LAPACKE_slaset(layout: c_int, uplo: u8, m: isize, n: isize, alpha: f32, beta: f32, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlaset(layout: c_int, uplo: u8, m: isize, n: isize, alpha: f64, beta: f64, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_claset(layout: c_int, uplo: u8, m: isize, n: isize, alpha: cf32, beta: cf32, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlaset(layout: c_int, uplo: u8, m: isize, n: isize, alpha: cf64, beta: cf64, a: [*c]cf64, lda: isize) isize;
pub const slaset = LAPACKE_slaset;
pub const dlaset = LAPACKE_dlaset;
pub const claset = LAPACKE_claset;
pub const zlaset = LAPACKE_zlaset;

extern fn LAPACKE_slasrt(id: u8, n: isize, d: [*c]f32) isize;
extern fn LAPACKE_dlasrt(id: u8, n: isize, d: [*c]f64) isize;
pub const slasrt = LAPACKE_slasrt;
pub const dlasrt = LAPACKE_dlasrt;

extern fn LAPACKE_slassq(n: isize, x: [*c]f32, incx: isize, scale: [*c]f32, sumsq: [*c]f32) isize;
extern fn LAPACKE_dlassq(n: isize, x: [*c]f64, incx: isize, scale: [*c]f64, sumsq: [*c]f64) isize;
extern fn LAPACKE_classq(n: isize, x: [*c]cf32, incx: isize, scale: [*c]f32, sumsq: [*c]f32) isize;
extern fn LAPACKE_zlassq(n: isize, x: [*c]cf64, incx: isize, scale: [*c]f64, sumsq: [*c]f64) isize;
pub const slassq = LAPACKE_slassq;
pub const dlassq = LAPACKE_dlassq;
pub const classq = LAPACKE_classq;
pub const zlassq = LAPACKE_zlassq;

extern fn LAPACKE_slaswp(layout: c_int, n: isize, a: [*c]f32, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
extern fn LAPACKE_dlaswp(layout: c_int, n: isize, a: [*c]f64, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
extern fn LAPACKE_claswp(layout: c_int, n: isize, a: [*c]cf32, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
extern fn LAPACKE_zlaswp(layout: c_int, n: isize, a: [*c]cf64, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
pub const slaswp = LAPACKE_slaswp;
pub const dlaswp = LAPACKE_dlaswp;
pub const claswp = LAPACKE_claswp;
pub const zlaswp = LAPACKE_zlaswp;

extern fn LAPACKE_slatms(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f32, mode: isize, cond: f32, dmax: f32, kl: isize, ku: isize, pack: u8, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlatms(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f64, mode: isize, cond: f64, dmax: f64, kl: isize, ku: isize, pack: u8, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_clatms(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f32, mode: isize, cond: f32, dmax: f32, kl: isize, ku: isize, pack: u8, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlatms(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f64, mode: isize, cond: f64, dmax: f64, kl: isize, ku: isize, pack: u8, a: [*c]cf64, lda: isize) isize;
pub const slatms = LAPACKE_slatms;
pub const dlatms = LAPACKE_dlatms;
pub const clatms = LAPACKE_clatms;
pub const zlatms = LAPACKE_zlatms;

extern fn LAPACKE_slauum(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlauum(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_clauum(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlauum(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const slauum = LAPACKE_slauum;
pub const dlauum = LAPACKE_dlauum;
pub const clauum = LAPACKE_clauum;
pub const zlauum = LAPACKE_zlauum;

extern fn LAPACKE_sopgtr(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, tau: [*c]const f32, q: [*c]f32, ldq: isize) isize;
extern fn LAPACKE_dopgtr(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, tau: [*c]const f64, q: [*c]f64, ldq: isize) isize;
pub const sopgtr = LAPACKE_sopgtr;
pub const dopgtr = LAPACKE_dopgtr;

extern fn LAPACKE_sopmtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const f32, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dopmtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const f64, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sopmtr = LAPACKE_sopmtr;
pub const dopmtr = LAPACKE_dopmtr;

extern fn LAPACKE_sorgbr(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorgbr(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorgbr = LAPACKE_sorgbr;
pub const dorgbr = LAPACKE_dorgbr;

extern fn LAPACKE_sorghr(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorghr(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorghr = LAPACKE_sorghr;
pub const dorghr = LAPACKE_dorghr;

extern fn LAPACKE_sorglq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorglq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorglq = LAPACKE_sorglq;
pub const dorglq = LAPACKE_dorglq;

extern fn LAPACKE_sorgql(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorgql(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorgql = LAPACKE_sorgql;
pub const dorgql = LAPACKE_dorgql;

extern fn LAPACKE_sorgqr(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorgqr(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorgqr = LAPACKE_sorgqr;
pub const dorgqr = LAPACKE_dorgqr;

extern fn LAPACKE_sorgrq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorgrq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorgrq = LAPACKE_sorgrq;
pub const dorgrq = LAPACKE_dorgrq;

extern fn LAPACKE_sorgtr(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, tau: [*c]const f32) isize;
extern fn LAPACKE_dorgtr(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, tau: [*c]const f64) isize;
pub const sorgtr = LAPACKE_sorgtr;
pub const dorgtr = LAPACKE_dorgtr;

extern fn LAPACKE_sorgtsqr_row(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]f32, lda: isize, t: [*c]const f32, ldt: isize) isize;
extern fn LAPACKE_dorgtsqr_row(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]f64, lda: isize, t: [*c]const f64, ldt: isize) isize;
pub const sorgtsqr_row = LAPACKE_sorgtsqr_row;
pub const dorgtsqr_row = LAPACKE_dorgtsqr_row;

extern fn LAPACKE_sormbr(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormbr(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormbr = LAPACKE_sormbr;
pub const dormbr = LAPACKE_dormbr;

extern fn LAPACKE_sormhr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormhr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormhr = LAPACKE_sormhr;
pub const dormhr = LAPACKE_dormhr;

extern fn LAPACKE_sormlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormlq = LAPACKE_sormlq;
pub const dormlq = LAPACKE_dormlq;

extern fn LAPACKE_sormql(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormql(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormql = LAPACKE_sormql;
pub const dormql = LAPACKE_dormql;

extern fn LAPACKE_sormqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormqr = LAPACKE_sormqr;
pub const dormqr = LAPACKE_dormqr;

extern fn LAPACKE_sormrq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormrq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormrq = LAPACKE_sormrq;
pub const dormrq = LAPACKE_dormrq;

extern fn LAPACKE_sormrz(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormrz(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormrz = LAPACKE_sormrz;
pub const dormrz = LAPACKE_dormrz;

extern fn LAPACKE_sormtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dormtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize) isize;
pub const sormtr = LAPACKE_sormtr;
pub const dormtr = LAPACKE_dormtr;

extern fn LAPACKE_spbcon(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f32, ldab: isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dpbcon(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f64, ldab: isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cpbcon(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf32, ldab: isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zpbcon(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf64, ldab: isize, anorm: f64, rcond: [*c]f64) isize;
pub const spbcon = LAPACKE_spbcon;
pub const dpbcon = LAPACKE_dpbcon;
pub const cpbcon = LAPACKE_cpbcon;
pub const zpbcon = LAPACKE_zpbcon;

extern fn LAPACKE_spbequ(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f32, ldab: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dpbequ(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f64, ldab: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cpbequ(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf32, ldab: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zpbequ(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf64, ldab: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const spbequ = LAPACKE_spbequ;
pub const dpbequ = LAPACKE_dpbequ;
pub const cpbequ = LAPACKE_cpbequ;
pub const zpbequ = LAPACKE_zpbequ;

extern fn LAPACKE_spbrfs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, afb: [*c]const f32, ldafb: isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dpbrfs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, afb: [*c]const f64, ldafb: isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cpbrfs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, afb: [*c]const cf32, ldafb: isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zpbrfs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, afb: [*c]const cf64, ldafb: isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const spbrfs = LAPACKE_spbrfs;
pub const dpbrfs = LAPACKE_dpbrfs;
pub const cpbrfs = LAPACKE_cpbrfs;
pub const zpbrfs = LAPACKE_zpbrfs;

extern fn LAPACKE_spbstf(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]f32, ldbb: isize) isize;
extern fn LAPACKE_dpbstf(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]f64, ldbb: isize) isize;
extern fn LAPACKE_cpbstf(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]cf32, ldbb: isize) isize;
extern fn LAPACKE_zpbstf(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]cf64, ldbb: isize) isize;
pub const spbstf = LAPACKE_spbstf;
pub const dpbstf = LAPACKE_dpbstf;
pub const cpbstf = LAPACKE_cpbstf;
pub const zpbstf = LAPACKE_zpbstf;

extern fn LAPACKE_spbsv(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f32, ldab: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpbsv(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f64, ldab: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpbsv(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpbsv(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, b: [*c]cf64, ldb: isize) isize;
pub const spbsv = LAPACKE_spbsv;
pub const dpbsv = LAPACKE_dpbsv;
pub const cpbsv = LAPACKE_cpbsv;
pub const zpbsv = LAPACKE_zpbsv;

extern fn LAPACKE_spbsvx(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f32, ldab: isize, afb: [*c]f32, ldafb: isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dpbsvx(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f64, ldab: isize, afb: [*c]f64, ldafb: isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cpbsvx(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, afb: [*c]cf32, ldafb: isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zpbsvx(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, afb: [*c]cf64, ldafb: isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const spbsvx = LAPACKE_spbsvx;
pub const dpbsvx = LAPACKE_dpbsvx;
pub const cpbsvx = LAPACKE_cpbsvx;
pub const zpbsvx = LAPACKE_zpbsvx;

extern fn LAPACKE_spbtrf(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize) isize;
extern fn LAPACKE_dpbtrf(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize) isize;
extern fn LAPACKE_cpbtrf(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize) isize;
extern fn LAPACKE_zpbtrf(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize) isize;
pub const spbtrf = LAPACKE_spbtrf;
pub const dpbtrf = LAPACKE_dpbtrf;
pub const cpbtrf = LAPACKE_cpbtrf;
pub const zpbtrf = LAPACKE_zpbtrf;

extern fn LAPACKE_spbtrs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpbtrs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpbtrs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpbtrs(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, b: [*c]cf64, ldb: isize) isize;
pub const spbtrs = LAPACKE_spbtrs;
pub const dpbtrs = LAPACKE_dpbtrs;
pub const cpbtrs = LAPACKE_cpbtrs;
pub const zpbtrs = LAPACKE_zpbtrs;

extern fn LAPACKE_spftrf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f32) isize;
extern fn LAPACKE_dpftrf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f64) isize;
extern fn LAPACKE_cpftrf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf32) isize;
extern fn LAPACKE_zpftrf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf64) isize;
pub const spftrf = LAPACKE_spftrf;
pub const dpftrf = LAPACKE_dpftrf;
pub const cpftrf = LAPACKE_cpftrf;
pub const zpftrf = LAPACKE_zpftrf;

extern fn LAPACKE_spftri(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f32) isize;
extern fn LAPACKE_dpftri(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f64) isize;
extern fn LAPACKE_cpftri(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf32) isize;
extern fn LAPACKE_zpftri(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf64) isize;
pub const spftri = LAPACKE_spftri;
pub const dpftri = LAPACKE_dpftri;
pub const cpftri = LAPACKE_cpftri;
pub const zpftri = LAPACKE_zpftri;

extern fn LAPACKE_spftrs(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpftrs(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpftrs(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpftrs(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const spftrs = LAPACKE_spftrs;
pub const dpftrs = LAPACKE_dpftrs;
pub const cpftrs = LAPACKE_cpftrs;
pub const zpftrs = LAPACKE_zpftrs;

extern fn LAPACKE_spocon(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dpocon(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cpocon(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zpocon(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, anorm: f64, rcond: [*c]f64) isize;
pub const spocon = LAPACKE_spocon;
pub const dpocon = LAPACKE_dpocon;
pub const cpocon = LAPACKE_cpocon;
pub const zpocon = LAPACKE_zpocon;

extern fn LAPACKE_spoequ(layout: c_int, n: isize, a: [*c]const f32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dpoequ(layout: c_int, n: isize, a: [*c]const f64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cpoequ(layout: c_int, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zpoequ(layout: c_int, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const spoequ = LAPACKE_spoequ;
pub const dpoequ = LAPACKE_dpoequ;
pub const cpoequ = LAPACKE_cpoequ;
pub const zpoequ = LAPACKE_zpoequ;

extern fn LAPACKE_spoequb(layout: c_int, n: isize, a: [*c]const f32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dpoequb(layout: c_int, n: isize, a: [*c]const f64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cpoequb(layout: c_int, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zpoequb(layout: c_int, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const spoequb = LAPACKE_spoequb;
pub const dpoequb = LAPACKE_dpoequb;
pub const cpoequb = LAPACKE_cpoequb;
pub const zpoequb = LAPACKE_zpoequb;

extern fn LAPACKE_sporfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dporfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cporfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zporfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sporfs = LAPACKE_sporfs;
pub const dporfs = LAPACKE_dporfs;
pub const cporfs = LAPACKE_cporfs;
pub const zporfs = LAPACKE_zporfs;

extern fn LAPACKE_sporfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, s: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dporfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, s: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_cporfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, s: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zporfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, s: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const sporfsx = LAPACKE_sporfsx;
pub const dporfsx = LAPACKE_dporfsx;
pub const cporfsx = LAPACKE_cporfsx;
pub const zporfsx = LAPACKE_zporfsx;

extern fn LAPACKE_sposv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dposv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cposv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zposv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
extern fn LAPACKE_dsposv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, iter: [*c]isize) isize;
extern fn LAPACKE_zcposv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, iter: [*c]isize) isize;
pub const sposv = LAPACKE_sposv;
pub const dposv = LAPACKE_dposv;
pub const cposv = LAPACKE_cposv;
pub const zposv = LAPACKE_zposv;
pub const dsposv = LAPACKE_dsposv;
pub const zcposv = LAPACKE_zcposv;

extern fn LAPACKE_sposvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dposvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cposvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zposvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sposvx = LAPACKE_sposvx;
pub const dposvx = LAPACKE_dposvx;
pub const cposvx = LAPACKE_cposvx;
pub const zposvx = LAPACKE_zposvx;

extern fn LAPACKE_sposvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dposvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_cposvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zposvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const sposvxx = LAPACKE_sposvxx;
pub const dposvxx = LAPACKE_dposvxx;
pub const cposvxx = LAPACKE_cposvxx;
pub const zposvxx = LAPACKE_zposvxx;

extern fn LAPACKE_spotrf2(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dpotrf2(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_cpotrf2(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zpotrf2(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const spotrf2 = LAPACKE_spotrf2;
pub const dpotrf2 = LAPACKE_dpotrf2;
pub const cpotrf2 = LAPACKE_cpotrf2;
pub const zpotrf2 = LAPACKE_zpotrf2;

extern fn LAPACKE_spotrf(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dpotrf(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_cpotrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zpotrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const spotrf = LAPACKE_spotrf;
pub const dpotrf = LAPACKE_dpotrf;
pub const cpotrf = LAPACKE_cpotrf;
pub const zpotrf = LAPACKE_zpotrf;

extern fn LAPACKE_spotri(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dpotri(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_cpotri(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zpotri(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const spotri = LAPACKE_spotri;
pub const dpotri = LAPACKE_dpotri;
pub const cpotri = LAPACKE_cpotri;
pub const zpotri = LAPACKE_zpotri;

extern fn LAPACKE_spotrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpotrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpotrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpotrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const spotrs = LAPACKE_spotrs;
pub const dpotrs = LAPACKE_dpotrs;
pub const cpotrs = LAPACKE_cpotrs;
pub const zpotrs = LAPACKE_zpotrs;

extern fn LAPACKE_sppcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dppcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cppcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zppcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, anorm: f64, rcond: [*c]f64) isize;
pub const sppcon = LAPACKE_sppcon;
pub const dppcon = LAPACKE_dppcon;
pub const cppcon = LAPACKE_cppcon;
pub const zppcon = LAPACKE_zppcon;

extern fn LAPACKE_sppequ(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dppequ(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cppequ(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zppequ(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const sppequ = LAPACKE_sppequ;
pub const dppequ = LAPACKE_dppequ;
pub const cppequ = LAPACKE_cppequ;
pub const zppequ = LAPACKE_zppequ;

extern fn LAPACKE_spprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, afp: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dpprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, afp: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cpprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]const cf32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zpprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]const cf64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const spprfs = LAPACKE_spprfs;
pub const dpprfs = LAPACKE_dpprfs;
pub const cpprfs = LAPACKE_cpprfs;
pub const zpprfs = LAPACKE_zpprfs;

extern fn LAPACKE_sppsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dppsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cppsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zppsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, b: [*c]cf64, ldb: isize) isize;
pub const sppsv = LAPACKE_sppsv;
pub const dppsv = LAPACKE_dppsv;
pub const cppsv = LAPACKE_cppsv;
pub const zppsv = LAPACKE_zppsv;

extern fn LAPACKE_sppsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]f32, afp: [*c]f32, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dppsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]f64, afp: [*c]f64, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cppsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, afp: [*c]cf32, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zppsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, afp: [*c]cf64, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sppsvx = LAPACKE_sppsvx;
pub const dppsvx = LAPACKE_dppsvx;
pub const cppsvx = LAPACKE_cppsvx;
pub const zppsvx = LAPACKE_zppsvx;

extern fn LAPACKE_spptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dpptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]f64) isize;
extern fn LAPACKE_cpptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_zpptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64) isize;
pub const spptrf = LAPACKE_spptrf;
pub const dpptrf = LAPACKE_dpptrf;
pub const cpptrf = LAPACKE_cpptrf;
pub const zpptrf = LAPACKE_zpptrf;

extern fn LAPACKE_spptri(layout: c_int, uplo: u8, n: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dpptri(layout: c_int, uplo: u8, n: isize, ap: [*c]f64) isize;
extern fn LAPACKE_cpptri(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_zpptri(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64) isize;
pub const spptri = LAPACKE_spptri;
pub const dpptri = LAPACKE_dpptri;
pub const cpptri = LAPACKE_cpptri;
pub const zpptri = LAPACKE_zpptri;

extern fn LAPACKE_spptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const spptrs = LAPACKE_spptrs;
pub const dpptrs = LAPACKE_dpptrs;
pub const cpptrs = LAPACKE_cpptrs;
pub const zpptrs = LAPACKE_zpptrs;

extern fn LAPACKE_spstrf(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f32) isize;
extern fn LAPACKE_dpstrf(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f64) isize;
extern fn LAPACKE_cpstrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f32) isize;
extern fn LAPACKE_zpstrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f64) isize;
pub const spstrf = LAPACKE_spstrf;
pub const dpstrf = LAPACKE_dpstrf;
pub const cpstrf = LAPACKE_cpstrf;
pub const zpstrf = LAPACKE_zpstrf;

extern fn LAPACKE_sptcon(n: isize, d: [*c]const f32, e: [*c]const f32, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dptcon(n: isize, d: [*c]const f64, e: [*c]const f64, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cptcon(n: isize, d: [*c]const f32, e: [*c]const cf32, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zptcon(n: isize, d: [*c]const f64, e: [*c]const cf64, anorm: f64, rcond: [*c]f64) isize;
pub const sptcon = LAPACKE_sptcon;
pub const dptcon = LAPACKE_dptcon;
pub const cptcon = LAPACKE_cptcon;
pub const zptcon = LAPACKE_zptcon;

extern fn LAPACKE_spteqr(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dpteqr(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_cpteqr(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zpteqr(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const spteqr = LAPACKE_spteqr;
pub const dpteqr = LAPACKE_dpteqr;
pub const cpteqr = LAPACKE_cpteqr;
pub const zpteqr = LAPACKE_zpteqr;

extern fn LAPACKE_sptrfs(layout: c_int, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const f32, df: [*c]const f32, ef: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dptrfs(layout: c_int, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const f64, df: [*c]const f64, ef: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cptrfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const cf32, df: [*c]const f32, ef: [*c]const cf32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zptrfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const cf64, df: [*c]const f64, ef: [*c]const cf64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sptrfs = LAPACKE_sptrfs;
pub const dptrfs = LAPACKE_dptrfs;
pub const cptrfs = LAPACKE_cptrfs;
pub const zptrfs = LAPACKE_zptrfs;

extern fn LAPACKE_sptsv(layout: c_int, n: isize, nrhs: isize, d: [*c]f32, e: [*c]f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dptsv(layout: c_int, n: isize, nrhs: isize, d: [*c]f64, e: [*c]f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cptsv(layout: c_int, n: isize, nrhs: isize, d: [*c]f32, e: [*c]cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zptsv(layout: c_int, n: isize, nrhs: isize, d: [*c]f64, e: [*c]cf64, b: [*c]cf64, ldb: isize) isize;
pub const sptsv = LAPACKE_sptsv;
pub const dptsv = LAPACKE_dptsv;
pub const cptsv = LAPACKE_cptsv;
pub const zptsv = LAPACKE_zptsv;

extern fn LAPACKE_sptsvx(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const f32, df: [*c]f32, ef: [*c]f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dptsvx(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const f64, df: [*c]f64, ef: [*c]f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cptsvx(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const cf32, df: [*c]f32, ef: [*c]cf32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zptsvx(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const cf64, df: [*c]f64, ef: [*c]cf64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sptsvx = LAPACKE_sptsvx;
pub const dptsvx = LAPACKE_dptsvx;
pub const cptsvx = LAPACKE_cptsvx;
pub const zptsvx = LAPACKE_zptsvx;

extern fn LAPACKE_spttrf(n: isize, d: [*c]f32, e: [*c]f32) isize;
extern fn LAPACKE_dpttrf(n: isize, d: [*c]f64, e: [*c]f64) isize;
extern fn LAPACKE_cpttrf(n: isize, d: [*c]f32, e: [*c]cf32) isize;
extern fn LAPACKE_zpttrf(n: isize, d: [*c]f64, e: [*c]cf64) isize;
pub const spttrf = LAPACKE_spttrf;
pub const dpttrf = LAPACKE_dpttrf;
pub const cpttrf = LAPACKE_cpttrf;
pub const zpttrf = LAPACKE_zpttrf;

extern fn LAPACKE_spttrs(layout: c_int, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpttrs(layout: c_int, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpttrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpttrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const spttrs = LAPACKE_spttrs;
pub const dpttrs = LAPACKE_dpttrs;
pub const cpttrs = LAPACKE_cpttrs;
pub const zpttrs = LAPACKE_zpttrs;

extern fn LAPACKE_ssbev(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsbev(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const ssbev = LAPACKE_ssbev;
pub const dsbev = LAPACKE_dsbev;

extern fn LAPACKE_ssbevd(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsbevd(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const ssbevd = LAPACKE_ssbevd;
pub const dsbevd = LAPACKE_dsbevd;

extern fn LAPACKE_ssbevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, q: [*c]f32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsbevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, q: [*c]f64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const ssbevx = LAPACKE_ssbevx;
pub const dsbevx = LAPACKE_dsbevx;

extern fn LAPACKE_ssbgst(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]const f32, ldbb: isize, x: [*c]f32, ldx: isize) isize;
extern fn LAPACKE_dsbgst(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]const f64, ldbb: isize, x: [*c]f64, ldx: isize) isize;
pub const ssbgst = LAPACKE_ssbgst;
pub const dsbgst = LAPACKE_dsbgst;

extern fn LAPACKE_ssbgv(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]f32, ldbb: isize, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsbgv(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]f64, ldbb: isize, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const ssbgv = LAPACKE_ssbgv;
pub const dsbgv = LAPACKE_dsbgv;

extern fn LAPACKE_ssbgvd(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]f32, ldbb: isize, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsbgvd(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]f64, ldbb: isize, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const ssbgvd = LAPACKE_ssbgvd;
pub const dsbgvd = LAPACKE_dsbgvd;

extern fn LAPACKE_ssbgvx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]f32, ldbb: isize, q: [*c]f32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsbgvx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]f64, ldbb: isize, q: [*c]f64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const ssbgvx = LAPACKE_ssbgvx;
pub const dsbgvx = LAPACKE_dsbgvx;

extern fn LAPACKE_ssbtrd(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]f32, ldq: isize) isize;
extern fn LAPACKE_dsbtrd(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]f64, ldq: isize) isize;
pub const ssbtrd = LAPACKE_ssbtrd;
pub const dsbtrd = LAPACKE_dsbtrd;

extern fn LAPACKE_ssfrk(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, beta: f32, c: [*c]f32) isize;
extern fn LAPACKE_dsfrk(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, beta: f64, c: [*c]f64) isize;
pub const ssfrk = LAPACKE_ssfrk;
pub const dsfrk = LAPACKE_dsfrk;

extern fn LAPACKE_sspcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dspcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_cspcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zspcon(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const sspcon = LAPACKE_sspcon;
pub const dspcon = LAPACKE_dspcon;
pub const cspcon = LAPACKE_cspcon;
pub const zspcon = LAPACKE_zspcon;

extern fn LAPACKE_sspev(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dspev(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const sspev = LAPACKE_sspev;
pub const dspev = LAPACKE_dspev;

extern fn LAPACKE_sspevd(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dspevd(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const sspevd = LAPACKE_sspevd;
pub const dspevd = LAPACKE_dspevd;

extern fn LAPACKE_sspevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dspevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const sspevx = LAPACKE_sspevx;
pub const dspevx = LAPACKE_dspevx;

extern fn LAPACKE_sspgst(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]const f32) isize;
extern fn LAPACKE_dspgst(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]const f64) isize;
pub const sspgst = LAPACKE_sspgst;
pub const dspgst = LAPACKE_dspgst;

extern fn LAPACKE_sspgv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dspgv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const sspgv = LAPACKE_sspgv;
pub const dspgv = LAPACKE_dspgv;

extern fn LAPACKE_sspgvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dspgvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const sspgvd = LAPACKE_sspgvd;
pub const dspgvd = LAPACKE_dspgvd;

extern fn LAPACKE_sspgvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dspgvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const sspgvx = LAPACKE_sspgvx;
pub const dspgvx = LAPACKE_dspgvx;

extern fn LAPACKE_ssprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, afp: [*c]const f32, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dsprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, afp: [*c]const f64, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_csprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]const cf32, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zsprfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]const cf64, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const ssprfs = LAPACKE_ssprfs;
pub const dsprfs = LAPACKE_dsprfs;
pub const csprfs = LAPACKE_csprfs;
pub const zsprfs = LAPACKE_zsprfs;

extern fn LAPACKE_sspsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f32, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dspsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f64, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cspsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zspsv(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const sspsv = LAPACKE_sspsv;
pub const dspsv = LAPACKE_dspsv;
pub const cspsv = LAPACKE_cspsv;
pub const zspsv = LAPACKE_zspsv;

extern fn LAPACKE_sspsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, afp: [*c]f32, ipiv: [*c]isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dspsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, afp: [*c]f64, ipiv: [*c]isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_cspsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]cf32, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zspsvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]cf64, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const sspsvx = LAPACKE_sspsvx;
pub const dspsvx = LAPACKE_dspsvx;
pub const cspsvx = LAPACKE_cspsvx;
pub const zspsvx = LAPACKE_zspsvx;

extern fn LAPACKE_ssptrd(layout: c_int, uplo: u8, n: isize, ap: [*c]f32, d: [*c]f32, e: [*c]f32, tau: [*c]f32) isize;
extern fn LAPACKE_dsptrd(layout: c_int, uplo: u8, n: isize, ap: [*c]f64, d: [*c]f64, e: [*c]f64, tau: [*c]f64) isize;
pub const ssptrd = LAPACKE_ssptrd;
pub const dsptrd = LAPACKE_dsptrd;

extern fn LAPACKE_ssptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]f32, ipiv: [*c]isize) isize;
extern fn LAPACKE_dsptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]f64, ipiv: [*c]isize) isize;
extern fn LAPACKE_csptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zsptrf(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]isize) isize;
pub const ssptrf = LAPACKE_ssptrf;
pub const dsptrf = LAPACKE_dsptrf;
pub const csptrf = LAPACKE_csptrf;
pub const zsptrf = LAPACKE_zsptrf;

extern fn LAPACKE_ssptri(layout: c_int, uplo: u8, n: isize, ap: [*c]f32, ipiv: [*c]const isize) isize;
extern fn LAPACKE_dsptri(layout: c_int, uplo: u8, n: isize, ap: [*c]f64, ipiv: [*c]const isize) isize;
extern fn LAPACKE_csptri(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zsptri(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]const isize) isize;
pub const ssptri = LAPACKE_ssptri;
pub const dsptri = LAPACKE_dsptri;
pub const csptri = LAPACKE_csptri;
pub const zsptri = LAPACKE_zsptri;

extern fn LAPACKE_ssptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsptrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const ssptrs = LAPACKE_ssptrs;
pub const dsptrs = LAPACKE_dsptrs;
pub const csptrs = LAPACKE_csptrs;
pub const zsptrs = LAPACKE_zsptrs;

extern fn LAPACKE_sstebz(range: u8, order: u8, n: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, d: [*c]const f32, e: [*c]const f32, m: [*c]isize, nsplit: [*c]isize, w: [*c]f32, iblock: [*c]isize, isplit: [*c]isize) isize;
extern fn LAPACKE_dstebz(range: u8, order: u8, n: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, d: [*c]const f64, e: [*c]const f64, m: [*c]isize, nsplit: [*c]isize, w: [*c]f64, iblock: [*c]isize, isplit: [*c]isize) isize;
pub const sstebz = LAPACKE_sstebz;
pub const dstebz = LAPACKE_dstebz;

extern fn LAPACKE_sstedc(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dstedc(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_cstedc(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zstedc(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const sstedc = LAPACKE_sstedc;
pub const dstedc = LAPACKE_dstedc;
pub const cstedc = LAPACKE_cstedc;
pub const zstedc = LAPACKE_zstedc;

extern fn LAPACKE_sstegr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_dstegr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_cstegr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_zstegr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, isuppz: [*c]isize) isize;
pub const sstegr = LAPACKE_sstegr;
pub const dstegr = LAPACKE_dstegr;
pub const cstegr = LAPACKE_cstegr;
pub const zstegr = LAPACKE_zstegr;

extern fn LAPACKE_sstein(layout: c_int, n: isize, d: [*c]const f32, e: [*c]const f32, m: isize, w: [*c]const f32, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]f32, ldz: isize, ifailv: [*c]isize) isize;
extern fn LAPACKE_dstein(layout: c_int, n: isize, d: [*c]const f64, e: [*c]const f64, m: isize, w: [*c]const f64, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]f64, ldz: isize, ifailv: [*c]isize) isize;
extern fn LAPACKE_cstein(layout: c_int, n: isize, d: [*c]const f32, e: [*c]const f32, m: isize, w: [*c]const f32, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]cf32, ldz: isize, ifailv: [*c]isize) isize;
extern fn LAPACKE_zstein(layout: c_int, n: isize, d: [*c]const f64, e: [*c]const f64, m: isize, w: [*c]const f64, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]cf64, ldz: isize, ifailv: [*c]isize) isize;
pub const sstein = LAPACKE_sstein;
pub const dstein = LAPACKE_dstein;
pub const cstein = LAPACKE_cstein;
pub const zstein = LAPACKE_zstein;

extern fn LAPACKE_sstemr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize) isize;
extern fn LAPACKE_dstemr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize) isize;
extern fn LAPACKE_cstemr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize) isize;
extern fn LAPACKE_zstemr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize) isize;
pub const sstemr = LAPACKE_sstemr;
pub const dstemr = LAPACKE_dstemr;
pub const cstemr = LAPACKE_cstemr;
pub const zstemr = LAPACKE_zstemr;

extern fn LAPACKE_ssteqr(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsteqr(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_csteqr(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zsteqr(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const ssteqr = LAPACKE_ssteqr;
pub const dsteqr = LAPACKE_dsteqr;
pub const csteqr = LAPACKE_csteqr;
pub const zsteqr = LAPACKE_zsteqr;

extern fn LAPACKE_ssterf(n: isize, d: [*c]f32, e: [*c]f32) isize;
extern fn LAPACKE_dsterf(n: isize, d: [*c]f64, e: [*c]f64) isize;
pub const ssterf = LAPACKE_ssterf;
pub const dsterf = LAPACKE_dsterf;

extern fn LAPACKE_sstev(layout: c_int, jobz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dstev(layout: c_int, jobz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const sstev = LAPACKE_sstev;
pub const dstev = LAPACKE_dstev;

extern fn LAPACKE_sstevd(layout: c_int, jobz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dstevd(layout: c_int, jobz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const sstevd = LAPACKE_sstevd;
pub const dstevd = LAPACKE_dstevd;

extern fn LAPACKE_sstevr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_dstevr(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize) isize;
pub const sstevr = LAPACKE_sstevr;
pub const dstevr = LAPACKE_dstevr;

extern fn LAPACKE_sstevx(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dstevx(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const sstevx = LAPACKE_sstevx;
pub const dstevx = LAPACKE_dstevx;

extern fn LAPACKE_ssycon(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dsycon(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_csycon(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zsycon(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const ssycon = LAPACKE_ssycon;
pub const dsycon = LAPACKE_dsycon;
pub const csycon = LAPACKE_csycon;
pub const zsycon = LAPACKE_zsycon;

extern fn LAPACKE_ssyequb(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dsyequb(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_csyequb(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zsyequb(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const ssyequb = LAPACKE_ssyequb;
pub const dsyequb = LAPACKE_dsyequb;
pub const csyequb = LAPACKE_csyequb;
pub const zsyequb = LAPACKE_zsyequb;

extern fn LAPACKE_ssyev(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsyev(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64) isize;
pub const ssyev = LAPACKE_ssyev;
pub const dsyev = LAPACKE_dsyev;

extern fn LAPACKE_ssyevd(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsyevd(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64) isize;
pub const ssyevd = LAPACKE_ssyevd;
pub const dsyevd = LAPACKE_dsyevd;

extern fn LAPACKE_ssyevr(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_dsyevr(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize) isize;
pub const ssyevr = LAPACKE_ssyevr;
pub const dsyevr = LAPACKE_dsyevr;

extern fn LAPACKE_ssyevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsyevx(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const ssyevx = LAPACKE_ssyevx;
pub const dsyevx = LAPACKE_dsyevx;

extern fn LAPACKE_ssygst(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]const f32, ldb: isize) isize;
extern fn LAPACKE_dsygst(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]const f64, ldb: isize) isize;
pub const ssygst = LAPACKE_ssygst;
pub const dsygst = LAPACKE_dsygst;

extern fn LAPACKE_ssygv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsygv(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, w: [*c]f64) isize;
pub const ssygv = LAPACKE_ssygv;
pub const dsygv = LAPACKE_dsygv;

extern fn LAPACKE_ssygvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsygvd(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, w: [*c]f64) isize;
pub const ssygvd = LAPACKE_ssygvd;
pub const dsygvd = LAPACKE_dsygvd;

extern fn LAPACKE_ssygvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsygvx(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const ssygvx = LAPACKE_ssygvx;
pub const dsygvx = LAPACKE_dsygvx;

extern fn LAPACKE_ssyrfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dsyrfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_csyrfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zsyrfs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const ssyrfs = LAPACKE_ssyrfs;
pub const dsyrfs = LAPACKE_dsyrfs;
pub const csyrfs = LAPACKE_csyrfs;
pub const zsyrfs = LAPACKE_zsyrfs;

extern fn LAPACKE_ssyrfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dsyrfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_csyrfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zsyrfsx(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const ssyrfsx = LAPACKE_ssyrfsx;
pub const dsyrfsx = LAPACKE_dsyrfsx;
pub const csyrfsx = LAPACKE_csyrfsx;
pub const zsyrfsx = LAPACKE_zsyrfsx;

extern fn LAPACKE_ssysv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsysv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csysv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsysv(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const ssysv = LAPACKE_ssysv;
pub const dsysv = LAPACKE_dsysv;
pub const csysv = LAPACKE_csysv;
pub const zsysv = LAPACKE_zsysv;

extern fn LAPACKE_ssysvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dsysvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_csysvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_zsysvx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64) isize;
pub const ssysvx = LAPACKE_ssysvx;
pub const dsysvx = LAPACKE_dsysvx;
pub const csysvx = LAPACKE_csysvx;
pub const zsysvx = LAPACKE_zsysvx;

extern fn LAPACKE_ssysvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_dsysvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
extern fn LAPACKE_csysvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32) isize;
extern fn LAPACKE_zsysvxx(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64) isize;
pub const ssysvxx = LAPACKE_ssysvxx;
pub const dsysvxx = LAPACKE_dsysvxx;
pub const csysvxx = LAPACKE_csysvxx;
pub const zsysvxx = LAPACKE_zsysvxx;

extern fn LAPACKE_ssytrd(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, d: [*c]f32, e: [*c]f32, tau: [*c]f32) isize;
extern fn LAPACKE_dsytrd(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, d: [*c]f64, e: [*c]f64, tau: [*c]f64) isize;
pub const ssytrd = LAPACKE_ssytrd;
pub const dsytrd = LAPACKE_dsytrd;

extern fn LAPACKE_ssytrf(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dsytrf(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_csytrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zsytrf(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const ssytrf = LAPACKE_ssytrf;
pub const dsytrf = LAPACKE_dsytrf;
pub const csytrf = LAPACKE_csytrf;
pub const zsytrf = LAPACKE_zsytrf;

extern fn LAPACKE_ssytri(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_dsytri(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_csytri(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zsytri(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize) isize;
pub const ssytri = LAPACKE_ssytri;
pub const dsytri = LAPACKE_dsytri;
pub const csytri = LAPACKE_csytri;
pub const zsytri = LAPACKE_zsytri;

extern fn LAPACKE_ssytrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsytrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csytrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsytrs(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const ssytrs = LAPACKE_ssytrs;
pub const dsytrs = LAPACKE_dsytrs;
pub const csytrs = LAPACKE_csytrs;
pub const zsytrs = LAPACKE_zsytrs;

extern fn LAPACKE_stbcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const f32, ldab: isize, rcond: [*c]f32) isize;
extern fn LAPACKE_dtbcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const f64, ldab: isize, rcond: [*c]f64) isize;
extern fn LAPACKE_ctbcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const cf32, ldab: isize, rcond: [*c]f32) isize;
extern fn LAPACKE_ztbcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const cf64, ldab: isize, rcond: [*c]f64) isize;
pub const stbcon = LAPACKE_stbcon;
pub const dtbcon = LAPACKE_dtbcon;
pub const ctbcon = LAPACKE_ctbcon;
pub const ztbcon = LAPACKE_ztbcon;

extern fn LAPACKE_stbrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, b: [*c]const f32, ldb: isize, x: [*c]const f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dtbrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, b: [*c]const f64, ldb: isize, x: [*c]const f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_ctbrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, b: [*c]const cf32, ldb: isize, x: [*c]const cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_ztbrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, b: [*c]const cf64, ldb: isize, x: [*c]const cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const stbrfs = LAPACKE_stbrfs;
pub const dtbrfs = LAPACKE_dtbrfs;
pub const ctbrfs = LAPACKE_ctbrfs;
pub const ztbrfs = LAPACKE_ztbrfs;

extern fn LAPACKE_stbtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtbtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctbtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztbtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, b: [*c]cf64, ldb: isize) isize;
pub const stbtrs = LAPACKE_stbtrs;
pub const dtbtrs = LAPACKE_dtbtrs;
pub const ctbtrs = LAPACKE_ctbtrs;
pub const ztbtrs = LAPACKE_ztbtrs;

extern fn LAPACKE_stfsm(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: f32, a: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtfsm(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: f64, a: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctfsm(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: cf32, a: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztfsm(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: cf64, a: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const stfsm = LAPACKE_stfsm;
pub const dtfsm = LAPACKE_dtfsm;
pub const ctfsm = LAPACKE_ctfsm;
pub const ztfsm = LAPACKE_ztfsm;

extern fn LAPACKE_stftri(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]f32) isize;
extern fn LAPACKE_dtftri(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]f64) isize;
extern fn LAPACKE_ctftri(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]cf32) isize;
extern fn LAPACKE_ztftri(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]cf64) isize;
pub const stftri = LAPACKE_stftri;
pub const dtftri = LAPACKE_dtftri;
pub const ctftri = LAPACKE_ctftri;
pub const ztftri = LAPACKE_ztftri;

extern fn LAPACKE_stfttp(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f32, ap: [*c]f32) isize;
extern fn LAPACKE_dtfttp(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f64, ap: [*c]f64) isize;
extern fn LAPACKE_ctfttp(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf32, ap: [*c]cf32) isize;
extern fn LAPACKE_ztfttp(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf64, ap: [*c]cf64) isize;
pub const stfttp = LAPACKE_stfttp;
pub const dtfttp = LAPACKE_dtfttp;
pub const ctfttp = LAPACKE_ctfttp;
pub const ztfttp = LAPACKE_ztfttp;

extern fn LAPACKE_stfttr(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f32, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dtfttr(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f64, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_ctfttr(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf32, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_ztfttr(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf64, a: [*c]cf64, lda: isize) isize;
pub const stfttr = LAPACKE_stfttr;
pub const dtfttr = LAPACKE_dtfttr;
pub const ctfttr = LAPACKE_ctfttr;
pub const ztfttr = LAPACKE_ztfttr;

extern fn LAPACKE_stgevc(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const f32, lds: isize, p: [*c]const f32, ldp: isize, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_dtgevc(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const f64, lds: isize, p: [*c]const f64, ldp: isize, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ctgevc(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const cf32, lds: isize, p: [*c]const cf32, ldp: isize, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ztgevc(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const cf64, lds: isize, p: [*c]const cf64, ldp: isize, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, mm: isize, m: [*c]isize) isize;
pub const stgevc = LAPACKE_stgevc;
pub const dtgevc = LAPACKE_dtgevc;
pub const ctgevc = LAPACKE_ctgevc;
pub const ztgevc = LAPACKE_ztgevc;

extern fn LAPACKE_stgexc(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize, ifst: [*c]isize, ilst: [*c]isize) isize;
extern fn LAPACKE_dtgexc(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize, ifst: [*c]isize, ilst: [*c]isize) isize;
extern fn LAPACKE_ctgexc(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize, ifst: isize, ilst: isize) isize;
extern fn LAPACKE_ztgexc(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize, ifst: isize, ilst: isize) isize;
pub const stgexc = LAPACKE_stgexc;
pub const dtgexc = LAPACKE_dtgexc;
pub const ctgexc = LAPACKE_ctgexc;
pub const ztgexc = LAPACKE_ztgexc;

extern fn LAPACKE_stgsen(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize, m: [*c]isize, pl: [*c]f32, pr: [*c]f32, dif: [*c]f32) isize;
extern fn LAPACKE_dtgsen(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize, m: [*c]isize, pl: [*c]f64, pr: [*c]f64, dif: [*c]f64) isize;
extern fn LAPACKE_ctgsen(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize, m: [*c]isize, pl: [*c]f32, pr: [*c]f32, dif: [*c]f32) isize;
extern fn LAPACKE_ztgsen(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize, m: [*c]isize, pl: [*c]f64, pr: [*c]f64, dif: [*c]f64) isize;
pub const stgsen = LAPACKE_stgsen;
pub const dtgsen = LAPACKE_dtgsen;
pub const ctgsen = LAPACKE_ctgsen;
pub const ztgsen = LAPACKE_ztgsen;

extern fn LAPACKE_stgsja(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, tola: f32, tolb: f32, alpha: [*c]f32, beta: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, ncycle: [*c]isize) isize;
extern fn LAPACKE_dtgsja(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, tola: f64, tolb: f64, alpha: [*c]f64, beta: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, ncycle: [*c]isize) isize;
extern fn LAPACKE_ctgsja(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, tola: f32, tolb: f32, alpha: [*c]f32, beta: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, ncycle: [*c]isize) isize;
extern fn LAPACKE_ztgsja(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, tola: f64, tolb: f64, alpha: [*c]f64, beta: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, ncycle: [*c]isize) isize;
pub const stgsja = LAPACKE_stgsja;
pub const dtgsja = LAPACKE_dtgsja;
pub const ctgsja = LAPACKE_ctgsja;
pub const ztgsja = LAPACKE_ztgsja;

extern fn LAPACKE_stgsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, vl: [*c]const f32, ldvl: isize, vr: [*c]const f32, ldvr: isize, s: [*c]f32, dif: [*c]f32, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_dtgsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, vl: [*c]const f64, ldvl: isize, vr: [*c]const f64, ldvr: isize, s: [*c]f64, dif: [*c]f64, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ctgsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, vl: [*c]const cf32, ldvl: isize, vr: [*c]const cf32, ldvr: isize, s: [*c]f32, dif: [*c]f32, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ztgsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, vl: [*c]const cf64, ldvl: isize, vr: [*c]const cf64, ldvr: isize, s: [*c]f64, dif: [*c]f64, mm: isize, m: [*c]isize) isize;
pub const stgsna = LAPACKE_stgsna;
pub const dtgsna = LAPACKE_dtgsna;
pub const ctgsna = LAPACKE_ctgsna;
pub const ztgsna = LAPACKE_ztgsna;

extern fn LAPACKE_stgsyl(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]f32, ldc: isize, d: [*c]const f32, ldd: isize, e: [*c]const f32, lde: isize, f: [*c]f32, ldf: isize, scale: [*c]f32, dif: [*c]f32) isize;
extern fn LAPACKE_dtgsyl(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]f64, ldc: isize, d: [*c]const f64, ldd: isize, e: [*c]const f64, lde: isize, f: [*c]f64, ldf: isize, scale: [*c]f64, dif: [*c]f64) isize;
extern fn LAPACKE_ctgsyl(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize, d: [*c]const cf32, ldd: isize, e: [*c]const cf32, lde: isize, f: [*c]cf32, ldf: isize, scale: [*c]f32, dif: [*c]f32) isize;
extern fn LAPACKE_ztgsyl(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, d: [*c]const cf64, ldd: isize, e: [*c]const cf64, lde: isize, f: [*c]cf64, ldf: isize, scale: [*c]f64, dif: [*c]f64) isize;
pub const stgsyl = LAPACKE_stgsyl;
pub const dtgsyl = LAPACKE_dtgsyl;
pub const ctgsyl = LAPACKE_ctgsyl;
pub const ztgsyl = LAPACKE_ztgsyl;

extern fn LAPACKE_stpcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dtpcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const f64, rcond: [*c]f64) isize;
extern fn LAPACKE_ctpcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const cf32, rcond: [*c]f32) isize;
extern fn LAPACKE_ztpcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const cf64, rcond: [*c]f64) isize;
pub const stpcon = LAPACKE_stpcon;
pub const dtpcon = LAPACKE_dtpcon;
pub const ctpcon = LAPACKE_ctpcon;
pub const ztpcon = LAPACKE_ztpcon;

extern fn LAPACKE_stprfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]const f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dtprfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]const f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_ctprfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf32, b: [*c]const cf32, ldb: isize, x: [*c]const cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_ztprfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf64, b: [*c]const cf64, ldb: isize, x: [*c]const cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const stprfs = LAPACKE_stprfs;
pub const dtprfs = LAPACKE_dtprfs;
pub const ctprfs = LAPACKE_ctprfs;
pub const ztprfs = LAPACKE_ztprfs;

extern fn LAPACKE_stptri(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dtptri(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]f64) isize;
extern fn LAPACKE_ctptri(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_ztptri(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]cf64) isize;
pub const stptri = LAPACKE_stptri;
pub const dtptri = LAPACKE_dtptri;
pub const ctptri = LAPACKE_ctptri;
pub const ztptri = LAPACKE_ztptri;

extern fn LAPACKE_stptrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtptrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctptrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztptrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const stptrs = LAPACKE_stptrs;
pub const dtptrs = LAPACKE_dtptrs;
pub const ctptrs = LAPACKE_ctptrs;
pub const ztptrs = LAPACKE_ztptrs;

extern fn LAPACKE_stpttf(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const f32, arf: [*c]f32) isize;
extern fn LAPACKE_dtpttf(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const f64, arf: [*c]f64) isize;
extern fn LAPACKE_ctpttf(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const cf32, arf: [*c]cf32) isize;
extern fn LAPACKE_ztpttf(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const cf64, arf: [*c]cf64) isize;
pub const stpttf = LAPACKE_stpttf;
pub const dtpttf = LAPACKE_dtpttf;
pub const ctpttf = LAPACKE_ctpttf;
pub const ztpttf = LAPACKE_ztpttf;

extern fn LAPACKE_stpttr(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dtpttr(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_ctpttr(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_ztpttr(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, a: [*c]cf64, lda: isize) isize;
pub const stpttr = LAPACKE_stpttr;
pub const dtpttr = LAPACKE_dtpttr;
pub const ctpttr = LAPACKE_ctpttr;
pub const ztpttr = LAPACKE_ztpttr;

extern fn LAPACKE_strcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const f32, lda: isize, rcond: [*c]f32) isize;
extern fn LAPACKE_dtrcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const f64, lda: isize, rcond: [*c]f64) isize;
extern fn LAPACKE_ctrcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const cf32, lda: isize, rcond: [*c]f32) isize;
extern fn LAPACKE_ztrcon(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const cf64, lda: isize, rcond: [*c]f64) isize;
pub const strcon = LAPACKE_strcon;
pub const dtrcon = LAPACKE_dtrcon;
pub const ctrcon = LAPACKE_ctrcon;
pub const ztrcon = LAPACKE_ztrcon;

extern fn LAPACKE_strevc(layout: c_int, side: u8, howmny: u8, select: [*c]isize, n: isize, t: [*c]const f32, ldt: isize, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_dtrevc(layout: c_int, side: u8, howmny: u8, select: [*c]isize, n: isize, t: [*c]const f64, ldt: isize, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ctrevc(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]cf32, ldt: isize, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ztrevc(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]cf64, ldt: isize, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, mm: isize, m: [*c]isize) isize;
pub const strevc = LAPACKE_strevc;
pub const dtrevc = LAPACKE_dtrevc;
pub const ctrevc = LAPACKE_ctrevc;
pub const ztrevc = LAPACKE_ztrevc;

extern fn LAPACKE_strexc(layout: c_int, compq: u8, n: isize, t: [*c]f32, ldt: isize, q: [*c]f32, ldq: isize, ifst: [*c]isize, ilst: [*c]isize) isize;
extern fn LAPACKE_dtrexc(layout: c_int, compq: u8, n: isize, t: [*c]f64, ldt: isize, q: [*c]f64, ldq: isize, ifst: [*c]isize, ilst: [*c]isize) isize;
extern fn LAPACKE_ctrexc(layout: c_int, compq: u8, n: isize, t: [*c]cf32, ldt: isize, q: [*c]cf32, ldq: isize, ifst: isize, ilst: isize) isize;
extern fn LAPACKE_ztrexc(layout: c_int, compq: u8, n: isize, t: [*c]cf64, ldt: isize, q: [*c]cf64, ldq: isize, ifst: isize, ilst: isize) isize;
pub const strexc = LAPACKE_strexc;
pub const dtrexc = LAPACKE_dtrexc;
pub const ctrexc = LAPACKE_ctrexc;
pub const ztrexc = LAPACKE_ztrexc;

extern fn LAPACKE_strrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, x: [*c]const f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_dtrrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, x: [*c]const f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
extern fn LAPACKE_ctrrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, x: [*c]const cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32) isize;
extern fn LAPACKE_ztrrfs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, x: [*c]const cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64) isize;
pub const strrfs = LAPACKE_strrfs;
pub const dtrrfs = LAPACKE_dtrrfs;
pub const ctrrfs = LAPACKE_ctrrfs;
pub const ztrrfs = LAPACKE_ztrrfs;

extern fn LAPACKE_strsen(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]f32, ldt: isize, q: [*c]f32, ldq: isize, wr: [*c]f32, wi: [*c]f32, m: [*c]isize, s: [*c]f32, sep: [*c]f32) isize;
extern fn LAPACKE_dtrsen(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]f64, ldt: isize, q: [*c]f64, ldq: isize, wr: [*c]f64, wi: [*c]f64, m: [*c]isize, s: [*c]f64, sep: [*c]f64) isize;
extern fn LAPACKE_ctrsen(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]cf32, ldt: isize, q: [*c]cf32, ldq: isize, w: [*c]cf32, m: [*c]isize, s: [*c]f32, sep: [*c]f32) isize;
extern fn LAPACKE_ztrsen(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]cf64, ldt: isize, q: [*c]cf64, ldq: isize, w: [*c]cf64, m: [*c]isize, s: [*c]f64, sep: [*c]f64) isize;
pub const strsen = LAPACKE_strsen;
pub const dtrsen = LAPACKE_dtrsen;
pub const ctrsen = LAPACKE_ctrsen;
pub const ztrsen = LAPACKE_ztrsen;

extern fn LAPACKE_strsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const f32, ldt: isize, vl: [*c]const f32, ldvl: isize, vr: [*c]const f32, ldvr: isize, s: [*c]f32, sep: [*c]f32, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_dtrsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const f64, ldt: isize, vl: [*c]const f64, ldvl: isize, vr: [*c]const f64, ldvr: isize, s: [*c]f64, sep: [*c]f64, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ctrsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const cf32, ldt: isize, vl: [*c]const cf32, ldvl: isize, vr: [*c]const cf32, ldvr: isize, s: [*c]f32, sep: [*c]f32, mm: isize, m: [*c]isize) isize;
extern fn LAPACKE_ztrsna(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const cf64, ldt: isize, vl: [*c]const cf64, ldvl: isize, vr: [*c]const cf64, ldvr: isize, s: [*c]f64, sep: [*c]f64, mm: isize, m: [*c]isize) isize;
pub const strsna = LAPACKE_strsna;
pub const dtrsna = LAPACKE_dtrsna;
pub const ctrsna = LAPACKE_ctrsna;
pub const ztrsna = LAPACKE_ztrsna;

extern fn LAPACKE_strsyl(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]f32, ldc: isize, scale: [*c]f32) isize;
extern fn LAPACKE_dtrsyl(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]f64, ldc: isize, scale: [*c]f64) isize;
extern fn LAPACKE_ctrsyl(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize, scale: [*c]f32) isize;
extern fn LAPACKE_ztrsyl(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, scale: [*c]f64) isize;
pub const strsyl = LAPACKE_strsyl;
pub const dtrsyl = LAPACKE_dtrsyl;
pub const ctrsyl = LAPACKE_ctrsyl;
pub const ztrsyl = LAPACKE_ztrsyl;

extern fn LAPACKE_strsyl3(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]f32, ldc: isize, scale: [*c]f32) isize;
extern fn LAPACKE_dtrsyl3(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]f64, ldc: isize, scale: [*c]f64) isize;
extern fn LAPACKE_ztrsyl3(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, scale: [*c]f64) isize;
pub const strsyl3 = LAPACKE_strsyl3;
pub const dtrsyl3 = LAPACKE_dtrsyl3;
pub const ztrsyl3 = LAPACKE_ztrsyl3;

extern fn LAPACKE_strtri(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dtrtri(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_ctrtri(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_ztrtri(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const strtri = LAPACKE_strtri;
pub const dtrtri = LAPACKE_dtrtri;
pub const ctrtri = LAPACKE_ctrtri;
pub const ztrtri = LAPACKE_ztrtri;

extern fn LAPACKE_strtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtrtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctrtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztrtrs(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const strtrs = LAPACKE_strtrs;
pub const dtrtrs = LAPACKE_dtrtrs;
pub const ctrtrs = LAPACKE_ctrtrs;
pub const ztrtrs = LAPACKE_ztrtrs;

extern fn LAPACKE_strttf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const f32, lda: isize, arf: [*c]f32) isize;
extern fn LAPACKE_dtrttf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const f64, lda: isize, arf: [*c]f64) isize;
extern fn LAPACKE_ctrttf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, arf: [*c]cf32) isize;
extern fn LAPACKE_ztrttf(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, arf: [*c]cf64) isize;
pub const strttf = LAPACKE_strttf;
pub const dtrttf = LAPACKE_dtrttf;
pub const ctrttf = LAPACKE_ctrttf;
pub const ztrttf = LAPACKE_ztrttf;

extern fn LAPACKE_strttp(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dtrttp(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, ap: [*c]f64) isize;
extern fn LAPACKE_ctrttp(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_ztrttp(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, ap: [*c]cf64) isize;
pub const strttp = LAPACKE_strttp;
pub const dtrttp = LAPACKE_dtrttp;
pub const ctrttp = LAPACKE_ctrttp;
pub const ztrttp = LAPACKE_ztrttp;

extern fn LAPACKE_stzrzf(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dtzrzf(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64) isize;
extern fn LAPACKE_ctzrzf(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_ztzrzf(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64) isize;
pub const stzrzf = LAPACKE_stzrzf;
pub const dtzrzf = LAPACKE_dtzrzf;
pub const ctzrzf = LAPACKE_ctzrzf;
pub const ztzrzf = LAPACKE_ztzrzf;

extern fn LAPACKE_cungbr(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zungbr(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cungbr = LAPACKE_cungbr;
pub const zungbr = LAPACKE_zungbr;

extern fn LAPACKE_cunghr(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zunghr(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cunghr = LAPACKE_cunghr;
pub const zunghr = LAPACKE_zunghr;

extern fn LAPACKE_cunglq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zunglq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cunglq = LAPACKE_cunglq;
pub const zunglq = LAPACKE_zunglq;

extern fn LAPACKE_cungql(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zungql(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cungql = LAPACKE_cungql;
pub const zungql = LAPACKE_zungql;

extern fn LAPACKE_cungqr(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zungqr(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cungqr = LAPACKE_cungqr;
pub const zungqr = LAPACKE_zungqr;

extern fn LAPACKE_cungrq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zungrq(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cungrq = LAPACKE_cungrq;
pub const zungrq = LAPACKE_zungrq;

extern fn LAPACKE_cungtr(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32) isize;
extern fn LAPACKE_zungtr(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64) isize;
pub const cungtr = LAPACKE_cungtr;
pub const zungtr = LAPACKE_zungtr;

extern fn LAPACKE_cungtsqr_row(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]cf32, lda: isize, t: [*c]const cf32, ldt: isize) isize;
extern fn LAPACKE_zungtsqr_row(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]cf64, lda: isize, t: [*c]const cf64, ldt: isize) isize;
pub const cungtsqr_row = LAPACKE_cungtsqr_row;
pub const zungtsqr_row = LAPACKE_zungtsqr_row;

extern fn LAPACKE_cunmbr(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmbr(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmbr = LAPACKE_cunmbr;
pub const zunmbr = LAPACKE_zunmbr;

extern fn LAPACKE_cunmhr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmhr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmhr = LAPACKE_cunmhr;
pub const zunmhr = LAPACKE_zunmhr;

extern fn LAPACKE_cunmlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmlq = LAPACKE_cunmlq;
pub const zunmlq = LAPACKE_zunmlq;

extern fn LAPACKE_cunmql(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmql(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmql = LAPACKE_cunmql;
pub const zunmql = LAPACKE_zunmql;

extern fn LAPACKE_cunmqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmqr = LAPACKE_cunmqr;
pub const zunmqr = LAPACKE_zunmqr;

extern fn LAPACKE_cunmrq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmrq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmrq = LAPACKE_cunmrq;
pub const zunmrq = LAPACKE_zunmrq;

extern fn LAPACKE_cunmrz(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmrz(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmrz = LAPACKE_cunmrz;
pub const zunmrz = LAPACKE_zunmrz;

extern fn LAPACKE_cunmtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zunmtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cunmtr = LAPACKE_cunmtr;
pub const zunmtr = LAPACKE_zunmtr;

extern fn LAPACKE_cupgtr(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, tau: [*c]const cf32, q: [*c]cf32, ldq: isize) isize;
extern fn LAPACKE_zupgtr(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, tau: [*c]const cf64, q: [*c]cf64, ldq: isize) isize;
pub const cupgtr = LAPACKE_cupgtr;
pub const zupgtr = LAPACKE_zupgtr;

extern fn LAPACKE_cupmtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const cf32, tau: [*c]const cf32, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zupmtr(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const cf64, tau: [*c]const cf64, c: [*c]cf64, ldc: isize) isize;
pub const cupmtr = LAPACKE_cupmtr;
pub const zupmtr = LAPACKE_zupmtr;

extern fn LAPACKE_sbdsdc_work(layout: c_int, uplo: u8, compq: u8, n: isize, d: [*c]f32, e: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, q: [*c]f32, iq: [*c]isize, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dbdsdc_work(layout: c_int, uplo: u8, compq: u8, n: isize, d: [*c]f64, e: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, q: [*c]f64, iq: [*c]isize, work: [*c]f64, iwork: [*c]isize) isize;
pub const sbdsdc_work = LAPACKE_sbdsdc_work;
pub const dbdsdc_work = LAPACKE_dbdsdc_work;

extern fn LAPACKE_sbdsvdx_work(layout: c_int, uplo: u8, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, ns: [*c]isize, s: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dbdsvdx_work(layout: c_int, uplo: u8, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, ns: [*c]isize, s: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize) isize;
pub const sbdsvdx_work = LAPACKE_sbdsvdx_work;
pub const dbdsvdx_work = LAPACKE_dbdsvdx_work;

extern fn LAPACKE_sbdsqr_work(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f32, e: [*c]f32, vt: [*c]f32, ldvt: isize, u: [*c]f32, ldu: isize, c: [*c]f32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_dbdsqr_work(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f64, e: [*c]f64, vt: [*c]f64, ldvt: isize, u: [*c]f64, ldu: isize, c: [*c]f64, ldc: isize, work: [*c]f64) isize;
extern fn LAPACKE_cbdsqr_work(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f32, e: [*c]f32, vt: [*c]cf32, ldvt: isize, u: [*c]cf32, ldu: isize, c: [*c]cf32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_zbdsqr_work(layout: c_int, uplo: u8, n: isize, ncvt: isize, nru: isize, ncc: isize, d: [*c]f64, e: [*c]f64, vt: [*c]cf64, ldvt: isize, u: [*c]cf64, ldu: isize, c: [*c]cf64, ldc: isize, work: [*c]f64) isize;
pub const sbdsqr_work = LAPACKE_sbdsqr_work;
pub const dbdsqr_work = LAPACKE_dbdsqr_work;
pub const cbdsqr_work = LAPACKE_cbdsqr_work;
pub const zbdsqr_work = LAPACKE_zbdsqr_work;

extern fn LAPACKE_sdisna_work(job: u8, m: isize, n: isize, d: [*c]const f32, sep: [*c]f32) isize;
extern fn LAPACKE_ddisna_work(job: u8, m: isize, n: isize, d: [*c]const f64, sep: [*c]f64) isize;
pub const sdisna_work = LAPACKE_sdisna_work;
pub const ddisna_work = LAPACKE_ddisna_work;

extern fn LAPACKE_sgbbrd_work(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]f32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]f32, ldq: isize, pt: [*c]f32, ldpt: isize, c: [*c]f32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_dgbbrd_work(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]f64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]f64, ldq: isize, pt: [*c]f64, ldpt: isize, c: [*c]f64, ldc: isize, work: [*c]f64) isize;
extern fn LAPACKE_cgbbrd_work(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]cf32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]cf32, ldq: isize, pt: [*c]cf32, ldpt: isize, c: [*c]cf32, ldc: isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgbbrd_work(layout: c_int, vect: u8, m: isize, n: isize, ncc: isize, kl: isize, ku: isize, ab: [*c]cf64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]cf64, ldq: isize, pt: [*c]cf64, ldpt: isize, c: [*c]cf64, ldc: isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgbbrd_work = LAPACKE_sgbbrd_work;
pub const dgbbrd_work = LAPACKE_dgbbrd_work;
pub const cgbbrd_work = LAPACKE_cgbbrd_work;
pub const zgbbrd_work = LAPACKE_zgbbrd_work;

extern fn LAPACKE_sgbcon_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgbcon_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgbcon_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgbcon_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgbcon_work = LAPACKE_sgbcon_work;
pub const dgbcon_work = LAPACKE_dgbcon_work;
pub const cgbcon_work = LAPACKE_cgbcon_work;
pub const zgbcon_work = LAPACKE_zgbcon_work;

extern fn LAPACKE_sgbequ_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgbequ_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgbequ_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgbequ_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgbequ_work = LAPACKE_sgbequ_work;
pub const dgbequ_work = LAPACKE_dgbequ_work;
pub const cgbequ_work = LAPACKE_cgbequ_work;
pub const zgbequ_work = LAPACKE_zgbequ_work;

extern fn LAPACKE_sgbequb_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgbequb_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgbequb_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgbequb_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgbequb_work = LAPACKE_sgbequb_work;
pub const dgbequb_work = LAPACKE_dgbequb_work;
pub const cgbequb_work = LAPACKE_cgbequb_work;
pub const zgbequb_work = LAPACKE_zgbequb_work;

extern fn LAPACKE_sgbrfs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, afb: [*c]const f32, ldafb: isize, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgbrfs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, afb: [*c]const f64, ldafb: isize, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgbrfs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, afb: [*c]const cf32, ldafb: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgbrfs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, afb: [*c]const cf64, ldafb: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgbrfs_work = LAPACKE_sgbrfs_work;
pub const dgbrfs_work = LAPACKE_dgbrfs_work;
pub const cgbrfs_work = LAPACKE_cgbrfs_work;
pub const zgbrfs_work = LAPACKE_zgbrfs_work;

extern fn LAPACKE_sgbrfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, afb: [*c]const f32, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgbrfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, afb: [*c]const f64, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgbrfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, afb: [*c]const cf32, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgbrfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, afb: [*c]const cf64, ldafb: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgbrfsx_work = LAPACKE_sgbrfsx_work;
pub const dgbrfsx_work = LAPACKE_dgbrfsx_work;
pub const cgbrfsx_work = LAPACKE_cgbrfsx_work;
pub const zgbrfsx_work = LAPACKE_zgbrfsx_work;

extern fn LAPACKE_sgbsv_work(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f32, ldab: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgbsv_work(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f64, ldab: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgbsv_work(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgbsv_work(layout: c_int, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const sgbsv_work = LAPACKE_sgbsv_work;
pub const dgbsv_work = LAPACKE_dgbsv_work;
pub const cgbsv_work = LAPACKE_cgbsv_work;
pub const zgbsv_work = LAPACKE_zgbsv_work;

extern fn LAPACKE_sgbsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f32, ldab: isize, afb: [*c]f32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgbsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f64, ldab: isize, afb: [*c]f64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgbsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, afb: [*c]cf32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgbsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, afb: [*c]cf64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgbsvx_work = LAPACKE_sgbsvx_work;
pub const dgbsvx_work = LAPACKE_dgbsvx_work;
pub const cgbsvx_work = LAPACKE_cgbsvx_work;
pub const zgbsvx_work = LAPACKE_zgbsvx_work;

extern fn LAPACKE_sgbsvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f32, ldab: isize, afb: [*c]f32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgbsvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]f64, ldab: isize, afb: [*c]f64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgbsvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, afb: [*c]cf32, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgbsvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, afb: [*c]cf64, ldafb: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgbsvxx_work = LAPACKE_sgbsvxx_work;
pub const dgbsvxx_work = LAPACKE_dgbsvxx_work;
pub const cgbsvxx_work = LAPACKE_cgbsvxx_work;
pub const zgbsvxx_work = LAPACKE_zgbsvxx_work;

extern fn LAPACKE_sgbtrf_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]f32, ldab: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgbtrf_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]f64, ldab: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgbtrf_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]cf32, ldab: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgbtrf_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, ab: [*c]cf64, ldab: isize, ipiv: [*c]isize) isize;
pub const sgbtrf_work = LAPACKE_sgbtrf_work;
pub const dgbtrf_work = LAPACKE_dgbtrf_work;
pub const cgbtrf_work = LAPACKE_cgbtrf_work;
pub const zgbtrf_work = LAPACKE_zgbtrf_work;

extern fn LAPACKE_sgbtrs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgbtrs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgbtrs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgbtrs_work(layout: c_int, trans: u8, n: isize, kl: isize, ku: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const sgbtrs_work = LAPACKE_sgbtrs_work;
pub const dgbtrs_work = LAPACKE_dgbtrs_work;
pub const cgbtrs_work = LAPACKE_cgbtrs_work;
pub const zgbtrs_work = LAPACKE_zgbtrs_work;

extern fn LAPACKE_sgebak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f32, m: isize, v: [*c]f32, ldv: isize) isize;
extern fn LAPACKE_dgebak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f64, m: isize, v: [*c]f64, ldv: isize) isize;
extern fn LAPACKE_cgebak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f32, m: isize, v: [*c]cf32, ldv: isize) isize;
extern fn LAPACKE_zgebak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, scale: [*c]const f64, m: isize, v: [*c]cf64, ldv: isize) isize;
pub const sgebak_work = LAPACKE_sgebak_work;
pub const dgebak_work = LAPACKE_dgebak_work;
pub const cgebak_work = LAPACKE_cgebak_work;
pub const zgebak_work = LAPACKE_zgebak_work;

extern fn LAPACKE_sgebal_work(layout: c_int, job: u8, n: isize, a: [*c]f32, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32) isize;
extern fn LAPACKE_dgebal_work(layout: c_int, job: u8, n: isize, a: [*c]f64, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64) isize;
extern fn LAPACKE_cgebal_work(layout: c_int, job: u8, n: isize, a: [*c]cf32, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32) isize;
extern fn LAPACKE_zgebal_work(layout: c_int, job: u8, n: isize, a: [*c]cf64, lda: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64) isize;
pub const sgebal_work = LAPACKE_sgebal_work;
pub const dgebal_work = LAPACKE_dgebal_work;
pub const cgebal_work = LAPACKE_cgebal_work;
pub const zgebal_work = LAPACKE_zgebal_work;

extern fn LAPACKE_sgebrd_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, d: [*c]f32, e: [*c]f32, tauq: [*c]f32, taup: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgebrd_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, d: [*c]f64, e: [*c]f64, tauq: [*c]f64, taup: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgebrd_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, d: [*c]f32, e: [*c]f32, tauq: [*c]cf32, taup: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgebrd_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, d: [*c]f64, e: [*c]f64, tauq: [*c]cf64, taup: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgebrd_work = LAPACKE_sgebrd_work;
pub const dgebrd_work = LAPACKE_dgebrd_work;
pub const cgebrd_work = LAPACKE_cgebrd_work;
pub const zgebrd_work = LAPACKE_zgebrd_work;

extern fn LAPACKE_sgecon_work(layout: c_int, norm: u8, n: isize, a: [*c]const f32, lda: isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgecon_work(layout: c_int, norm: u8, n: isize, a: [*c]const f64, lda: isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgecon_work(layout: c_int, norm: u8, n: isize, a: [*c]const cf32, lda: isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgecon_work(layout: c_int, norm: u8, n: isize, a: [*c]const cf64, lda: isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgecon_work = LAPACKE_sgecon_work;
pub const dgecon_work = LAPACKE_dgecon_work;
pub const cgecon_work = LAPACKE_cgecon_work;
pub const zgecon_work = LAPACKE_zgecon_work;

extern fn LAPACKE_sgeequ_work(layout: c_int, m: isize, n: isize, a: [*c]const f32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgeequ_work(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgeequ_work(layout: c_int, m: isize, n: isize, a: [*c]const cf32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgeequ_work(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgeequ_work = LAPACKE_sgeequ_work;
pub const dgeequ_work = LAPACKE_dgeequ_work;
pub const cgeequ_work = LAPACKE_cgeequ_work;
pub const zgeequ_work = LAPACKE_zgeequ_work;

extern fn LAPACKE_sgeequb_work(layout: c_int, m: isize, n: isize, a: [*c]const f32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dgeequb_work(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cgeequb_work(layout: c_int, m: isize, n: isize, a: [*c]const cf32, lda: isize, r: [*c]f32, c: [*c]f32, rowcnd: [*c]f32, colcnd: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zgeequb_work(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, r: [*c]f64, c: [*c]f64, rowcnd: [*c]f64, colcnd: [*c]f64, amax: [*c]f64) isize;
pub const sgeequb_work = LAPACKE_sgeequb_work;
pub const dgeequb_work = LAPACKE_dgeequb_work;
pub const cgeequb_work = LAPACKE_cgeequb_work;
pub const zgeequb_work = LAPACKE_zgeequb_work;

extern fn LAPACKE_sgees_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f32, [*c]const f32) isize, n: isize, a: [*c]f32, lda: isize, sdim: [*c]isize, wr: [*c]f32, wi: [*c]f32, vs: [*c]f32, ldvs: isize, work: [*c]f32, lwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_dgees_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f64, [*c]const f64) isize, n: isize, a: [*c]f64, lda: isize, sdim: [*c]isize, wr: [*c]f64, wi: [*c]f64, vs: [*c]f64, ldvs: isize, work: [*c]f64, lwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_cgees_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf32) isize, n: isize, a: [*c]cf32, lda: isize, sdim: [*c]isize, w: [*c]cf32, vs: [*c]cf32, ldvs: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, bwork: [*c]isize) isize;
extern fn LAPACKE_zgees_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf64) isize, n: isize, a: [*c]cf64, lda: isize, sdim: [*c]isize, w: [*c]cf64, vs: [*c]cf64, ldvs: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, bwork: [*c]isize) isize;
pub const sgees_work = LAPACKE_sgees_work;
pub const dgees_work = LAPACKE_dgees_work;
pub const cgees_work = LAPACKE_cgees_work;
pub const zgees_work = LAPACKE_zgees_work;

extern fn LAPACKE_sgeesx_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f32, [*c]const f32) isize, sense: u8, n: isize, a: [*c]f32, lda: isize, sdim: [*c]isize, wr: [*c]f32, wi: [*c]f32, vs: [*c]f32, ldvs: isize, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_dgeesx_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const f64, [*c]const f64) isize, sense: u8, n: isize, a: [*c]f64, lda: isize, sdim: [*c]isize, wr: [*c]f64, wi: [*c]f64, vs: [*c]f64, ldvs: isize, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_cgeesx_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf32) isize, sense: u8, n: isize, a: [*c]cf32, lda: isize, sdim: [*c]isize, w: [*c]cf32, vs: [*c]cf32, ldvs: isize, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32, bwork: [*c]isize) isize;
extern fn LAPACKE_zgeesx_work(layout: c_int, jobvs: u8, sort: u8, select: *const fn ([*c]const cf64) isize, sense: u8, n: isize, a: [*c]cf64, lda: isize, sdim: [*c]isize, w: [*c]cf64, vs: [*c]cf64, ldvs: isize, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64, bwork: [*c]isize) isize;
pub const sgeesx_work = LAPACKE_sgeesx_work;
pub const dgeesx_work = LAPACKE_dgeesx_work;
pub const cgeesx_work = LAPACKE_cgeesx_work;
pub const zgeesx_work = LAPACKE_zgeesx_work;

extern fn LAPACKE_sgeev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f32, lda: isize, wr: [*c]f32, wi: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgeev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f64, lda: isize, wr: [*c]f64, wi: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgeev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zgeev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sgeev_work = LAPACKE_sgeev_work;
pub const dgeev_work = LAPACKE_dgeev_work;
pub const cgeev_work = LAPACKE_cgeev_work;
pub const zgeev_work = LAPACKE_zgeev_work;

extern fn LAPACKE_sgeevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f32, lda: isize, wr: [*c]f32, wi: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32, abnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dgeevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f64, lda: isize, wr: [*c]f64, wi: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64, abnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cgeevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f32, abnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zgeevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, scale: [*c]f64, abnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sgeevx_work = LAPACKE_sgeevx_work;
pub const dgeevx_work = LAPACKE_dgeevx_work;
pub const cgeevx_work = LAPACKE_cgeevx_work;
pub const zgeevx_work = LAPACKE_zgeevx_work;

extern fn LAPACKE_sgehrd_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgehrd_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgehrd_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgehrd_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgehrd_work = LAPACKE_sgehrd_work;
pub const dgehrd_work = LAPACKE_dgehrd_work;
pub const cgehrd_work = LAPACKE_cgehrd_work;
pub const zgehrd_work = LAPACKE_zgehrd_work;

extern fn LAPACKE_sgejsv_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]f32, lda: isize, sva: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dgejsv_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]f64, lda: isize, sva: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cgejsv_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, sva: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, cwork: [*c]cf32, lwork: isize, work: [*c]f32, lrwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_zgejsv_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, jobr: u8, jobt: u8, jobp: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, sva: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, cwork: [*c]cf64, lwork: isize, work: [*c]f64, lrwork: isize, iwork: [*c]isize) isize;
pub const sgejsv_work = LAPACKE_sgejsv_work;
pub const dgejsv_work = LAPACKE_dgejsv_work;
pub const cgejsv_work = LAPACKE_cgejsv_work;
pub const zgejsv_work = LAPACKE_zgejsv_work;

extern fn LAPACKE_sgelq2_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dgelq2_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cgelq2_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32) isize;
extern fn LAPACKE_zgelq2_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64) isize;
pub const sgelq2_work = LAPACKE_sgelq2_work;
pub const dgelq2_work = LAPACKE_dgelq2_work;
pub const cgelq2_work = LAPACKE_cgelq2_work;
pub const zgelq2_work = LAPACKE_zgelq2_work;

extern fn LAPACKE_sgelqf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgelqf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgelqf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgelqf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgelqf_work = LAPACKE_sgelqf_work;
pub const dgelqf_work = LAPACKE_dgelqf_work;
pub const cgelqf_work = LAPACKE_cgelqf_work;
pub const zgelqf_work = LAPACKE_zgelqf_work;

extern fn LAPACKE_sgels_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgels_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgels_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgels_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgels_work = LAPACKE_sgels_work;
pub const dgels_work = LAPACKE_dgels_work;
pub const cgels_work = LAPACKE_cgels_work;
pub const zgels_work = LAPACKE_zgels_work;

extern fn LAPACKE_sgelsd_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dgelsd_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cgelsd_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_zgelsd_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize) isize;
pub const sgelsd_work = LAPACKE_sgelsd_work;
pub const dgelsd_work = LAPACKE_dgelsd_work;
pub const cgelsd_work = LAPACKE_cgelsd_work;
pub const zgelsd_work = LAPACKE_zgelsd_work;

extern fn LAPACKE_sgelss_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgelss_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgelss_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, s: [*c]f32, rcond: f32, rank: [*c]isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zgelss_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, s: [*c]f64, rcond: f64, rank: [*c]isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sgelss_work = LAPACKE_sgelss_work;
pub const dgelss_work = LAPACKE_dgelss_work;
pub const cgelss_work = LAPACKE_cgelss_work;
pub const zgelss_work = LAPACKE_zgelss_work;

extern fn LAPACKE_sgelsy_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, jpvt: [*c]isize, rcond: f32, rank: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgelsy_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, jpvt: [*c]isize, rcond: f64, rank: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgelsy_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, jpvt: [*c]isize, rcond: f32, rank: [*c]isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zgelsy_work(layout: c_int, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, jpvt: [*c]isize, rcond: f64, rank: [*c]isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sgelsy_work = LAPACKE_sgelsy_work;
pub const dgelsy_work = LAPACKE_dgelsy_work;
pub const cgelsy_work = LAPACKE_cgelsy_work;
pub const zgelsy_work = LAPACKE_zgelsy_work;

extern fn LAPACKE_sgeqlf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgeqlf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgeqlf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgeqlf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgeqlf_work = LAPACKE_sgeqlf_work;
pub const dgeqlf_work = LAPACKE_dgeqlf_work;
pub const cgeqlf_work = LAPACKE_cgeqlf_work;
pub const zgeqlf_work = LAPACKE_zgeqlf_work;

extern fn LAPACKE_sgeqp3_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, jpvt: [*c]isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgeqp3_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, jpvt: [*c]isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgeqp3_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, jpvt: [*c]isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zgeqp3_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, jpvt: [*c]isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sgeqp3_work = LAPACKE_sgeqp3_work;
pub const dgeqp3_work = LAPACKE_dgeqp3_work;
pub const cgeqp3_work = LAPACKE_cgeqp3_work;
pub const zgeqp3_work = LAPACKE_zgeqp3_work;

extern fn LAPACKE_sgeqpf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, jpvt: [*c]isize, tau: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dgeqpf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, jpvt: [*c]isize, tau: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cgeqpf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, jpvt: [*c]isize, tau: [*c]cf32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgeqpf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, jpvt: [*c]isize, tau: [*c]cf64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgeqpf_work = LAPACKE_sgeqpf_work;
pub const dgeqpf_work = LAPACKE_dgeqpf_work;
pub const cgeqpf_work = LAPACKE_cgeqpf_work;
pub const zgeqpf_work = LAPACKE_zgeqpf_work;

extern fn LAPACKE_sgeqr2_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dgeqr2_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cgeqr2_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32) isize;
extern fn LAPACKE_zgeqr2_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64) isize;
pub const sgeqr2_work = LAPACKE_sgeqr2_work;
pub const dgeqr2_work = LAPACKE_dgeqr2_work;
pub const cgeqr2_work = LAPACKE_cgeqr2_work;
pub const zgeqr2_work = LAPACKE_zgeqr2_work;

extern fn LAPACKE_sgeqrf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgeqrf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgeqrf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgeqrf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgeqrf_work = LAPACKE_sgeqrf_work;
pub const dgeqrf_work = LAPACKE_dgeqrf_work;
pub const cgeqrf_work = LAPACKE_cgeqrf_work;
pub const zgeqrf_work = LAPACKE_zgeqrf_work;

extern fn LAPACKE_sgeqrfp_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgeqrfp_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgeqrfp_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgeqrfp_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgeqrfp_work = LAPACKE_sgeqrfp_work;
pub const dgeqrfp_work = LAPACKE_dgeqrfp_work;
pub const cgeqrfp_work = LAPACKE_cgeqrfp_work;
pub const zgeqrfp_work = LAPACKE_zgeqrfp_work;

extern fn LAPACKE_sgerfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgerfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgerfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgerfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgerfs_work = LAPACKE_sgerfs_work;
pub const dgerfs_work = LAPACKE_dgerfs_work;
pub const cgerfs_work = LAPACKE_cgerfs_work;
pub const zgerfs_work = LAPACKE_zgerfs_work;

extern fn LAPACKE_sgerfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgerfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgerfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f32, c: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgerfsx_work(layout: c_int, trans: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, r: [*c]const f64, c: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgerfsx_work = LAPACKE_sgerfsx_work;
pub const dgerfsx_work = LAPACKE_dgerfsx_work;
pub const cgerfsx_work = LAPACKE_cgerfsx_work;
pub const zgerfsx_work = LAPACKE_zgerfsx_work;

extern fn LAPACKE_sgerqf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgerqf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgerqf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgerqf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgerqf_work = LAPACKE_sgerqf_work;
pub const dgerqf_work = LAPACKE_dgerqf_work;
pub const cgerqf_work = LAPACKE_cgerqf_work;
pub const zgerqf_work = LAPACKE_zgerqf_work;

extern fn LAPACKE_sgesdd_work(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]f32, lda: isize, s: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dgesdd_work(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]f64, lda: isize, s: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cgesdd_work(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, s: [*c]f32, u: [*c]cf32, ldu: isize, vt: [*c]cf32, ldvt: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_zgesdd_work(layout: c_int, jobz: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, s: [*c]f64, u: [*c]cf64, ldu: isize, vt: [*c]cf64, ldvt: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize) isize;
pub const sgesdd_work = LAPACKE_sgesdd_work;
pub const dgesdd_work = LAPACKE_dgesdd_work;
pub const cgesdd_work = LAPACKE_cgesdd_work;
pub const zgesdd_work = LAPACKE_zgesdd_work;

extern fn LAPACKE_sgedmd_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, x: [*c]f32, ldx: isize, y: [*c]f32, ldy: isize, nrnk: isize, tol: [*c]f32, k: isize, reig: [*c]f32, imeig: [*c]f32, z: [*c]f32, ldz: isize, res: [*c]f32, b: [*c]f32, ldb: isize, w: [*c]f32, ldw: isize, s: [*c]f32, lds: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dgedmd_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, x: [*c]f64, ldx: isize, y: [*c]f64, ldy: isize, nrnk: isize, tol: [*c]f64, k: isize, reig: [*c]f64, imeig: [*c]f64, z: [*c]f64, ldz: isize, res: [*c]f64, b: [*c]f64, ldb: isize, w: [*c]f64, ldw: isize, s: [*c]f64, lds: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_cgedmd_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, x: [*c]cf32, ldx: isize, y: [*c]cf32, ldy: isize, nrnk: isize, tol: [*c]f32, k: isize, eigs: [*c]cf32, z: [*c]cf32, ldz: isize, res: [*c]f32, b: [*c]cf32, ldb: isize, w: [*c]cf32, ldw: isize, s: [*c]cf32, lds: isize, zwork: [*c]cf32, lzwork: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zgedmd_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, x: [*c]cf64, ldx: isize, y: [*c]cf64, ldy: isize, nrnk: isize, tol: [*c]f64, k: isize, eigs: [*c]cf64, z: [*c]cf64, ldz: isize, res: [*c]f64, b: [*c]cf64, ldb: isize, w: [*c]cf64, ldw: isize, s: [*c]cf64, lds: isize, zwork: [*c]cf64, lzwork: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sgedmd_work = LAPACKE_sgedmd_work;
pub const dgedmd_work = LAPACKE_dgedmd_work;
pub const cgedmd_work = LAPACKE_cgedmd_work;
pub const zgedmd_work = LAPACKE_zgedmd_work;

extern fn LAPACKE_sgedmdq_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobq: u8, jobt: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, f: [*c]f32, ldf: isize, x: [*c]f32, ldx: isize, y: [*c]f32, ldy: isize, nrnk: isize, tol: [*c]f32, k: isize, reig: [*c]f32, imeig: [*c]f32, z: [*c]f32, ldz: isize, res: [*c]f32, b: [*c]f32, ldb: isize, v: [*c]f32, ldv: isize, s: [*c]f32, lds: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dgedmdq_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobq: u8, jobt: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, f: [*c]f64, ldf: isize, x: [*c]f64, ldx: isize, y: [*c]f64, ldy: isize, nrnk: isize, tol: [*c]f64, k: isize, reig: [*c]f64, imeig: [*c]f64, z: [*c]f64, ldz: isize, res: [*c]f64, b: [*c]f64, ldb: isize, v: [*c]f64, ldv: isize, s: [*c]f64, lds: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_cgedmdq_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobq: u8, jobt: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, f: [*c]cf32, ldf: isize, x: [*c]cf32, ldx: isize, y: [*c]cf32, ldy: isize, nrnk: isize, tol: [*c]f32, k: isize, eigs: [*c]cf32, z: [*c]cf32, ldz: isize, res: [*c]f32, b: [*c]cf32, ldb: isize, v: [*c]cf32, ldv: isize, s: [*c]cf32, lds: isize, zwork: [*c]cf32, lzwork: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zgedmdq_work(layout: c_int, jobs: u8, jobz: u8, jobr: u8, jobq: u8, jobt: u8, jobf: u8, whtsvd: isize, m: isize, n: isize, f: [*c]cf64, ldf: isize, x: [*c]cf64, ldx: isize, y: [*c]cf64, ldy: isize, nrnk: isize, tol: [*c]f64, k: isize, eigs: [*c]cf64, z: [*c]cf64, ldz: isize, res: [*c]f64, b: [*c]cf64, ldb: isize, v: [*c]cf64, ldv: isize, s: [*c]cf64, lds: isize, zwork: [*c]cf64, lzwork: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sgedmdq_work = LAPACKE_sgedmdq_work;
pub const dgedmdq_work = LAPACKE_dgedmdq_work;
pub const cgedmdq_work = LAPACKE_cgedmdq_work;
pub const zgedmdq_work = LAPACKE_zgedmdq_work;

extern fn LAPACKE_sgesv_work(layout: c_int, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgesv_work(layout: c_int, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgesv_work(layout: c_int, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgesv_work(layout: c_int, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
extern fn LAPACKE_dsgesv_work(layout: c_int, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, work: [*c]f64, swork: [*c]f32, iter: [*c]isize) isize;
extern fn LAPACKE_zcgesv_work(layout: c_int, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, work: [*c]cf64, swork: [*c]cf32, rwork: [*c]f64, iter: [*c]isize) isize;
pub const sgesv_work = LAPACKE_sgesv_work;
pub const dgesv_work = LAPACKE_dgesv_work;
pub const cgesv_work = LAPACKE_cgesv_work;
pub const zgesv_work = LAPACKE_zgesv_work;
pub const dsgesv_work = LAPACKE_dsgesv_work;
pub const zcgesv_work = LAPACKE_zcgesv_work;

extern fn LAPACKE_sgesvd_work(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]f32, lda: isize, s: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgesvd_work(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]f64, lda: isize, s: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgesvd_work(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, s: [*c]f32, u: [*c]cf32, ldu: isize, vt: [*c]cf32, ldvt: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zgesvd_work(layout: c_int, jobu: u8, jobvt: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, s: [*c]f64, u: [*c]cf64, ldu: isize, vt: [*c]cf64, ldvt: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sgesvd_work = LAPACKE_sgesvd_work;
pub const dgesvd_work = LAPACKE_dgesvd_work;
pub const cgesvd_work = LAPACKE_cgesvd_work;
pub const zgesvd_work = LAPACKE_zgesvd_work;

extern fn LAPACKE_sgesvdx_work(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, ns: [*c]isize, s: [*c]f32, u: [*c]f32, ldu: isize, vt: [*c]f32, ldvt: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dgesvdx_work(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, ns: [*c]isize, s: [*c]f64, u: [*c]f64, ldu: isize, vt: [*c]f64, ldvt: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cgesvdx_work(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, ns: [*c]isize, s: [*c]f32, u: [*c]cf32, ldu: isize, vt: [*c]cf32, ldvt: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_zgesvdx_work(layout: c_int, jobu: u8, jobvt: u8, range: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, ns: [*c]isize, s: [*c]f64, u: [*c]cf64, ldu: isize, vt: [*c]cf64, ldvt: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize) isize;
pub const sgesvdx_work = LAPACKE_sgesvdx_work;
pub const dgesvdx_work = LAPACKE_dgesvdx_work;
pub const cgesvdx_work = LAPACKE_cgesvdx_work;
pub const zgesvdx_work = LAPACKE_zgesvdx_work;

extern fn LAPACKE_sgesvdq_work(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f32, lda: isize, s: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, numrank: [*c]isize, iwork: [*c]isize, liwork: isize, work: [*c]f32, lwork: isize, rwork: [*c]f32, lrwork: isize) isize;
extern fn LAPACKE_dgesvdq_work(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f64, lda: isize, s: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, numrank: [*c]isize, iwork: [*c]isize, liwork: isize, work: [*c]f64, lwork: isize, rwork: [*c]f64, lrwork: isize) isize;
extern fn LAPACKE_cgesvdq_work(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, s: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, numrank: [*c]isize, iwork: [*c]isize, liwork: isize, cwork: [*c]cf32, lcwork: isize, rwork: [*c]f32, lrwork: isize) isize;
extern fn LAPACKE_zgesvdq_work(layout: c_int, joba: u8, jobp: u8, jobr: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, s: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, numrank: [*c]isize, iwork: [*c]isize, liwork: isize, cwork: [*c]cf64, lcwork: isize, rwork: [*c]f64, lrwork: isize) isize;
pub const sgesvdq_work = LAPACKE_sgesvdq_work;
pub const dgesvdq_work = LAPACKE_dgesvdq_work;
pub const cgesvdq_work = LAPACKE_cgesvdq_work;
pub const zgesvdq_work = LAPACKE_zgesvdq_work;

extern fn LAPACKE_sgesvj_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f32, lda: isize, sva: [*c]f32, mv: isize, v: [*c]f32, ldv: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgesvj_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]f64, lda: isize, sva: [*c]f64, mv: isize, v: [*c]f64, ldv: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgesvj_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf32, lda: isize, sva: [*c]f32, mv: isize, v: [*c]cf32, ldv: isize, cwork: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize) isize;
extern fn LAPACKE_zgesvj_work(layout: c_int, joba: u8, jobu: u8, jobv: u8, m: isize, n: isize, a: [*c]cf64, lda: isize, sva: [*c]f64, mv: isize, v: [*c]cf64, ldv: isize, cwork: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize) isize;
pub const sgesvj_work = LAPACKE_sgesvj_work;
pub const dgesvj_work = LAPACKE_dgesvj_work;
pub const cgesvj_work = LAPACKE_cgesvj_work;
pub const zgesvj_work = LAPACKE_zgesvj_work;

extern fn LAPACKE_sgesvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgesvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgesvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgesvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgesvx_work = LAPACKE_sgesvx_work;
pub const dgesvx_work = LAPACKE_dgesvx_work;
pub const cgesvx_work = LAPACKE_cgesvx_work;
pub const zgesvx_work = LAPACKE_zgesvx_work;

extern fn LAPACKE_sgesvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgesvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgesvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f32, c: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgesvxx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, r: [*c]f64, c: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgesvxx_work = LAPACKE_sgesvxx_work;
pub const dgesvxx_work = LAPACKE_dgesvxx_work;
pub const cgesvxx_work = LAPACKE_cgesvxx_work;
pub const zgesvxx_work = LAPACKE_zgesvxx_work;

extern fn LAPACKE_sgetf2_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgetf2_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgetf2_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgetf2_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const sgetf2_work = LAPACKE_sgetf2_work;
pub const dgetf2_work = LAPACKE_dgetf2_work;
pub const cgetf2_work = LAPACKE_cgetf2_work;
pub const zgetf2_work = LAPACKE_zgetf2_work;

extern fn LAPACKE_sgetrf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgetrf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgetrf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgetrf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const sgetrf_work = LAPACKE_sgetrf_work;
pub const dgetrf_work = LAPACKE_dgetrf_work;
pub const cgetrf_work = LAPACKE_cgetrf_work;
pub const zgetrf_work = LAPACKE_zgetrf_work;

extern fn LAPACKE_sgetrf2_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgetrf2_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgetrf2_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgetrf2_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const sgetrf2_work = LAPACKE_sgetrf2_work;
pub const dgetrf2_work = LAPACKE_dgetrf2_work;
pub const cgetrf2_work = LAPACKE_cgetrf2_work;
pub const zgetrf2_work = LAPACKE_zgetrf2_work;

extern fn LAPACKE_sgetri_work(layout: c_int, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgetri_work(layout: c_int, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgetri_work(layout: c_int, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgetri_work(layout: c_int, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64, lwork: isize) isize;
pub const sgetri_work = LAPACKE_sgetri_work;
pub const dgetri_work = LAPACKE_dgetri_work;
pub const cgetri_work = LAPACKE_cgetri_work;
pub const zgetri_work = LAPACKE_zgetri_work;

extern fn LAPACKE_sgetrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgetrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgetrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgetrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const sgetrs_work = LAPACKE_sgetrs_work;
pub const dgetrs_work = LAPACKE_dgetrs_work;
pub const cgetrs_work = LAPACKE_cgetrs_work;
pub const zgetrs_work = LAPACKE_zgetrs_work;

extern fn LAPACKE_sggbak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f32, rscale: [*c]const f32, m: isize, v: [*c]f32, ldv: isize) isize;
extern fn LAPACKE_dggbak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f64, rscale: [*c]const f64, m: isize, v: [*c]f64, ldv: isize) isize;
extern fn LAPACKE_cggbak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f32, rscale: [*c]const f32, m: isize, v: [*c]cf32, ldv: isize) isize;
extern fn LAPACKE_zggbak_work(layout: c_int, job: u8, side: u8, n: isize, ilo: isize, ihi: isize, lscale: [*c]const f64, rscale: [*c]const f64, m: isize, v: [*c]cf64, ldv: isize) isize;
pub const sggbak_work = LAPACKE_sggbak_work;
pub const dggbak_work = LAPACKE_dggbak_work;
pub const cggbak_work = LAPACKE_cggbak_work;
pub const zggbak_work = LAPACKE_zggbak_work;

extern fn LAPACKE_sggbal_work(layout: c_int, job: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dggbal_work(layout: c_int, job: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cggbal_work(layout: c_int, job: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_zggbal_work(layout: c_int, job: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64, work: [*c]f64) isize;
pub const sggbal_work = LAPACKE_sggbal_work;
pub const dggbal_work = LAPACKE_dggbal_work;
pub const cggbal_work = LAPACKE_cggbal_work;
pub const zggbal_work = LAPACKE_zggbal_work;

extern fn LAPACKE_sgges_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f32, [*c]const f32, [*c]const f32) isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, sdim: [*c]isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vsl: [*c]f32, ldvsl: isize, vsr: [*c]f32, ldvsr: isize, work: [*c]f32, lwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_dgges_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f64, [*c]const f64, [*c]const f64) isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, sdim: [*c]isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vsl: [*c]f64, ldvsl: isize, vsr: [*c]f64, ldvsr: isize, work: [*c]f64, lwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_cgges_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf32, [*c]const cf32) isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, sdim: [*c]isize, alpha: [*c]cf32, beta: [*c]cf32, vsl: [*c]cf32, ldvsl: isize, vsr: [*c]cf32, ldvsr: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, bwork: [*c]isize) isize;
extern fn LAPACKE_zgges_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf64, [*c]const cf64) isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, sdim: [*c]isize, alpha: [*c]cf64, beta: [*c]cf64, vsl: [*c]cf64, ldvsl: isize, vsr: [*c]cf64, ldvsr: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, bwork: [*c]isize) isize;
pub const sgges_work = LAPACKE_sgges_work;
pub const dgges_work = LAPACKE_dgges_work;
pub const cgges_work = LAPACKE_cgges_work;
pub const zgges_work = LAPACKE_zgges_work;

extern fn LAPACKE_sgges3_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f32, [*c]const f32, [*c]const f32) isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, sdim: [*c]isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vsl: [*c]f32, ldvsl: isize, vsr: [*c]f32, ldvsr: isize, work: [*c]f32, lwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_dgges3_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f64, [*c]const f64, [*c]const f64) isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, sdim: [*c]isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vsl: [*c]f64, ldvsl: isize, vsr: [*c]f64, ldvsr: isize, work: [*c]f64, lwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_cgges3_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf32, [*c]const cf32) isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, sdim: [*c]isize, alpha: [*c]cf32, beta: [*c]cf32, vsl: [*c]cf32, ldvsl: isize, vsr: [*c]cf32, ldvsr: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, bwork: [*c]isize) isize;
extern fn LAPACKE_zgges3_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf64, [*c]const cf64) isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, sdim: [*c]isize, alpha: [*c]cf64, beta: [*c]cf64, vsl: [*c]cf64, ldvsl: isize, vsr: [*c]cf64, ldvsr: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, bwork: [*c]isize) isize;
pub const sgges3_work = LAPACKE_sgges3_work;
pub const dgges3_work = LAPACKE_dgges3_work;
pub const cgges3_work = LAPACKE_cgges3_work;
pub const zgges3_work = LAPACKE_zgges3_work;

extern fn LAPACKE_sggesx_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f32, [*c]const f32, [*c]const f32) isize, sense: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, sdim: [*c]isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vsl: [*c]f32, ldvsl: isize, vsr: [*c]f32, ldvsr: isize, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_dggesx_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const f64, [*c]const f64, [*c]const f64) isize, sense: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, sdim: [*c]isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vsl: [*c]f64, ldvsl: isize, vsr: [*c]f64, ldvsr: isize, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_cggesx_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf32, [*c]const cf32) isize, sense: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, sdim: [*c]isize, alpha: [*c]cf32, beta: [*c]cf32, vsl: [*c]cf32, ldvsl: isize, vsr: [*c]cf32, ldvsr: isize, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize, liwork: isize, bwork: [*c]isize) isize;
extern fn LAPACKE_zggesx_work(layout: c_int, jobvsl: u8, jobvsr: u8, sort: u8, selctg: *const fn ([*c]const cf64, [*c]const cf64) isize, sense: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, sdim: [*c]isize, alpha: [*c]cf64, beta: [*c]cf64, vsl: [*c]cf64, ldvsl: isize, vsr: [*c]cf64, ldvsr: isize, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize, liwork: isize, bwork: [*c]isize) isize;
pub const sggesx_work = LAPACKE_sggesx_work;
pub const dggesx_work = LAPACKE_dggesx_work;
pub const cggesx_work = LAPACKE_cggesx_work;
pub const zggesx_work = LAPACKE_zggesx_work;

extern fn LAPACKE_sggev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dggev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cggev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zggev_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sggev_work = LAPACKE_sggev_work;
pub const dggev_work = LAPACKE_dggev_work;
pub const cggev_work = LAPACKE_cggev_work;
pub const zggev_work = LAPACKE_zggev_work;

extern fn LAPACKE_sggev3_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dggev3_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cggev3_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zggev3_work(layout: c_int, jobvl: u8, jobvr: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const sggev3_work = LAPACKE_sggev3_work;
pub const dggev3_work = LAPACKE_dggev3_work;
pub const cggev3_work = LAPACKE_cggev3_work;
pub const zggev3_work = LAPACKE_zggev3_work;

extern fn LAPACKE_sggevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32, abnrm: [*c]f32, bbnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, bwork: [*c]isize) isize;
extern fn LAPACKE_dggevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64, abnrm: [*c]f64, bbnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, bwork: [*c]isize) isize;
extern fn LAPACKE_cggevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f32, rscale: [*c]f32, abnrm: [*c]f32, bbnrm: [*c]f32, rconde: [*c]f32, rcondv: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize, bwork: [*c]isize) isize;
extern fn LAPACKE_zggevx_work(layout: c_int, balanc: u8, jobvl: u8, jobvr: u8, sense: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, ilo: [*c]isize, ihi: [*c]isize, lscale: [*c]f64, rscale: [*c]f64, abnrm: [*c]f64, bbnrm: [*c]f64, rconde: [*c]f64, rcondv: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize, bwork: [*c]isize) isize;
pub const sggevx_work = LAPACKE_sggevx_work;
pub const dggevx_work = LAPACKE_dggevx_work;
pub const cggevx_work = LAPACKE_cggevx_work;
pub const zggevx_work = LAPACKE_zggevx_work;

extern fn LAPACKE_sggglm_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, d: [*c]f32, x: [*c]f32, y: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dggglm_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, d: [*c]f64, x: [*c]f64, y: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cggglm_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, d: [*c]cf32, x: [*c]cf32, y: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zggglm_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, d: [*c]cf64, x: [*c]cf64, y: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sggglm_work = LAPACKE_sggglm_work;
pub const dggglm_work = LAPACKE_dggglm_work;
pub const cggglm_work = LAPACKE_cggglm_work;
pub const zggglm_work = LAPACKE_zggglm_work;

extern fn LAPACKE_sgghrd_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dgghrd_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize) isize;
extern fn LAPACKE_cgghrd_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zgghrd_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize) isize;
pub const sgghrd_work = LAPACKE_sgghrd_work;
pub const dgghrd_work = LAPACKE_dgghrd_work;
pub const cgghrd_work = LAPACKE_cgghrd_work;
pub const zgghrd_work = LAPACKE_zgghrd_work;

extern fn LAPACKE_sgghd3_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgghd3_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgghd3_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgghd3_work(layout: c_int, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgghd3_work = LAPACKE_sgghd3_work;
pub const dgghd3_work = LAPACKE_dgghd3_work;
pub const cgghd3_work = LAPACKE_cgghd3_work;
pub const zgghd3_work = LAPACKE_zgghd3_work;

extern fn LAPACKE_sgglse_work(layout: c_int, m: isize, n: isize, p: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, c: [*c]f32, d: [*c]f32, x: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgglse_work(layout: c_int, m: isize, n: isize, p: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, c: [*c]f64, d: [*c]f64, x: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgglse_work(layout: c_int, m: isize, n: isize, p: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, c: [*c]cf32, d: [*c]cf32, x: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgglse_work(layout: c_int, m: isize, n: isize, p: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, c: [*c]cf64, d: [*c]cf64, x: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sgglse_work = LAPACKE_sgglse_work;
pub const dgglse_work = LAPACKE_dgglse_work;
pub const cgglse_work = LAPACKE_cgglse_work;
pub const zgglse_work = LAPACKE_zgglse_work;

extern fn LAPACKE_sggqrf_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f32, lda: isize, taua: [*c]f32, b: [*c]f32, ldb: isize, taub: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dggqrf_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]f64, lda: isize, taua: [*c]f64, b: [*c]f64, ldb: isize, taub: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cggqrf_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf32, lda: isize, taua: [*c]cf32, b: [*c]cf32, ldb: isize, taub: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zggqrf_work(layout: c_int, n: isize, m: isize, p: isize, a: [*c]cf64, lda: isize, taua: [*c]cf64, b: [*c]cf64, ldb: isize, taub: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sggqrf_work = LAPACKE_sggqrf_work;
pub const dggqrf_work = LAPACKE_dggqrf_work;
pub const cggqrf_work = LAPACKE_cggqrf_work;
pub const zggqrf_work = LAPACKE_zggqrf_work;

extern fn LAPACKE_sggrqf_work(layout: c_int, m: isize, p: isize, n: isize, a: [*c]f32, lda: isize, taua: [*c]f32, b: [*c]f32, ldb: isize, taub: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dggrqf_work(layout: c_int, m: isize, p: isize, n: isize, a: [*c]f64, lda: isize, taua: [*c]f64, b: [*c]f64, ldb: isize, taub: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cggrqf_work(layout: c_int, m: isize, p: isize, n: isize, a: [*c]cf32, lda: isize, taua: [*c]cf32, b: [*c]cf32, ldb: isize, taub: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zggrqf_work(layout: c_int, m: isize, p: isize, n: isize, a: [*c]cf64, lda: isize, taua: [*c]cf64, b: [*c]cf64, ldb: isize, taub: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sggrqf_work = LAPACKE_sggrqf_work;
pub const dggrqf_work = LAPACKE_dggrqf_work;
pub const cggrqf_work = LAPACKE_cggrqf_work;
pub const zggrqf_work = LAPACKE_zggrqf_work;

extern fn LAPACKE_sggsvd_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dggsvd_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cggsvd_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, work: [*c]cf32, rwork: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_zggsvd_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, work: [*c]cf64, rwork: [*c]f64, iwork: [*c]isize) isize;
pub const sggsvd_work = LAPACKE_sggsvd_work;
pub const dggsvd_work = LAPACKE_dggsvd_work;
pub const cggsvd_work = LAPACKE_cggsvd_work;
pub const zggsvd_work = LAPACKE_zggsvd_work;

extern fn LAPACKE_sggsvd3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dggsvd3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_cggsvd3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]f32, beta: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_zggsvd3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, n: isize, p: isize, k: [*c]isize, l: [*c]isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]f64, beta: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize) isize;
pub const sggsvd3_work = LAPACKE_sggsvd3_work;
pub const dggsvd3_work = LAPACKE_dggsvd3_work;
pub const cggsvd3_work = LAPACKE_cggsvd3_work;
pub const zggsvd3_work = LAPACKE_zggsvd3_work;

extern fn LAPACKE_sggsvp_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, iwork: [*c]isize, tau: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dggsvp_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, iwork: [*c]isize, tau: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cggsvp_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, iwork: [*c]isize, rwork: [*c]f32, tau: [*c]cf32, work: [*c]cf32) isize;
extern fn LAPACKE_zggsvp_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, iwork: [*c]isize, rwork: [*c]f64, tau: [*c]cf64, work: [*c]cf64) isize;
pub const sggsvp_work = LAPACKE_sggsvp_work;
pub const dggsvp_work = LAPACKE_dggsvp_work;
pub const cggsvp_work = LAPACKE_cggsvp_work;
pub const zggsvp_work = LAPACKE_zggsvp_work;

extern fn LAPACKE_sggsvp3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, iwork: [*c]isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dggsvp3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, iwork: [*c]isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cggsvp3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, tola: f32, tolb: f32, k: [*c]isize, l: [*c]isize, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, iwork: [*c]isize, rwork: [*c]f32, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zggsvp3_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, tola: f64, tolb: f64, k: [*c]isize, l: [*c]isize, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, iwork: [*c]isize, rwork: [*c]f64, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const sggsvp3_work = LAPACKE_sggsvp3_work;
pub const dggsvp3_work = LAPACKE_dggsvp3_work;
pub const cggsvp3_work = LAPACKE_cggsvp3_work;
pub const zggsvp3_work = LAPACKE_zggsvp3_work;

extern fn LAPACKE_sgtcon_work(norm: u8, n: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, du2: [*c]const f32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgtcon_work(norm: u8, n: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, du2: [*c]const f64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgtcon_work(norm: u8, n: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, du2: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zgtcon_work(norm: u8, n: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, du2: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const sgtcon_work = LAPACKE_sgtcon_work;
pub const dgtcon_work = LAPACKE_dgtcon_work;
pub const cgtcon_work = LAPACKE_cgtcon_work;
pub const zgtcon_work = LAPACKE_zgtcon_work;

extern fn LAPACKE_sgtrfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, dlf: [*c]const f32, df: [*c]const f32, duf: [*c]const f32, du2: [*c]const f32, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgtrfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, dlf: [*c]const f64, df: [*c]const f64, duf: [*c]const f64, du2: [*c]const f64, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgtrfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, dlf: [*c]const cf32, df: [*c]const cf32, duf: [*c]const cf32, du2: [*c]const cf32, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgtrfs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, dlf: [*c]const cf64, df: [*c]const cf64, duf: [*c]const cf64, du2: [*c]const cf64, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgtrfs_work = LAPACKE_sgtrfs_work;
pub const dgtrfs_work = LAPACKE_dgtrfs_work;
pub const cgtrfs_work = LAPACKE_cgtrfs_work;
pub const zgtrfs_work = LAPACKE_zgtrfs_work;

extern fn LAPACKE_sgtsv_work(layout: c_int, n: isize, nrhs: isize, dl: [*c]f32, d: [*c]f32, du: [*c]f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgtsv_work(layout: c_int, n: isize, nrhs: isize, dl: [*c]f64, d: [*c]f64, du: [*c]f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgtsv_work(layout: c_int, n: isize, nrhs: isize, dl: [*c]cf32, d: [*c]cf32, du: [*c]cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgtsv_work(layout: c_int, n: isize, nrhs: isize, dl: [*c]cf64, d: [*c]cf64, du: [*c]cf64, b: [*c]cf64, ldb: isize) isize;
pub const sgtsv_work = LAPACKE_sgtsv_work;
pub const dgtsv_work = LAPACKE_dgtsv_work;
pub const cgtsv_work = LAPACKE_cgtsv_work;
pub const zgtsv_work = LAPACKE_zgtsv_work;

extern fn LAPACKE_sgtsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, dlf: [*c]f32, df: [*c]f32, duf: [*c]f32, du2: [*c]f32, ipiv: [*c]isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dgtsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, dlf: [*c]f64, df: [*c]f64, duf: [*c]f64, du2: [*c]f64, ipiv: [*c]isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cgtsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, dlf: [*c]cf32, df: [*c]cf32, duf: [*c]cf32, du2: [*c]cf32, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zgtsvx_work(layout: c_int, fact: u8, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, dlf: [*c]cf64, df: [*c]cf64, duf: [*c]cf64, du2: [*c]cf64, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sgtsvx_work = LAPACKE_sgtsvx_work;
pub const dgtsvx_work = LAPACKE_dgtsvx_work;
pub const cgtsvx_work = LAPACKE_cgtsvx_work;
pub const zgtsvx_work = LAPACKE_zgtsvx_work;

extern fn LAPACKE_sgttrf_work(n: isize, dl: [*c]f32, d: [*c]f32, du: [*c]f32, du2: [*c]f32, ipiv: [*c]isize) isize;
extern fn LAPACKE_dgttrf_work(n: isize, dl: [*c]f64, d: [*c]f64, du: [*c]f64, du2: [*c]f64, ipiv: [*c]isize) isize;
extern fn LAPACKE_cgttrf_work(n: isize, dl: [*c]cf32, d: [*c]cf32, du: [*c]cf32, du2: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zgttrf_work(n: isize, dl: [*c]cf64, d: [*c]cf64, du: [*c]cf64, du2: [*c]cf64, ipiv: [*c]isize) isize;
pub const sgttrf_work = LAPACKE_sgttrf_work;
pub const dgttrf_work = LAPACKE_dgttrf_work;
pub const cgttrf_work = LAPACKE_cgttrf_work;
pub const zgttrf_work = LAPACKE_zgttrf_work;

extern fn LAPACKE_sgttrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f32, d: [*c]const f32, du: [*c]const f32, du2: [*c]const f32, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgttrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const f64, d: [*c]const f64, du: [*c]const f64, du2: [*c]const f64, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgttrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf32, d: [*c]const cf32, du: [*c]const cf32, du2: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgttrs_work(layout: c_int, trans: u8, n: isize, nrhs: isize, dl: [*c]const cf64, d: [*c]const cf64, du: [*c]const cf64, du2: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const sgttrs_work = LAPACKE_sgttrs_work;
pub const dgttrs_work = LAPACKE_dgttrs_work;
pub const cgttrs_work = LAPACKE_cgttrs_work;
pub const zgttrs_work = LAPACKE_zgttrs_work;

extern fn LAPACKE_chbev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhbev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chbev_work = LAPACKE_chbev_work;
pub const zhbev_work = LAPACKE_zhbev_work;

extern fn LAPACKE_chbevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zhbevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const chbevd_work = LAPACKE_chbevd_work;
pub const zhbevd_work = LAPACKE_zhbevd_work;

extern fn LAPACKE_chbevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, q: [*c]cf32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhbevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, q: [*c]cf64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const chbevx_work = LAPACKE_chbevx_work;
pub const zhbevx_work = LAPACKE_zhbevx_work;

extern fn LAPACKE_chbgst_work(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]const cf32, ldbb: isize, x: [*c]cf32, ldx: isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhbgst_work(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]const cf64, ldbb: isize, x: [*c]cf64, ldx: isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chbgst_work = LAPACKE_chbgst_work;
pub const zhbgst_work = LAPACKE_zhbgst_work;

extern fn LAPACKE_chbgv_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]cf32, ldbb: isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhbgv_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]cf64, ldbb: isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chbgv_work = LAPACKE_chbgv_work;
pub const zhbgv_work = LAPACKE_zhbgv_work;

extern fn LAPACKE_chbgvd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]cf32, ldbb: isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zhbgvd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]cf64, ldbb: isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const chbgvd_work = LAPACKE_chbgvd_work;
pub const zhbgvd_work = LAPACKE_zhbgvd_work;

extern fn LAPACKE_chbgvx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf32, ldab: isize, bb: [*c]cf32, ldbb: isize, q: [*c]cf32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhbgvx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]cf64, ldab: isize, bb: [*c]cf64, ldbb: isize, q: [*c]cf64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const chbgvx_work = LAPACKE_chbgvx_work;
pub const zhbgvx_work = LAPACKE_zhbgvx_work;

extern fn LAPACKE_chbtrd_work(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]cf32, ldq: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zhbtrd_work(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]cf64, ldq: isize, work: [*c]cf64) isize;
pub const chbtrd_work = LAPACKE_chbtrd_work;
pub const zhbtrd_work = LAPACKE_zhbtrd_work;

extern fn LAPACKE_checon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zhecon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const checon_work = LAPACKE_checon_work;
pub const zhecon_work = LAPACKE_zhecon_work;

extern fn LAPACKE_cheequb_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zheequb_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64, work: [*c]cf64) isize;
pub const cheequb_work = LAPACKE_cheequb_work;
pub const zheequb_work = LAPACKE_zheequb_work;

extern fn LAPACKE_cheev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zheev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const cheev_work = LAPACKE_cheev_work;
pub const zheev_work = LAPACKE_zheev_work;

extern fn LAPACKE_cheevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zheevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const cheevd_work = LAPACKE_cheevd_work;
pub const zheevd_work = LAPACKE_zheevd_work;

extern fn LAPACKE_cheevr_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, isuppz: [*c]isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zheevr_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, isuppz: [*c]isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const cheevr_work = LAPACKE_cheevr_work;
pub const zheevr_work = LAPACKE_zheevr_work;

extern fn LAPACKE_cheevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zheevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const cheevx_work = LAPACKE_cheevx_work;
pub const zheevx_work = LAPACKE_zheevx_work;

extern fn LAPACKE_chegst_work(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]const cf32, ldb: isize) isize;
extern fn LAPACKE_zhegst_work(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]const cf64, ldb: isize) isize;
pub const chegst_work = LAPACKE_chegst_work;
pub const zhegst_work = LAPACKE_zhegst_work;

extern fn LAPACKE_chegv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zhegv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const chegv_work = LAPACKE_chegv_work;
pub const zhegv_work = LAPACKE_zhegv_work;

extern fn LAPACKE_chegvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zhegvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const chegvd_work = LAPACKE_chegvd_work;
pub const zhegvd_work = LAPACKE_zhegvd_work;

extern fn LAPACKE_chegvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhegvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const chegvx_work = LAPACKE_chegvx_work;
pub const zhegvx_work = LAPACKE_zhegvx_work;

extern fn LAPACKE_cherfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zherfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const cherfs_work = LAPACKE_cherfs_work;
pub const zherfs_work = LAPACKE_zherfs_work;

extern fn LAPACKE_cherfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zherfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const cherfsx_work = LAPACKE_cherfsx_work;
pub const zherfsx_work = LAPACKE_zherfsx_work;

extern fn LAPACKE_chesv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhesv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const chesv_work = LAPACKE_chesv_work;
pub const zhesv_work = LAPACKE_zhesv_work;

extern fn LAPACKE_chesvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zhesvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const chesvx_work = LAPACKE_chesvx_work;
pub const zhesvx_work = LAPACKE_zhesvx_work;

extern fn LAPACKE_chesvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhesvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chesvxx_work = LAPACKE_chesvxx_work;
pub const zhesvxx_work = LAPACKE_zhesvxx_work;

extern fn LAPACKE_chetrd_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, d: [*c]f32, e: [*c]f32, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhetrd_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, d: [*c]f64, e: [*c]f64, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const chetrd_work = LAPACKE_chetrd_work;
pub const zhetrd_work = LAPACKE_zhetrd_work;

extern fn LAPACKE_chetrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhetrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const chetrf_work = LAPACKE_chetrf_work;
pub const zhetrf_work = LAPACKE_zhetrf_work;

extern fn LAPACKE_chetri_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32) isize;
extern fn LAPACKE_zhetri_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64) isize;
pub const chetri_work = LAPACKE_chetri_work;
pub const zhetri_work = LAPACKE_zhetri_work;

extern fn LAPACKE_chetrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhetrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const chetrs_work = LAPACKE_chetrs_work;
pub const zhetrs_work = LAPACKE_zhetrs_work;

extern fn LAPACKE_chfrk_work(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f32, a: [*c]const cf32, lda: isize, beta: f32, c: [*c]cf32) isize;
extern fn LAPACKE_zhfrk_work(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f64, a: [*c]const cf64, lda: isize, beta: f64, c: [*c]cf64) isize;
pub const chfrk_work = LAPACKE_chfrk_work;
pub const zhfrk_work = LAPACKE_zhfrk_work;

extern fn LAPACKE_shgeqz_work(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f32, ldh: isize, t: [*c]f32, ldt: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dhgeqz_work(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f64, ldh: isize, t: [*c]f64, ldt: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_chgeqz_work(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf32, ldh: isize, t: [*c]cf32, ldt: isize, alpha: [*c]cf32, beta: [*c]cf32, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zhgeqz_work(layout: c_int, job: u8, compq: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf64, ldh: isize, t: [*c]cf64, ldt: isize, alpha: [*c]cf64, beta: [*c]cf64, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const shgeqz_work = LAPACKE_shgeqz_work;
pub const dhgeqz_work = LAPACKE_dhgeqz_work;
pub const chgeqz_work = LAPACKE_chgeqz_work;
pub const zhgeqz_work = LAPACKE_zhgeqz_work;

extern fn LAPACKE_chpcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zhpcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const chpcon_work = LAPACKE_chpcon_work;
pub const zhpcon_work = LAPACKE_zhpcon_work;

extern fn LAPACKE_chpev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhpev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chpev_work = LAPACKE_chpev_work;
pub const zhpev_work = LAPACKE_zhpev_work;

extern fn LAPACKE_chpevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zhpevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const chpevd_work = LAPACKE_chpevd_work;
pub const zhpevd_work = LAPACKE_zhpevd_work;

extern fn LAPACKE_chpevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhpevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const chpevx_work = LAPACKE_chpevx_work;
pub const zhpevx_work = LAPACKE_zhpevx_work;

extern fn LAPACKE_chpgst_work(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]const cf32) isize;
extern fn LAPACKE_zhpgst_work(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]const cf64) isize;
pub const chpgst_work = LAPACKE_chpgst_work;
pub const zhpgst_work = LAPACKE_zhpgst_work;

extern fn LAPACKE_chpgv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhpgv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chpgv_work = LAPACKE_chpgv_work;
pub const zhpgv_work = LAPACKE_zhpgv_work;

extern fn LAPACKE_chpgvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]cf32, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zhpgvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]cf64, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const chpgvd_work = LAPACKE_chpgvd_work;
pub const zhpgvd_work = LAPACKE_zhpgvd_work;

extern fn LAPACKE_chpgvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf32, bp: [*c]cf32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhpgvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]cf64, bp: [*c]cf64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const chpgvx_work = LAPACKE_chpgvx_work;
pub const zhpgvx_work = LAPACKE_zhpgvx_work;

extern fn LAPACKE_chprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]const cf32, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]const cf64, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chprfs_work = LAPACKE_chprfs_work;
pub const zhprfs_work = LAPACKE_zhprfs_work;

extern fn LAPACKE_chpsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhpsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const chpsv_work = LAPACKE_chpsv_work;
pub const zhpsv_work = LAPACKE_zhpsv_work;

extern fn LAPACKE_chpsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]cf32, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zhpsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]cf64, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const chpsvx_work = LAPACKE_chpsvx_work;
pub const zhpsvx_work = LAPACKE_zhpsvx_work;

extern fn LAPACKE_chptrd_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, d: [*c]f32, e: [*c]f32, tau: [*c]cf32) isize;
extern fn LAPACKE_zhptrd_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, d: [*c]f64, e: [*c]f64, tau: [*c]cf64) isize;
pub const chptrd_work = LAPACKE_chptrd_work;
pub const zhptrd_work = LAPACKE_zhptrd_work;

extern fn LAPACKE_chptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zhptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]isize) isize;
pub const chptrf_work = LAPACKE_chptrf_work;
pub const zhptrf_work = LAPACKE_zhptrf_work;

extern fn LAPACKE_chptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]const isize, work: [*c]cf32) isize;
extern fn LAPACKE_zhptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]const isize, work: [*c]cf64) isize;
pub const chptri_work = LAPACKE_chptri_work;
pub const zhptri_work = LAPACKE_zhptri_work;

extern fn LAPACKE_chptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const chptrs_work = LAPACKE_chptrs_work;
pub const zhptrs_work = LAPACKE_zhptrs_work;

extern fn LAPACKE_shsein_work(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]isize, n: isize, h: [*c]const f32, ldh: isize, wr: [*c]f32, wi: [*c]const f32, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]f32, ifaill: [*c]isize, ifailr: [*c]isize) isize;
extern fn LAPACKE_dhsein_work(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]isize, n: isize, h: [*c]const f64, ldh: isize, wr: [*c]f64, wi: [*c]const f64, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]f64, ifaill: [*c]isize, ifailr: [*c]isize) isize;
extern fn LAPACKE_chsein_work(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]const isize, n: isize, h: [*c]const cf32, ldh: isize, w: [*c]cf32, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]cf32, rwork: [*c]f32, ifaill: [*c]isize, ifailr: [*c]isize) isize;
extern fn LAPACKE_zhsein_work(layout: c_int, job: u8, eigsrc: u8, initv: u8, select: [*c]const isize, n: isize, h: [*c]const cf64, ldh: isize, w: [*c]cf64, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]cf64, rwork: [*c]f64, ifaill: [*c]isize, ifailr: [*c]isize) isize;
pub const shsein_work = LAPACKE_shsein_work;
pub const dhsein_work = LAPACKE_dhsein_work;
pub const chsein_work = LAPACKE_chsein_work;
pub const zhsein_work = LAPACKE_zhsein_work;

extern fn LAPACKE_shseqr_work(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f32, ldh: isize, wr: [*c]f32, wi: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dhseqr_work(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]f64, ldh: isize, wr: [*c]f64, wi: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_chseqr_work(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf32, ldh: isize, w: [*c]cf32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhseqr_work(layout: c_int, job: u8, compz: u8, n: isize, ilo: isize, ihi: isize, h: [*c]cf64, ldh: isize, w: [*c]cf64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize) isize;
pub const shseqr_work = LAPACKE_shseqr_work;
pub const dhseqr_work = LAPACKE_dhseqr_work;
pub const chseqr_work = LAPACKE_chseqr_work;
pub const zhseqr_work = LAPACKE_zhseqr_work;

extern fn LAPACKE_clacgv_work(n: isize, x: [*c]cf32, incx: isize) isize;
extern fn LAPACKE_zlacgv_work(n: isize, x: [*c]cf64, incx: isize) isize;
pub const clacgv_work = LAPACKE_clacgv_work;
pub const zlacgv_work = LAPACKE_zlacgv_work;

extern fn LAPACKE_slacn2_work(n: isize, v: [*c]f32, x: [*c]f32, isgn: [*c]isize, est: [*c]f32, kase: [*c]isize, isave: [*c]isize) isize;
extern fn LAPACKE_dlacn2_work(n: isize, v: [*c]f64, x: [*c]f64, isgn: [*c]isize, est: [*c]f64, kase: [*c]isize, isave: [*c]isize) isize;
extern fn LAPACKE_clacn2_work(n: isize, v: [*c]cf32, x: [*c]cf32, est: [*c]f32, kase: [*c]isize, isave: [*c]isize) isize;
extern fn LAPACKE_zlacn2_work(n: isize, v: [*c]cf64, x: [*c]cf64, est: [*c]f64, kase: [*c]isize, isave: [*c]isize) isize;
pub const slacn2_work = LAPACKE_slacn2_work;
pub const dlacn2_work = LAPACKE_dlacn2_work;
pub const clacn2_work = LAPACKE_clacn2_work;
pub const zlacn2_work = LAPACKE_zlacn2_work;

extern fn LAPACKE_slacpy_work(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dlacpy_work(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_clacpy_work(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zlacpy_work(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const slacpy_work = LAPACKE_slacpy_work;
pub const dlacpy_work = LAPACKE_dlacpy_work;
pub const clacpy_work = LAPACKE_clacpy_work;
pub const zlacpy_work = LAPACKE_zlacpy_work;

extern fn LAPACKE_clacp2_work(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zlacp2_work(layout: c_int, uplo: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const clacp2_work = LAPACKE_clacp2_work;
pub const zlacp2_work = LAPACKE_zlacp2_work;

extern fn LAPACKE_slag2d_work(layout: c_int, m: isize, n: isize, sa: [*c]const f32, ldsa: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_dlag2s_work(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, sa: [*c]f32, ldsa: isize) isize;
extern fn LAPACKE_clag2z_work(layout: c_int, m: isize, n: isize, sa: [*c]const cf32, ldsa: isize, a: [*c]cf64, lda: isize) isize;
extern fn LAPACKE_zlag2c_work(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, sa: [*c]cf32, ldsa: isize) isize;
pub const slag2d_work = LAPACKE_slag2d_work;
pub const dlag2s_work = LAPACKE_dlag2s_work;
pub const clag2z_work = LAPACKE_clag2z_work;
pub const zlag2c_work = LAPACKE_zlag2c_work;

extern fn LAPACKE_slagge_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f32, a: [*c]f32, lda: isize, iseed: [*c]isize, work: [*c]f32) isize;
extern fn LAPACKE_dlagge_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f64, a: [*c]f64, lda: isize, iseed: [*c]isize, work: [*c]f64) isize;
extern fn LAPACKE_clagge_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f32, a: [*c]cf32, lda: isize, iseed: [*c]isize, work: [*c]cf32) isize;
extern fn LAPACKE_zlagge_work(layout: c_int, m: isize, n: isize, kl: isize, ku: isize, d: [*c]const f64, a: [*c]cf64, lda: isize, iseed: [*c]isize, work: [*c]cf64) isize;
pub const slagge_work = LAPACKE_slagge_work;
pub const dlagge_work = LAPACKE_dlagge_work;
pub const clagge_work = LAPACKE_clagge_work;
pub const zlagge_work = LAPACKE_zlagge_work;

extern fn LAPACKE_claghe_work(layout: c_int, n: isize, k: isize, d: [*c]const f32, a: [*c]cf32, lda: isize, iseed: [*c]isize, work: [*c]cf32) isize;
extern fn LAPACKE_zlaghe_work(layout: c_int, n: isize, k: isize, d: [*c]const f64, a: [*c]cf64, lda: isize, iseed: [*c]isize, work: [*c]cf64) isize;
pub const claghe_work = LAPACKE_claghe_work;
pub const zlaghe_work = LAPACKE_zlaghe_work;

extern fn LAPACKE_slagsy_work(layout: c_int, n: isize, k: isize, d: [*c]const f32, a: [*c]f32, lda: isize, iseed: [*c]isize, work: [*c]f32) isize;
extern fn LAPACKE_dlagsy_work(layout: c_int, n: isize, k: isize, d: [*c]const f64, a: [*c]f64, lda: isize, iseed: [*c]isize, work: [*c]f64) isize;
extern fn LAPACKE_clagsy_work(layout: c_int, n: isize, k: isize, d: [*c]const f32, a: [*c]cf32, lda: isize, iseed: [*c]isize, work: [*c]cf32) isize;
extern fn LAPACKE_zlagsy_work(layout: c_int, n: isize, k: isize, d: [*c]const f64, a: [*c]cf64, lda: isize, iseed: [*c]isize, work: [*c]cf64) isize;
pub const slagsy_work = LAPACKE_slagsy_work;
pub const dlagsy_work = LAPACKE_dlagsy_work;
pub const clagsy_work = LAPACKE_clagsy_work;
pub const zlagsy_work = LAPACKE_zlagsy_work;

extern fn LAPACKE_slapmr_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_dlapmr_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f64, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_clapmr_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_zlapmr_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf64, ldx: isize, k: [*c]isize) isize;
pub const slapmr_work = LAPACKE_slapmr_work;
pub const dlapmr_work = LAPACKE_dlapmr_work;
pub const clapmr_work = LAPACKE_clapmr_work;
pub const zlapmr_work = LAPACKE_zlapmr_work;

extern fn LAPACKE_slapmt_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_dlapmt_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f64, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_clapmt_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_zlapmt_work(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf64, ldx: isize, k: [*c]isize) isize;
pub const slapmt_work = LAPACKE_slapmt_work;
pub const dlapmt_work = LAPACKE_dlapmt_work;
pub const clapmt_work = LAPACKE_clapmt_work;
pub const zlapmt_work = LAPACKE_zlapmt_work;

extern fn LAPACKE_slartgp_work(f: f32, g: f32, cs: [*c]f32, sn: [*c]f32, r: [*c]f32) isize;
extern fn LAPACKE_dlartgp_work(f: f64, g: f64, cs: [*c]f64, sn: [*c]f64, r: [*c]f64) isize;
pub const slartgp_work = LAPACKE_slartgp_work;
pub const dlartgp_work = LAPACKE_dlartgp_work;

extern fn LAPACKE_slartgs_work(x: f32, y: f32, sigma: f32, cs: [*c]f32, sn: [*c]f32) isize;
extern fn LAPACKE_dlartgs_work(x: f64, y: f64, sigma: f64, cs: [*c]f64, sn: [*c]f64) isize;
pub const slartgs_work = LAPACKE_slartgs_work;
pub const dlartgs_work = LAPACKE_dlartgs_work;

extern fn LAPACKE_slapy2_work(x: f32, y: f32) f32;
extern fn LAPACKE_dlapy2_work(x: f64, y: f64) f64;
pub const slapy2_work = LAPACKE_slapy2_work;
pub const dlapy2_work = LAPACKE_dlapy2_work;

extern fn LAPACKE_slapy3_work(x: f32, y: f32, z: f32) f32;
extern fn LAPACKE_dlapy3_work(x: f64, y: f64, z: f64) f64;
pub const slapy3_work = LAPACKE_slapy3_work;
pub const dlapy3_work = LAPACKE_dlapy3_work;

extern fn LAPACKE_slamch_work(cmach: u8) f32;
extern fn LAPACKE_dlamch_work(cmach: u8) f64;
pub const slamch_work = LAPACKE_slamch_work;
pub const dlamch_work = LAPACKE_dlamch_work;

extern fn LAPACKE_slangb_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f32, ldab: isize, work: [*c]f32) f32;
extern fn LAPACKE_dlangb_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const f64, ldab: isize, work: [*c]f64) f64;
extern fn LAPACKE_clangb_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf32, ldab: isize, work: [*c]f32) f32;
extern fn LAPACKE_zlangb_work(layout: c_int, norm: u8, n: isize, kl: isize, ku: isize, ab: [*c]const cf64, ldab: isize, work: [*c]f64) f64;
pub const slangb_work = LAPACKE_slangb_work;
pub const dlangb_work = LAPACKE_dlangb_work;
pub const clangb_work = LAPACKE_clangb_work;
pub const zlangb_work = LAPACKE_zlangb_work;

extern fn LAPACKE_slange_work(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_dlange_work(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, work: [*c]f64) f64;
extern fn LAPACKE_clange_work(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_zlange_work(layout: c_int, norm: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize, work: [*c]f64) f64;
pub const slange_work = LAPACKE_slange_work;
pub const dlange_work = LAPACKE_dlange_work;
pub const clange_work = LAPACKE_clange_work;
pub const zlange_work = LAPACKE_zlange_work;

extern fn LAPACKE_clanhe_work(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_zlanhe_work(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, work: [*c]f64) f64;
pub const clanhe_work = LAPACKE_clanhe_work;
pub const zlanhe_work = LAPACKE_zlanhe_work;

extern fn LAPACKE_clacrm_work(layout: c_int, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]cf32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_zlacrm_work(layout: c_int, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]cf64, ldc: isize, work: [*c]f64) isize;
pub const clacrm_work = LAPACKE_clacrm_work;
pub const zlacrm_work = LAPACKE_zlacrm_work;

extern fn LAPACKE_clarcm_work(layout: c_int, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_zlarcm_work(layout: c_int, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, work: [*c]f64) isize;
pub const clarcm_work = LAPACKE_clarcm_work;
pub const zlarcm_work = LAPACKE_zlarcm_work;

extern fn LAPACKE_slansy_work(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const f32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_dlansy_work(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const f64, lda: isize, work: [*c]f64) f64;
extern fn LAPACKE_clansy_work(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_zlansy_work(layout: c_int, norm: u8, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, work: [*c]f64) f64;
pub const slansy_work = LAPACKE_slansy_work;
pub const dlansy_work = LAPACKE_dlansy_work;
pub const clansy_work = LAPACKE_clansy_work;
pub const zlansy_work = LAPACKE_zlansy_work;

extern fn LAPACKE_slantr_work(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_dlantr_work(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, work: [*c]f64) f64;
extern fn LAPACKE_clantr_work(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize, work: [*c]f32) f32;
extern fn LAPACKE_zlantr_work(layout: c_int, norm: u8, uplo: u8, diag: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize, work: [*c]f64) f64;
pub const slantr_work = LAPACKE_slantr_work;
pub const dlantr_work = LAPACKE_dlantr_work;
pub const clantr_work = LAPACKE_clantr_work;
pub const zlantr_work = LAPACKE_zlantr_work;

extern fn LAPACKE_slarfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, c: [*c]f32, ldc: isize, work: [*c]f32, ldwork: isize) isize;
extern fn LAPACKE_dlarfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, c: [*c]f64, ldc: isize, work: [*c]f64, ldwork: isize) isize;
extern fn LAPACKE_clarfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, c: [*c]cf32, ldc: isize, work: [*c]cf32, ldwork: isize) isize;
extern fn LAPACKE_zlarfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, c: [*c]cf64, ldc: isize, work: [*c]cf64, ldwork: isize) isize;
pub const slarfb_work = LAPACKE_slarfb_work;
pub const dlarfb_work = LAPACKE_dlarfb_work;
pub const clarfb_work = LAPACKE_clarfb_work;
pub const zlarfb_work = LAPACKE_zlarfb_work;

extern fn LAPACKE_slarfg_work(n: isize, alpha: [*c]f32, x: [*c]f32, incx: isize, tau: [*c]f32) isize;
extern fn LAPACKE_dlarfg_work(n: isize, alpha: [*c]f64, x: [*c]f64, incx: isize, tau: [*c]f64) isize;
extern fn LAPACKE_clarfg_work(n: isize, alpha: [*c]cf32, x: [*c]cf32, incx: isize, tau: [*c]cf32) isize;
extern fn LAPACKE_zlarfg_work(n: isize, alpha: [*c]cf64, x: [*c]cf64, incx: isize, tau: [*c]cf64) isize;
pub const slarfg_work = LAPACKE_slarfg_work;
pub const dlarfg_work = LAPACKE_dlarfg_work;
pub const clarfg_work = LAPACKE_clarfg_work;
pub const zlarfg_work = LAPACKE_zlarfg_work;

extern fn LAPACKE_slarft_work(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const f32, ldv: isize, tau: [*c]const f32, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dlarft_work(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const f64, ldv: isize, tau: [*c]const f64, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_clarft_work(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const cf32, ldv: isize, tau: [*c]const cf32, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zlarft_work(layout: c_int, direct: u8, storev: u8, n: isize, k: isize, v: [*c]const cf64, ldv: isize, tau: [*c]const cf64, t: [*c]cf64, ldt: isize) isize;
pub const slarft_work = LAPACKE_slarft_work;
pub const dlarft_work = LAPACKE_dlarft_work;
pub const clarft_work = LAPACKE_clarft_work;
pub const zlarft_work = LAPACKE_zlarft_work;

extern fn LAPACKE_slarfx_work(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const f32, tau: f32, c: [*c]f32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_dlarfx_work(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const f64, tau: f64, c: [*c]f64, ldc: isize, work: [*c]f64) isize;
extern fn LAPACKE_clarfx_work(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const cf32, tau: cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zlarfx_work(layout: c_int, side: u8, m: isize, n: isize, v: [*c]const cf64, tau: cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64) isize;
pub const slarfx_work = LAPACKE_slarfx_work;
pub const dlarfx_work = LAPACKE_dlarfx_work;
pub const clarfx_work = LAPACKE_clarfx_work;
pub const zlarfx_work = LAPACKE_zlarfx_work;

extern fn LAPACKE_slarnv_work(idist: isize, iseed: [*c]isize, n: isize, x: [*c]f32) isize;
extern fn LAPACKE_dlarnv_work(idist: isize, iseed: [*c]isize, n: isize, x: [*c]f64) isize;
extern fn LAPACKE_clarnv_work(idist: isize, iseed: [*c]isize, n: isize, x: [*c]cf32) isize;
extern fn LAPACKE_zlarnv_work(idist: isize, iseed: [*c]isize, n: isize, x: [*c]cf64) isize;
pub const slarnv_work = LAPACKE_slarnv_work;
pub const dlarnv_work = LAPACKE_dlarnv_work;
pub const clarnv_work = LAPACKE_clarnv_work;
pub const zlarnv_work = LAPACKE_zlarnv_work;

extern fn LAPACKE_slascl_work(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f32, cto: f32, m: isize, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlascl_work(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f64, cto: f64, m: isize, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_clascl_work(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f32, cto: f32, m: isize, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlascl_work(layout: c_int, @"type": u8, kl: isize, ku: isize, cfrom: f64, cto: f64, m: isize, n: isize, a: [*c]cf64, lda: isize) isize;
pub const slascl_work = LAPACKE_slascl_work;
pub const dlascl_work = LAPACKE_dlascl_work;
pub const clascl_work = LAPACKE_clascl_work;
pub const zlascl_work = LAPACKE_zlascl_work;

extern fn LAPACKE_slaset_work(layout: c_int, uplo: u8, m: isize, n: isize, alpha: f32, beta: f32, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlaset_work(layout: c_int, uplo: u8, m: isize, n: isize, alpha: f64, beta: f64, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_claset_work(layout: c_int, uplo: u8, m: isize, n: isize, alpha: cf32, beta: cf32, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlaset_work(layout: c_int, uplo: u8, m: isize, n: isize, alpha: cf64, beta: cf64, a: [*c]cf64, lda: isize) isize;
pub const slaset_work = LAPACKE_slaset_work;
pub const dlaset_work = LAPACKE_dlaset_work;
pub const claset_work = LAPACKE_claset_work;
pub const zlaset_work = LAPACKE_zlaset_work;

extern fn LAPACKE_slasrt_work(id: u8, n: isize, d: [*c]f32) isize;
extern fn LAPACKE_dlasrt_work(id: u8, n: isize, d: [*c]f64) isize;
pub const slasrt_work = LAPACKE_slasrt_work;
pub const dlasrt_work = LAPACKE_dlasrt_work;

extern fn LAPACKE_slassq_work(n: isize, x: [*c]f32, incx: isize, scale: [*c]f32, sumsq: [*c]f32) isize;
extern fn LAPACKE_dlassq_work(n: isize, x: [*c]f64, incx: isize, scale: [*c]f64, sumsq: [*c]f64) isize;
extern fn LAPACKE_classq_work(n: isize, x: [*c]cf32, incx: isize, scale: [*c]f32, sumsq: [*c]f32) isize;
extern fn LAPACKE_zlassq_work(n: isize, x: [*c]cf64, incx: isize, scale: [*c]f64, sumsq: [*c]f64) isize;
pub const slassq_work = LAPACKE_slassq_work;
pub const dlassq_work = LAPACKE_dlassq_work;
pub const classq_work = LAPACKE_classq_work;
pub const zlassq_work = LAPACKE_zlassq_work;

extern fn LAPACKE_slaswp_work(layout: c_int, n: isize, a: [*c]f32, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
extern fn LAPACKE_dlaswp_work(layout: c_int, n: isize, a: [*c]f64, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
extern fn LAPACKE_claswp_work(layout: c_int, n: isize, a: [*c]cf32, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
extern fn LAPACKE_zlaswp_work(layout: c_int, n: isize, a: [*c]cf64, lda: isize, k1: isize, k2: isize, ipiv: [*c]const isize, incx: isize) isize;
pub const slaswp_work = LAPACKE_slaswp_work;
pub const dlaswp_work = LAPACKE_dlaswp_work;
pub const claswp_work = LAPACKE_claswp_work;
pub const zlaswp_work = LAPACKE_zlaswp_work;

extern fn LAPACKE_slatms_work(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f32, mode: isize, cond: f32, dmax: f32, kl: isize, ku: isize, pack: u8, a: [*c]f32, lda: isize, work: [*c]f32) isize;
extern fn LAPACKE_dlatms_work(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f64, mode: isize, cond: f64, dmax: f64, kl: isize, ku: isize, pack: u8, a: [*c]f64, lda: isize, work: [*c]f64) isize;
extern fn LAPACKE_clatms_work(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f32, mode: isize, cond: f32, dmax: f32, kl: isize, ku: isize, pack: u8, a: [*c]cf32, lda: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zlatms_work(layout: c_int, m: isize, n: isize, dist: u8, iseed: [*c]isize, sym: u8, d: [*c]f64, mode: isize, cond: f64, dmax: f64, kl: isize, ku: isize, pack: u8, a: [*c]cf64, lda: isize, work: [*c]cf64) isize;
pub const slatms_work = LAPACKE_slatms_work;
pub const dlatms_work = LAPACKE_dlatms_work;
pub const clatms_work = LAPACKE_clatms_work;
pub const zlatms_work = LAPACKE_zlatms_work;

extern fn LAPACKE_slauum_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dlauum_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_clauum_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zlauum_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const slauum_work = LAPACKE_slauum_work;
pub const dlauum_work = LAPACKE_dlauum_work;
pub const clauum_work = LAPACKE_clauum_work;
pub const zlauum_work = LAPACKE_zlauum_work;

extern fn LAPACKE_sopgtr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, tau: [*c]const f32, q: [*c]f32, ldq: isize, work: [*c]f32) isize;
extern fn LAPACKE_dopgtr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, tau: [*c]const f64, q: [*c]f64, ldq: isize, work: [*c]f64) isize;
pub const sopgtr_work = LAPACKE_sopgtr_work;
pub const dopgtr_work = LAPACKE_dopgtr_work;

extern fn LAPACKE_sopmtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const f32, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_dopmtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const f64, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64) isize;
pub const sopmtr_work = LAPACKE_sopmtr_work;
pub const dopmtr_work = LAPACKE_dopmtr_work;

extern fn LAPACKE_sorgbr_work(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorgbr_work(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorgbr_work = LAPACKE_sorgbr_work;
pub const dorgbr_work = LAPACKE_dorgbr_work;

extern fn LAPACKE_sorghr_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorghr_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorghr_work = LAPACKE_sorghr_work;
pub const dorghr_work = LAPACKE_dorghr_work;

extern fn LAPACKE_sorglq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorglq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorglq_work = LAPACKE_sorglq_work;
pub const dorglq_work = LAPACKE_dorglq_work;

extern fn LAPACKE_sorgql_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorgql_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorgql_work = LAPACKE_sorgql_work;
pub const dorgql_work = LAPACKE_dorgql_work;

extern fn LAPACKE_sorgqr_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorgqr_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorgqr_work = LAPACKE_sorgqr_work;
pub const dorgqr_work = LAPACKE_dorgqr_work;

extern fn LAPACKE_sorgrq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorgrq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorgrq_work = LAPACKE_sorgrq_work;
pub const dorgrq_work = LAPACKE_dorgrq_work;

extern fn LAPACKE_sorgtr_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, tau: [*c]const f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorgtr_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, tau: [*c]const f64, work: [*c]f64, lwork: isize) isize;
pub const sorgtr_work = LAPACKE_sorgtr_work;
pub const dorgtr_work = LAPACKE_dorgtr_work;

extern fn LAPACKE_sorgtsqr_row_work(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]f32, lda: isize, t: [*c]const f32, ldt: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dorgtsqr_row_work(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]f64, lda: isize, t: [*c]const f64, ldt: isize, work: [*c]f64, lwork: isize) isize;
pub const sorgtsqr_row_work = LAPACKE_sorgtsqr_row_work;
pub const dorgtsqr_row_work = LAPACKE_dorgtsqr_row_work;

extern fn LAPACKE_sormbr_work(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormbr_work(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormbr_work = LAPACKE_sormbr_work;
pub const dormbr_work = LAPACKE_dormbr_work;

extern fn LAPACKE_sormhr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormhr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormhr_work = LAPACKE_sormhr_work;
pub const dormhr_work = LAPACKE_dormhr_work;

extern fn LAPACKE_sormlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormlq_work = LAPACKE_sormlq_work;
pub const dormlq_work = LAPACKE_dormlq_work;

extern fn LAPACKE_sormql_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormql_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormql_work = LAPACKE_sormql_work;
pub const dormql_work = LAPACKE_dormql_work;

extern fn LAPACKE_sormqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormqr_work = LAPACKE_sormqr_work;
pub const dormqr_work = LAPACKE_dormqr_work;

extern fn LAPACKE_sormrq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormrq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormrq_work = LAPACKE_sormrq_work;
pub const dormrq_work = LAPACKE_dormrq_work;

extern fn LAPACKE_sormrz_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormrz_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormrz_work = LAPACKE_sormrz_work;
pub const dormrz_work = LAPACKE_dormrz_work;

extern fn LAPACKE_sormtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const f32, lda: isize, tau: [*c]const f32, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dormtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const f64, lda: isize, tau: [*c]const f64, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
pub const sormtr_work = LAPACKE_sormtr_work;
pub const dormtr_work = LAPACKE_dormtr_work;

extern fn LAPACKE_spbcon_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f32, ldab: isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dpbcon_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f64, ldab: isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cpbcon_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf32, ldab: isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zpbcon_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf64, ldab: isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const spbcon_work = LAPACKE_spbcon_work;
pub const dpbcon_work = LAPACKE_dpbcon_work;
pub const cpbcon_work = LAPACKE_cpbcon_work;
pub const zpbcon_work = LAPACKE_zpbcon_work;

extern fn LAPACKE_spbequ_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f32, ldab: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dpbequ_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const f64, ldab: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cpbequ_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf32, ldab: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zpbequ_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]const cf64, ldab: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const spbequ_work = LAPACKE_spbequ_work;
pub const dpbequ_work = LAPACKE_dpbequ_work;
pub const cpbequ_work = LAPACKE_cpbequ_work;
pub const zpbequ_work = LAPACKE_zpbequ_work;

extern fn LAPACKE_spbrfs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, afb: [*c]const f32, ldafb: isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dpbrfs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, afb: [*c]const f64, ldafb: isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cpbrfs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, afb: [*c]const cf32, ldafb: isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zpbrfs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, afb: [*c]const cf64, ldafb: isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const spbrfs_work = LAPACKE_spbrfs_work;
pub const dpbrfs_work = LAPACKE_dpbrfs_work;
pub const cpbrfs_work = LAPACKE_cpbrfs_work;
pub const zpbrfs_work = LAPACKE_zpbrfs_work;

extern fn LAPACKE_spbstf_work(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]f32, ldbb: isize) isize;
extern fn LAPACKE_dpbstf_work(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]f64, ldbb: isize) isize;
extern fn LAPACKE_cpbstf_work(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]cf32, ldbb: isize) isize;
extern fn LAPACKE_zpbstf_work(layout: c_int, uplo: u8, n: isize, kb: isize, bb: [*c]cf64, ldbb: isize) isize;
pub const spbstf_work = LAPACKE_spbstf_work;
pub const dpbstf_work = LAPACKE_dpbstf_work;
pub const cpbstf_work = LAPACKE_cpbstf_work;
pub const zpbstf_work = LAPACKE_zpbstf_work;

extern fn LAPACKE_spbsv_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f32, ldab: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpbsv_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f64, ldab: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpbsv_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpbsv_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, b: [*c]cf64, ldb: isize) isize;
pub const spbsv_work = LAPACKE_spbsv_work;
pub const dpbsv_work = LAPACKE_dpbsv_work;
pub const cpbsv_work = LAPACKE_cpbsv_work;
pub const zpbsv_work = LAPACKE_zpbsv_work;

extern fn LAPACKE_spbsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f32, ldab: isize, afb: [*c]f32, ldafb: isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dpbsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]f64, ldab: isize, afb: [*c]f64, ldafb: isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cpbsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf32, ldab: isize, afb: [*c]cf32, ldafb: isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zpbsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]cf64, ldab: isize, afb: [*c]cf64, ldafb: isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const spbsvx_work = LAPACKE_spbsvx_work;
pub const dpbsvx_work = LAPACKE_dpbsvx_work;
pub const cpbsvx_work = LAPACKE_cpbsvx_work;
pub const zpbsvx_work = LAPACKE_zpbsvx_work;

extern fn LAPACKE_spbtrf_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize) isize;
extern fn LAPACKE_dpbtrf_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize) isize;
extern fn LAPACKE_cpbtrf_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize) isize;
extern fn LAPACKE_zpbtrf_work(layout: c_int, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize) isize;
pub const spbtrf_work = LAPACKE_spbtrf_work;
pub const dpbtrf_work = LAPACKE_dpbtrf_work;
pub const cpbtrf_work = LAPACKE_cpbtrf_work;
pub const zpbtrf_work = LAPACKE_zpbtrf_work;

extern fn LAPACKE_spbtrs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpbtrs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpbtrs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpbtrs_work(layout: c_int, uplo: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, b: [*c]cf64, ldb: isize) isize;
pub const spbtrs_work = LAPACKE_spbtrs_work;
pub const dpbtrs_work = LAPACKE_dpbtrs_work;
pub const cpbtrs_work = LAPACKE_cpbtrs_work;
pub const zpbtrs_work = LAPACKE_zpbtrs_work;

extern fn LAPACKE_spftrf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f32) isize;
extern fn LAPACKE_dpftrf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f64) isize;
extern fn LAPACKE_cpftrf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf32) isize;
extern fn LAPACKE_zpftrf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf64) isize;
pub const spftrf_work = LAPACKE_spftrf_work;
pub const dpftrf_work = LAPACKE_dpftrf_work;
pub const cpftrf_work = LAPACKE_cpftrf_work;
pub const zpftrf_work = LAPACKE_zpftrf_work;

extern fn LAPACKE_spftri_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f32) isize;
extern fn LAPACKE_dpftri_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]f64) isize;
extern fn LAPACKE_cpftri_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf32) isize;
extern fn LAPACKE_zpftri_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]cf64) isize;
pub const spftri_work = LAPACKE_spftri_work;
pub const dpftri_work = LAPACKE_dpftri_work;
pub const cpftri_work = LAPACKE_cpftri_work;
pub const zpftri_work = LAPACKE_zpftri_work;

extern fn LAPACKE_spftrs_work(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpftrs_work(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpftrs_work(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpftrs_work(layout: c_int, transr: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const spftrs_work = LAPACKE_spftrs_work;
pub const dpftrs_work = LAPACKE_dpftrs_work;
pub const cpftrs_work = LAPACKE_cpftrs_work;
pub const zpftrs_work = LAPACKE_zpftrs_work;

extern fn LAPACKE_spocon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dpocon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cpocon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zpocon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const spocon_work = LAPACKE_spocon_work;
pub const dpocon_work = LAPACKE_dpocon_work;
pub const cpocon_work = LAPACKE_cpocon_work;
pub const zpocon_work = LAPACKE_zpocon_work;

extern fn LAPACKE_spoequ_work(layout: c_int, n: isize, a: [*c]const f32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dpoequ_work(layout: c_int, n: isize, a: [*c]const f64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cpoequ_work(layout: c_int, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zpoequ_work(layout: c_int, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const spoequ_work = LAPACKE_spoequ_work;
pub const dpoequ_work = LAPACKE_dpoequ_work;
pub const cpoequ_work = LAPACKE_cpoequ_work;
pub const zpoequ_work = LAPACKE_zpoequ_work;

extern fn LAPACKE_spoequb_work(layout: c_int, n: isize, a: [*c]const f32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dpoequb_work(layout: c_int, n: isize, a: [*c]const f64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cpoequb_work(layout: c_int, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zpoequb_work(layout: c_int, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const spoequb_work = LAPACKE_spoequb_work;
pub const dpoequb_work = LAPACKE_dpoequb_work;
pub const cpoequb_work = LAPACKE_cpoequb_work;
pub const zpoequb_work = LAPACKE_zpoequb_work;

extern fn LAPACKE_sporfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dporfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cporfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zporfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sporfs_work = LAPACKE_sporfs_work;
pub const dporfs_work = LAPACKE_dporfs_work;
pub const cporfs_work = LAPACKE_cporfs_work;
pub const zporfs_work = LAPACKE_zporfs_work;

extern fn LAPACKE_sporfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, s: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dporfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, s: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cporfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, s: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zporfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, s: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sporfsx_work = LAPACKE_sporfsx_work;
pub const dporfsx_work = LAPACKE_dporfsx_work;
pub const cporfsx_work = LAPACKE_cporfsx_work;
pub const zporfsx_work = LAPACKE_zporfsx_work;

extern fn LAPACKE_sposv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dposv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cposv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zposv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
extern fn LAPACKE_dsposv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, work: [*c]f64, swork: [*c]f32, iter: [*c]isize) isize;
extern fn LAPACKE_zcposv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, work: [*c]cf64, swork: [*c]cf32, rwork: [*c]f64, iter: [*c]isize) isize;
pub const sposv_work = LAPACKE_sposv_work;
pub const dposv_work = LAPACKE_dposv_work;
pub const cposv_work = LAPACKE_cposv_work;
pub const zposv_work = LAPACKE_zposv_work;
pub const dsposv_work = LAPACKE_dsposv_work;
pub const zcposv_work = LAPACKE_zcposv_work;

extern fn LAPACKE_sposvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dposvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cposvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zposvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sposvx_work = LAPACKE_sposvx_work;
pub const dposvx_work = LAPACKE_dposvx_work;
pub const cposvx_work = LAPACKE_cposvx_work;
pub const zposvx_work = LAPACKE_zposvx_work;

extern fn LAPACKE_sposvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dposvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cposvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zposvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sposvxx_work = LAPACKE_sposvxx_work;
pub const dposvxx_work = LAPACKE_dposvxx_work;
pub const cposvxx_work = LAPACKE_cposvxx_work;
pub const zposvxx_work = LAPACKE_zposvxx_work;

extern fn LAPACKE_spotrf2_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dpotrf2_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_cpotrf2_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zpotrf2_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const spotrf2_work = LAPACKE_spotrf2_work;
pub const dpotrf2_work = LAPACKE_dpotrf2_work;
pub const cpotrf2_work = LAPACKE_cpotrf2_work;
pub const zpotrf2_work = LAPACKE_zpotrf2_work;

extern fn LAPACKE_spotrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dpotrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_cpotrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zpotrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const spotrf_work = LAPACKE_spotrf_work;
pub const dpotrf_work = LAPACKE_dpotrf_work;
pub const cpotrf_work = LAPACKE_cpotrf_work;
pub const zpotrf_work = LAPACKE_zpotrf_work;

extern fn LAPACKE_spotri_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dpotri_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_cpotri_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zpotri_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const spotri_work = LAPACKE_spotri_work;
pub const dpotri_work = LAPACKE_dpotri_work;
pub const cpotri_work = LAPACKE_cpotri_work;
pub const zpotri_work = LAPACKE_zpotri_work;

extern fn LAPACKE_spotrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpotrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpotrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpotrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const spotrs_work = LAPACKE_spotrs_work;
pub const dpotrs_work = LAPACKE_dpotrs_work;
pub const cpotrs_work = LAPACKE_cpotrs_work;
pub const zpotrs_work = LAPACKE_zpotrs_work;

extern fn LAPACKE_sppcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dppcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cppcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, anorm: f32, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zppcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, anorm: f64, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sppcon_work = LAPACKE_sppcon_work;
pub const dppcon_work = LAPACKE_dppcon_work;
pub const cppcon_work = LAPACKE_cppcon_work;
pub const zppcon_work = LAPACKE_zppcon_work;

extern fn LAPACKE_sppequ_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_dppequ_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
extern fn LAPACKE_cppequ_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, s: [*c]f32, scond: [*c]f32, amax: [*c]f32) isize;
extern fn LAPACKE_zppequ_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, s: [*c]f64, scond: [*c]f64, amax: [*c]f64) isize;
pub const sppequ_work = LAPACKE_sppequ_work;
pub const dppequ_work = LAPACKE_dppequ_work;
pub const cppequ_work = LAPACKE_cppequ_work;
pub const zppequ_work = LAPACKE_zppequ_work;

extern fn LAPACKE_spprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, afp: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dpprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, afp: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cpprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]const cf32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zpprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]const cf64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const spprfs_work = LAPACKE_spprfs_work;
pub const dpprfs_work = LAPACKE_dpprfs_work;
pub const cpprfs_work = LAPACKE_cpprfs_work;
pub const zpprfs_work = LAPACKE_zpprfs_work;

extern fn LAPACKE_sppsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dppsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cppsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zppsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, b: [*c]cf64, ldb: isize) isize;
pub const sppsv_work = LAPACKE_sppsv_work;
pub const dppsv_work = LAPACKE_dppsv_work;
pub const cppsv_work = LAPACKE_cppsv_work;
pub const zppsv_work = LAPACKE_zppsv_work;

extern fn LAPACKE_sppsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]f32, afp: [*c]f32, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dppsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]f64, afp: [*c]f64, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cppsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, afp: [*c]cf32, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zppsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, afp: [*c]cf64, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sppsvx_work = LAPACKE_sppsvx_work;
pub const dppsvx_work = LAPACKE_dppsvx_work;
pub const cppsvx_work = LAPACKE_cppsvx_work;
pub const zppsvx_work = LAPACKE_zppsvx_work;

extern fn LAPACKE_spptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dpptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f64) isize;
extern fn LAPACKE_cpptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_zpptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64) isize;
pub const spptrf_work = LAPACKE_spptrf_work;
pub const dpptrf_work = LAPACKE_dpptrf_work;
pub const cpptrf_work = LAPACKE_cpptrf_work;
pub const zpptrf_work = LAPACKE_zpptrf_work;

extern fn LAPACKE_spptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dpptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f64) isize;
extern fn LAPACKE_cpptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_zpptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64) isize;
pub const spptri_work = LAPACKE_spptri_work;
pub const dpptri_work = LAPACKE_dpptri_work;
pub const cpptri_work = LAPACKE_cpptri_work;
pub const zpptri_work = LAPACKE_zpptri_work;

extern fn LAPACKE_spptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const spptrs_work = LAPACKE_spptrs_work;
pub const dpptrs_work = LAPACKE_dpptrs_work;
pub const cpptrs_work = LAPACKE_cpptrs_work;
pub const zpptrs_work = LAPACKE_zpptrs_work;

extern fn LAPACKE_spstrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f32, work: [*c]f32) isize;
extern fn LAPACKE_dpstrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f64, work: [*c]f64) isize;
extern fn LAPACKE_cpstrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f32, work: [*c]f32) isize;
extern fn LAPACKE_zpstrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, piv: [*c]isize, rank: [*c]isize, tol: f64, work: [*c]f64) isize;
pub const spstrf_work = LAPACKE_spstrf_work;
pub const dpstrf_work = LAPACKE_dpstrf_work;
pub const cpstrf_work = LAPACKE_cpstrf_work;
pub const zpstrf_work = LAPACKE_zpstrf_work;

extern fn LAPACKE_sptcon_work(n: isize, d: [*c]const f32, e: [*c]const f32, anorm: f32, rcond: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dptcon_work(n: isize, d: [*c]const f64, e: [*c]const f64, anorm: f64, rcond: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cptcon_work(n: isize, d: [*c]const f32, e: [*c]const cf32, anorm: f32, rcond: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_zptcon_work(n: isize, d: [*c]const f64, e: [*c]const cf64, anorm: f64, rcond: [*c]f64, work: [*c]f64) isize;
pub const sptcon_work = LAPACKE_sptcon_work;
pub const dptcon_work = LAPACKE_dptcon_work;
pub const cptcon_work = LAPACKE_cptcon_work;
pub const zptcon_work = LAPACKE_zptcon_work;

extern fn LAPACKE_spteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dpteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
extern fn LAPACKE_cpteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_zpteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]f64) isize;
pub const spteqr_work = LAPACKE_spteqr_work;
pub const dpteqr_work = LAPACKE_dpteqr_work;
pub const cpteqr_work = LAPACKE_cpteqr_work;
pub const zpteqr_work = LAPACKE_zpteqr_work;

extern fn LAPACKE_sptrfs_work(layout: c_int, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const f32, df: [*c]const f32, ef: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dptrfs_work(layout: c_int, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const f64, df: [*c]const f64, ef: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cptrfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const cf32, df: [*c]const f32, ef: [*c]const cf32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zptrfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const cf64, df: [*c]const f64, ef: [*c]const cf64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sptrfs_work = LAPACKE_sptrfs_work;
pub const dptrfs_work = LAPACKE_dptrfs_work;
pub const cptrfs_work = LAPACKE_cptrfs_work;
pub const zptrfs_work = LAPACKE_zptrfs_work;

extern fn LAPACKE_sptsv_work(layout: c_int, n: isize, nrhs: isize, d: [*c]f32, e: [*c]f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dptsv_work(layout: c_int, n: isize, nrhs: isize, d: [*c]f64, e: [*c]f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cptsv_work(layout: c_int, n: isize, nrhs: isize, d: [*c]f32, e: [*c]cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zptsv_work(layout: c_int, n: isize, nrhs: isize, d: [*c]f64, e: [*c]cf64, b: [*c]cf64, ldb: isize) isize;
pub const sptsv_work = LAPACKE_sptsv_work;
pub const dptsv_work = LAPACKE_dptsv_work;
pub const cptsv_work = LAPACKE_cptsv_work;
pub const zptsv_work = LAPACKE_zptsv_work;

extern fn LAPACKE_sptsvx_work(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const f32, df: [*c]f32, ef: [*c]f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dptsvx_work(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const f64, df: [*c]f64, ef: [*c]f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_cptsvx_work(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const cf32, df: [*c]f32, ef: [*c]cf32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zptsvx_work(layout: c_int, fact: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const cf64, df: [*c]f64, ef: [*c]cf64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sptsvx_work = LAPACKE_sptsvx_work;
pub const dptsvx_work = LAPACKE_dptsvx_work;
pub const cptsvx_work = LAPACKE_cptsvx_work;
pub const zptsvx_work = LAPACKE_zptsvx_work;

extern fn LAPACKE_spttrf_work(n: isize, d: [*c]f32, e: [*c]f32) isize;
extern fn LAPACKE_dpttrf_work(n: isize, d: [*c]f64, e: [*c]f64) isize;
extern fn LAPACKE_cpttrf_work(n: isize, d: [*c]f32, e: [*c]cf32) isize;
extern fn LAPACKE_zpttrf_work(n: isize, d: [*c]f64, e: [*c]cf64) isize;
pub const spttrf_work = LAPACKE_spttrf_work;
pub const dpttrf_work = LAPACKE_dpttrf_work;
pub const cpttrf_work = LAPACKE_cpttrf_work;
pub const zpttrf_work = LAPACKE_zpttrf_work;

extern fn LAPACKE_spttrs_work(layout: c_int, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dpttrs_work(layout: c_int, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cpttrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f32, e: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zpttrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, d: [*c]const f64, e: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const spttrs_work = LAPACKE_spttrs_work;
pub const dpttrs_work = LAPACKE_dpttrs_work;
pub const cpttrs_work = LAPACKE_cpttrs_work;
pub const zpttrs_work = LAPACKE_zpttrs_work;

extern fn LAPACKE_ssbev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dsbev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
pub const ssbev_work = LAPACKE_ssbev_work;
pub const dsbev_work = LAPACKE_dsbev_work;

extern fn LAPACKE_ssbevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsbevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssbevd_work = LAPACKE_ssbevd_work;
pub const dsbevd_work = LAPACKE_dsbevd_work;

extern fn LAPACKE_ssbevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, q: [*c]f32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsbevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, q: [*c]f64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const ssbevx_work = LAPACKE_ssbevx_work;
pub const dsbevx_work = LAPACKE_dsbevx_work;

extern fn LAPACKE_ssbgst_work(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]const f32, ldbb: isize, x: [*c]f32, ldx: isize, work: [*c]f32) isize;
extern fn LAPACKE_dsbgst_work(layout: c_int, vect: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]const f64, ldbb: isize, x: [*c]f64, ldx: isize, work: [*c]f64) isize;
pub const ssbgst_work = LAPACKE_ssbgst_work;
pub const dsbgst_work = LAPACKE_dsbgst_work;

extern fn LAPACKE_ssbgv_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]f32, ldbb: isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dsbgv_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]f64, ldbb: isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
pub const ssbgv_work = LAPACKE_ssbgv_work;
pub const dsbgv_work = LAPACKE_dsbgv_work;

extern fn LAPACKE_ssbgvd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]f32, ldbb: isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsbgvd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]f64, ldbb: isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssbgvd_work = LAPACKE_ssbgvd_work;
pub const dsbgvd_work = LAPACKE_dsbgvd_work;

extern fn LAPACKE_ssbgvx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f32, ldab: isize, bb: [*c]f32, ldbb: isize, q: [*c]f32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsbgvx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ka: isize, kb: isize, ab: [*c]f64, ldab: isize, bb: [*c]f64, ldbb: isize, q: [*c]f64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const ssbgvx_work = LAPACKE_ssbgvx_work;
pub const dsbgvx_work = LAPACKE_dsbgvx_work;

extern fn LAPACKE_ssbtrd_work(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, d: [*c]f32, e: [*c]f32, q: [*c]f32, ldq: isize, work: [*c]f32) isize;
extern fn LAPACKE_dsbtrd_work(layout: c_int, vect: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, d: [*c]f64, e: [*c]f64, q: [*c]f64, ldq: isize, work: [*c]f64) isize;
pub const ssbtrd_work = LAPACKE_ssbtrd_work;
pub const dsbtrd_work = LAPACKE_dsbtrd_work;

extern fn LAPACKE_ssfrk_work(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f32, a: [*c]const f32, lda: isize, beta: f32, c: [*c]f32) isize;
extern fn LAPACKE_dsfrk_work(layout: c_int, transr: u8, uplo: u8, trans: u8, n: isize, k: isize, alpha: f64, a: [*c]const f64, lda: isize, beta: f64, c: [*c]f64) isize;
pub const ssfrk_work = LAPACKE_ssfrk_work;
pub const dsfrk_work = LAPACKE_dsfrk_work;

extern fn LAPACKE_sspcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dspcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cspcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zspcon_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const sspcon_work = LAPACKE_sspcon_work;
pub const dspcon_work = LAPACKE_dspcon_work;
pub const cspcon_work = LAPACKE_cspcon_work;
pub const zspcon_work = LAPACKE_zspcon_work;

extern fn LAPACKE_sspev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dspev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
pub const sspev_work = LAPACKE_sspev_work;
pub const dspev_work = LAPACKE_dspev_work;

extern fn LAPACKE_sspevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dspevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sspevd_work = LAPACKE_sspevd_work;
pub const dspevd_work = LAPACKE_dspevd_work;

extern fn LAPACKE_sspevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dspevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const sspevx_work = LAPACKE_sspevx_work;
pub const dspevx_work = LAPACKE_dspevx_work;

extern fn LAPACKE_sspgst_work(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]const f32) isize;
extern fn LAPACKE_dspgst_work(layout: c_int, itype: isize, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]const f64) isize;
pub const sspgst_work = LAPACKE_sspgst_work;
pub const dspgst_work = LAPACKE_dspgst_work;

extern fn LAPACKE_sspgv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dspgv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
pub const sspgv_work = LAPACKE_sspgv_work;
pub const dspgv_work = LAPACKE_dspgv_work;

extern fn LAPACKE_sspgvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]f32, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dspgvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]f64, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sspgvd_work = LAPACKE_sspgvd_work;
pub const dspgvd_work = LAPACKE_dspgvd_work;

extern fn LAPACKE_sspgvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f32, bp: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dspgvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, ap: [*c]f64, bp: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const sspgvx_work = LAPACKE_sspgvx_work;
pub const dspgvx_work = LAPACKE_dspgvx_work;

extern fn LAPACKE_ssprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, afp: [*c]const f32, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dsprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, afp: [*c]const f64, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_csprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]const cf32, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zsprfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]const cf64, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const ssprfs_work = LAPACKE_ssprfs_work;
pub const dsprfs_work = LAPACKE_dsprfs_work;
pub const csprfs_work = LAPACKE_csprfs_work;
pub const zsprfs_work = LAPACKE_zsprfs_work;

extern fn LAPACKE_sspsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f32, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dspsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]f64, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cspsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zspsv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const sspsv_work = LAPACKE_sspsv_work;
pub const dspsv_work = LAPACKE_dspsv_work;
pub const cspsv_work = LAPACKE_cspsv_work;
pub const zspsv_work = LAPACKE_zspsv_work;

extern fn LAPACKE_sspsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, afp: [*c]f32, ipiv: [*c]isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dspsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, afp: [*c]f64, ipiv: [*c]isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_cspsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, afp: [*c]cf32, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zspsvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, afp: [*c]cf64, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const sspsvx_work = LAPACKE_sspsvx_work;
pub const dspsvx_work = LAPACKE_dspsvx_work;
pub const cspsvx_work = LAPACKE_cspsvx_work;
pub const zspsvx_work = LAPACKE_zspsvx_work;

extern fn LAPACKE_ssptrd_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f32, d: [*c]f32, e: [*c]f32, tau: [*c]f32) isize;
extern fn LAPACKE_dsptrd_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f64, d: [*c]f64, e: [*c]f64, tau: [*c]f64) isize;
pub const ssptrd_work = LAPACKE_ssptrd_work;
pub const dsptrd_work = LAPACKE_dsptrd_work;

extern fn LAPACKE_ssptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f32, ipiv: [*c]isize) isize;
extern fn LAPACKE_dsptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f64, ipiv: [*c]isize) isize;
extern fn LAPACKE_csptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zsptrf_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]isize) isize;
pub const ssptrf_work = LAPACKE_ssptrf_work;
pub const dsptrf_work = LAPACKE_dsptrf_work;
pub const csptrf_work = LAPACKE_csptrf_work;
pub const zsptrf_work = LAPACKE_zsptrf_work;

extern fn LAPACKE_ssptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f32, ipiv: [*c]const isize, work: [*c]f32) isize;
extern fn LAPACKE_dsptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]f64, ipiv: [*c]const isize, work: [*c]f64) isize;
extern fn LAPACKE_csptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf32, ipiv: [*c]const isize, work: [*c]cf32) isize;
extern fn LAPACKE_zsptri_work(layout: c_int, uplo: u8, n: isize, ap: [*c]cf64, ipiv: [*c]const isize, work: [*c]cf64) isize;
pub const ssptri_work = LAPACKE_ssptri_work;
pub const dsptri_work = LAPACKE_dsptri_work;
pub const csptri_work = LAPACKE_csptri_work;
pub const zsptri_work = LAPACKE_zsptri_work;

extern fn LAPACKE_ssptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f32, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const f64, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsptrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, ap: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const ssptrs_work = LAPACKE_ssptrs_work;
pub const dsptrs_work = LAPACKE_dsptrs_work;
pub const csptrs_work = LAPACKE_csptrs_work;
pub const zsptrs_work = LAPACKE_zsptrs_work;

extern fn LAPACKE_sstebz_work(range: u8, order: u8, n: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, d: [*c]const f32, e: [*c]const f32, m: [*c]isize, nsplit: [*c]isize, w: [*c]f32, iblock: [*c]isize, isplit: [*c]isize, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dstebz_work(range: u8, order: u8, n: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, d: [*c]const f64, e: [*c]const f64, m: [*c]isize, nsplit: [*c]isize, w: [*c]f64, iblock: [*c]isize, isplit: [*c]isize, work: [*c]f64, iwork: [*c]isize) isize;
pub const sstebz_work = LAPACKE_sstebz_work;
pub const dstebz_work = LAPACKE_dstebz_work;

extern fn LAPACKE_sstedc_work(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dstedc_work(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_cstedc_work(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zstedc_work(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sstedc_work = LAPACKE_sstedc_work;
pub const dstedc_work = LAPACKE_dstedc_work;
pub const cstedc_work = LAPACKE_cstedc_work;
pub const zstedc_work = LAPACKE_zstedc_work;

extern fn LAPACKE_sstegr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dstegr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_cstegr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, isuppz: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zstegr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, isuppz: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sstegr_work = LAPACKE_sstegr_work;
pub const dstegr_work = LAPACKE_dstegr_work;
pub const cstegr_work = LAPACKE_cstegr_work;
pub const zstegr_work = LAPACKE_zstegr_work;

extern fn LAPACKE_sstein_work(layout: c_int, n: isize, d: [*c]const f32, e: [*c]const f32, m: isize, w: [*c]const f32, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifailv: [*c]isize) isize;
extern fn LAPACKE_dstein_work(layout: c_int, n: isize, d: [*c]const f64, e: [*c]const f64, m: isize, w: [*c]const f64, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifailv: [*c]isize) isize;
extern fn LAPACKE_cstein_work(layout: c_int, n: isize, d: [*c]const f32, e: [*c]const f32, m: isize, w: [*c]const f32, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]cf32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifailv: [*c]isize) isize;
extern fn LAPACKE_zstein_work(layout: c_int, n: isize, d: [*c]const f64, e: [*c]const f64, m: isize, w: [*c]const f64, iblock: [*c]const isize, isplit: [*c]const isize, z: [*c]cf64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifailv: [*c]isize) isize;
pub const sstein_work = LAPACKE_sstein_work;
pub const dstein_work = LAPACKE_dstein_work;
pub const cstein_work = LAPACKE_cstein_work;
pub const zstein_work = LAPACKE_zstein_work;

extern fn LAPACKE_sstemr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dstemr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_cstemr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zstemr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, nzc: isize, isuppz: [*c]isize, tryrac: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sstemr_work = LAPACKE_sstemr_work;
pub const dstemr_work = LAPACKE_dstemr_work;
pub const cstemr_work = LAPACKE_cstemr_work;
pub const zstemr_work = LAPACKE_zstemr_work;

extern fn LAPACKE_ssteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dsteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
extern fn LAPACKE_csteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_zsteqr_work(layout: c_int, compz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]f64) isize;
pub const ssteqr_work = LAPACKE_ssteqr_work;
pub const dsteqr_work = LAPACKE_dsteqr_work;
pub const csteqr_work = LAPACKE_csteqr_work;
pub const zsteqr_work = LAPACKE_zsteqr_work;

extern fn LAPACKE_ssterf_work(n: isize, d: [*c]f32, e: [*c]f32) isize;
extern fn LAPACKE_dsterf_work(n: isize, d: [*c]f64, e: [*c]f64) isize;
pub const ssterf_work = LAPACKE_ssterf_work;
pub const dsterf_work = LAPACKE_dsterf_work;

extern fn LAPACKE_sstev_work(layout: c_int, jobz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32) isize;
extern fn LAPACKE_dstev_work(layout: c_int, jobz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64) isize;
pub const sstev_work = LAPACKE_sstev_work;
pub const dstev_work = LAPACKE_dstev_work;

extern fn LAPACKE_sstevd_work(layout: c_int, jobz: u8, n: isize, d: [*c]f32, e: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dstevd_work(layout: c_int, jobz: u8, n: isize, d: [*c]f64, e: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sstevd_work = LAPACKE_sstevd_work;
pub const dstevd_work = LAPACKE_dstevd_work;

extern fn LAPACKE_sstevr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dstevr_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const sstevr_work = LAPACKE_sstevr_work;
pub const dstevr_work = LAPACKE_dstevr_work;

extern fn LAPACKE_sstevx_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f32, e: [*c]f32, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dstevx_work(layout: c_int, jobz: u8, range: u8, n: isize, d: [*c]f64, e: [*c]f64, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const sstevx_work = LAPACKE_sstevx_work;
pub const dstevx_work = LAPACKE_dstevx_work;

extern fn LAPACKE_ssycon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dsycon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_csycon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zsycon_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const ssycon_work = LAPACKE_ssycon_work;
pub const dsycon_work = LAPACKE_dsycon_work;
pub const csycon_work = LAPACKE_csycon_work;
pub const zsycon_work = LAPACKE_zsycon_work;

extern fn LAPACKE_ssyequb_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32, work: [*c]f32) isize;
extern fn LAPACKE_dsyequb_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64, work: [*c]f64) isize;
extern fn LAPACKE_csyequb_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, s: [*c]f32, scond: [*c]f32, amax: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zsyequb_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, s: [*c]f64, scond: [*c]f64, amax: [*c]f64, work: [*c]cf64) isize;
pub const ssyequb_work = LAPACKE_ssyequb_work;
pub const dsyequb_work = LAPACKE_dsyequb_work;
pub const csyequb_work = LAPACKE_csyequb_work;
pub const zsyequb_work = LAPACKE_zsyequb_work;

extern fn LAPACKE_ssyev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsyev_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const ssyev_work = LAPACKE_ssyev_work;
pub const dsyev_work = LAPACKE_dsyev_work;

extern fn LAPACKE_ssyevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsyevd_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssyevd_work = LAPACKE_ssyevd_work;
pub const dsyevd_work = LAPACKE_dsyevd_work;

extern fn LAPACKE_ssyevr_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsyevr_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssyevr_work = LAPACKE_ssyevr_work;
pub const dsyevr_work = LAPACKE_dsyevr_work;

extern fn LAPACKE_ssyevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsyevx_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const ssyevx_work = LAPACKE_ssyevx_work;
pub const dsyevx_work = LAPACKE_dsyevx_work;

extern fn LAPACKE_ssygst_work(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]const f32, ldb: isize) isize;
extern fn LAPACKE_dsygst_work(layout: c_int, itype: isize, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]const f64, ldb: isize) isize;
pub const ssygst_work = LAPACKE_ssygst_work;
pub const dsygst_work = LAPACKE_dsygst_work;

extern fn LAPACKE_ssygv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, w: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsygv_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, w: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const ssygv_work = LAPACKE_ssygv_work;
pub const dsygv_work = LAPACKE_dsygv_work;

extern fn LAPACKE_ssygvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, w: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsygvd_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, w: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssygvd_work = LAPACKE_ssygvd_work;
pub const dsygvd_work = LAPACKE_dsygvd_work;

extern fn LAPACKE_ssygvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsygvx_work(layout: c_int, itype: isize, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const ssygvx_work = LAPACKE_ssygvx_work;
pub const dsygvx_work = LAPACKE_dsygvx_work;

extern fn LAPACKE_ssyrfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dsyrfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_csyrfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zsyrfs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const ssyrfs_work = LAPACKE_ssyrfs_work;
pub const dsyrfs_work = LAPACKE_dsyrfs_work;
pub const csyrfs_work = LAPACKE_csyrfs_work;
pub const zsyrfs_work = LAPACKE_zsyrfs_work;

extern fn LAPACKE_ssyrfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]const f32, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dsyrfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]const f64, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_csyrfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]const cf32, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f32, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zsyrfsx_work(layout: c_int, uplo: u8, equed: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]const cf64, ldaf: isize, ipiv: [*c]const isize, s: [*c]const f64, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const ssyrfsx_work = LAPACKE_ssyrfsx_work;
pub const dsyrfsx_work = LAPACKE_dsyrfsx_work;
pub const csyrfsx_work = LAPACKE_csyrfsx_work;
pub const zsyrfsx_work = LAPACKE_zsyrfsx_work;

extern fn LAPACKE_ssysv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsysv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csysv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsysv_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const ssysv_work = LAPACKE_ssysv_work;
pub const dsysv_work = LAPACKE_dsysv_work;
pub const csysv_work = LAPACKE_csysv_work;
pub const zsysv_work = LAPACKE_zsysv_work;

extern fn LAPACKE_ssysvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, b: [*c]const f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dsysvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, b: [*c]const f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_csysvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zsysvx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, b: [*c]const cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const ssysvx_work = LAPACKE_ssysvx_work;
pub const dsysvx_work = LAPACKE_dsysvx_work;
pub const csysvx_work = LAPACKE_csysvx_work;
pub const zsysvx_work = LAPACKE_zsysvx_work;

extern fn LAPACKE_ssysvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, af: [*c]f32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f32, b: [*c]f32, ldb: isize, x: [*c]f32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dsysvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, af: [*c]f64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f64, b: [*c]f64, ldb: isize, x: [*c]f64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_csysvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, af: [*c]cf32, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f32, b: [*c]cf32, ldb: isize, x: [*c]cf32, ldx: isize, rcond: [*c]f32, rpvgrw: [*c]f32, berr: [*c]f32, n_err_bnds: isize, err_bnds_norm: [*c]f32, err_bnds_comp: [*c]f32, nparams: isize, params: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_zsysvxx_work(layout: c_int, fact: u8, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, af: [*c]cf64, ldaf: isize, ipiv: [*c]isize, equed: [*c]u8, s: [*c]f64, b: [*c]cf64, ldb: isize, x: [*c]cf64, ldx: isize, rcond: [*c]f64, rpvgrw: [*c]f64, berr: [*c]f64, n_err_bnds: isize, err_bnds_norm: [*c]f64, err_bnds_comp: [*c]f64, nparams: isize, params: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const ssysvxx_work = LAPACKE_ssysvxx_work;
pub const dsysvxx_work = LAPACKE_dsysvxx_work;
pub const csysvxx_work = LAPACKE_csysvxx_work;
pub const zsysvxx_work = LAPACKE_zsysvxx_work;

extern fn LAPACKE_ssytrd_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, d: [*c]f32, e: [*c]f32, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsytrd_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, d: [*c]f64, e: [*c]f64, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const ssytrd_work = LAPACKE_ssytrd_work;
pub const dsytrd_work = LAPACKE_dsytrd_work;

extern fn LAPACKE_ssytrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsytrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csytrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsytrf_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const ssytrf_work = LAPACKE_ssytrf_work;
pub const dsytrf_work = LAPACKE_dsytrf_work;
pub const csytrf_work = LAPACKE_csytrf_work;
pub const zsytrf_work = LAPACKE_zsytrf_work;

extern fn LAPACKE_ssytri_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, work: [*c]f32) isize;
extern fn LAPACKE_dsytri_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, work: [*c]f64) isize;
extern fn LAPACKE_csytri_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32) isize;
extern fn LAPACKE_zsytri_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64) isize;
pub const ssytri_work = LAPACKE_ssytri_work;
pub const dsytri_work = LAPACKE_dsytri_work;
pub const csytri_work = LAPACKE_csytri_work;
pub const zsytri_work = LAPACKE_zsytri_work;

extern fn LAPACKE_ssytrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsytrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csytrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsytrs_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const ssytrs_work = LAPACKE_ssytrs_work;
pub const dsytrs_work = LAPACKE_dsytrs_work;
pub const csytrs_work = LAPACKE_csytrs_work;
pub const zsytrs_work = LAPACKE_zsytrs_work;

extern fn LAPACKE_stbcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const f32, ldab: isize, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dtbcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const f64, ldab: isize, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_ctbcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const cf32, ldab: isize, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztbcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, kd: isize, ab: [*c]const cf64, ldab: isize, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const stbcon_work = LAPACKE_stbcon_work;
pub const dtbcon_work = LAPACKE_dtbcon_work;
pub const ctbcon_work = LAPACKE_ctbcon_work;
pub const ztbcon_work = LAPACKE_ztbcon_work;

extern fn LAPACKE_stbrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, b: [*c]const f32, ldb: isize, x: [*c]const f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dtbrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, b: [*c]const f64, ldb: isize, x: [*c]const f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_ctbrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, b: [*c]const cf32, ldb: isize, x: [*c]const cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztbrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, b: [*c]const cf64, ldb: isize, x: [*c]const cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const stbrfs_work = LAPACKE_stbrfs_work;
pub const dtbrfs_work = LAPACKE_dtbrfs_work;
pub const ctbrfs_work = LAPACKE_ctbrfs_work;
pub const ztbrfs_work = LAPACKE_ztbrfs_work;

extern fn LAPACKE_stbtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f32, ldab: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtbtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const f64, ldab: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctbtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf32, ldab: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztbtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, kd: isize, nrhs: isize, ab: [*c]const cf64, ldab: isize, b: [*c]cf64, ldb: isize) isize;
pub const stbtrs_work = LAPACKE_stbtrs_work;
pub const dtbtrs_work = LAPACKE_dtbtrs_work;
pub const ctbtrs_work = LAPACKE_ctbtrs_work;
pub const ztbtrs_work = LAPACKE_ztbtrs_work;

extern fn LAPACKE_stfsm_work(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: f32, a: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtfsm_work(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: f64, a: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctfsm_work(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: cf32, a: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztfsm_work(layout: c_int, transr: u8, side: u8, uplo: u8, trans: u8, diag: u8, m: isize, n: isize, alpha: cf64, a: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const stfsm_work = LAPACKE_stfsm_work;
pub const dtfsm_work = LAPACKE_dtfsm_work;
pub const ctfsm_work = LAPACKE_ctfsm_work;
pub const ztfsm_work = LAPACKE_ztfsm_work;

extern fn LAPACKE_stftri_work(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]f32) isize;
extern fn LAPACKE_dtftri_work(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]f64) isize;
extern fn LAPACKE_ctftri_work(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]cf32) isize;
extern fn LAPACKE_ztftri_work(layout: c_int, transr: u8, uplo: u8, diag: u8, n: isize, a: [*c]cf64) isize;
pub const stftri_work = LAPACKE_stftri_work;
pub const dtftri_work = LAPACKE_dtftri_work;
pub const ctftri_work = LAPACKE_ctftri_work;
pub const ztftri_work = LAPACKE_ztftri_work;

extern fn LAPACKE_stfttp_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f32, ap: [*c]f32) isize;
extern fn LAPACKE_dtfttp_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f64, ap: [*c]f64) isize;
extern fn LAPACKE_ctfttp_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf32, ap: [*c]cf32) isize;
extern fn LAPACKE_ztfttp_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf64, ap: [*c]cf64) isize;
pub const stfttp_work = LAPACKE_stfttp_work;
pub const dtfttp_work = LAPACKE_dtfttp_work;
pub const ctfttp_work = LAPACKE_ctfttp_work;
pub const ztfttp_work = LAPACKE_ztfttp_work;

extern fn LAPACKE_stfttr_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f32, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dtfttr_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const f64, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_ctfttr_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf32, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_ztfttr_work(layout: c_int, transr: u8, uplo: u8, n: isize, arf: [*c]const cf64, a: [*c]cf64, lda: isize) isize;
pub const stfttr_work = LAPACKE_stfttr_work;
pub const dtfttr_work = LAPACKE_dtfttr_work;
pub const ctfttr_work = LAPACKE_ctfttr_work;
pub const ztfttr_work = LAPACKE_ztfttr_work;

extern fn LAPACKE_stgevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const f32, lds: isize, p: [*c]const f32, ldp: isize, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]f32) isize;
extern fn LAPACKE_dtgevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const f64, lds: isize, p: [*c]const f64, ldp: isize, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]f64) isize;
extern fn LAPACKE_ctgevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const cf32, lds: isize, p: [*c]const cf32, ldp: isize, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztgevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, s: [*c]const cf64, lds: isize, p: [*c]const cf64, ldp: isize, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const stgevc_work = LAPACKE_stgevc_work;
pub const dtgevc_work = LAPACKE_dtgevc_work;
pub const ctgevc_work = LAPACKE_ctgevc_work;
pub const ztgevc_work = LAPACKE_ztgevc_work;

extern fn LAPACKE_stgexc_work(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize, ifst: [*c]isize, ilst: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dtgexc_work(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize, ifst: [*c]isize, ilst: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_ctgexc_work(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize, ifst: isize, ilst: isize) isize;
extern fn LAPACKE_ztgexc_work(layout: c_int, wantq: isize, wantz: isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize, ifst: isize, ilst: isize) isize;
pub const stgexc_work = LAPACKE_stgexc_work;
pub const dtgexc_work = LAPACKE_dtgexc_work;
pub const ctgexc_work = LAPACKE_ctgexc_work;
pub const ztgexc_work = LAPACKE_ztgexc_work;

extern fn LAPACKE_stgsen_work(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, alphar: [*c]f32, alphai: [*c]f32, beta: [*c]f32, q: [*c]f32, ldq: isize, z: [*c]f32, ldz: isize, m: [*c]isize, pl: [*c]f32, pr: [*c]f32, dif: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dtgsen_work(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, alphar: [*c]f64, alphai: [*c]f64, beta: [*c]f64, q: [*c]f64, ldq: isize, z: [*c]f64, ldz: isize, m: [*c]isize, pl: [*c]f64, pr: [*c]f64, dif: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_ctgsen_work(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, alpha: [*c]cf32, beta: [*c]cf32, q: [*c]cf32, ldq: isize, z: [*c]cf32, ldz: isize, m: [*c]isize, pl: [*c]f32, pr: [*c]f32, dif: [*c]f32, work: [*c]cf32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_ztgsen_work(layout: c_int, ijob: isize, wantq: isize, wantz: isize, select: [*c]const isize, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, alpha: [*c]cf64, beta: [*c]cf64, q: [*c]cf64, ldq: isize, z: [*c]cf64, ldz: isize, m: [*c]isize, pl: [*c]f64, pr: [*c]f64, dif: [*c]f64, work: [*c]cf64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const stgsen_work = LAPACKE_stgsen_work;
pub const dtgsen_work = LAPACKE_dtgsen_work;
pub const ctgsen_work = LAPACKE_ctgsen_work;
pub const ztgsen_work = LAPACKE_ztgsen_work;

extern fn LAPACKE_stgsja_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, tola: f32, tolb: f32, alpha: [*c]f32, beta: [*c]f32, u: [*c]f32, ldu: isize, v: [*c]f32, ldv: isize, q: [*c]f32, ldq: isize, work: [*c]f32, ncycle: [*c]isize) isize;
extern fn LAPACKE_dtgsja_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, tola: f64, tolb: f64, alpha: [*c]f64, beta: [*c]f64, u: [*c]f64, ldu: isize, v: [*c]f64, ldv: isize, q: [*c]f64, ldq: isize, work: [*c]f64, ncycle: [*c]isize) isize;
extern fn LAPACKE_ctgsja_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, tola: f32, tolb: f32, alpha: [*c]f32, beta: [*c]f32, u: [*c]cf32, ldu: isize, v: [*c]cf32, ldv: isize, q: [*c]cf32, ldq: isize, work: [*c]cf32, ncycle: [*c]isize) isize;
extern fn LAPACKE_ztgsja_work(layout: c_int, jobu: u8, jobv: u8, jobq: u8, m: isize, p: isize, n: isize, k: isize, l: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, tola: f64, tolb: f64, alpha: [*c]f64, beta: [*c]f64, u: [*c]cf64, ldu: isize, v: [*c]cf64, ldv: isize, q: [*c]cf64, ldq: isize, work: [*c]cf64, ncycle: [*c]isize) isize;
pub const stgsja_work = LAPACKE_stgsja_work;
pub const dtgsja_work = LAPACKE_dtgsja_work;
pub const ctgsja_work = LAPACKE_ctgsja_work;
pub const ztgsja_work = LAPACKE_ztgsja_work;

extern fn LAPACKE_stgsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, vl: [*c]const f32, ldvl: isize, vr: [*c]const f32, ldvr: isize, s: [*c]f32, dif: [*c]f32, mm: isize, m: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dtgsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, vl: [*c]const f64, ldvl: isize, vr: [*c]const f64, ldvr: isize, s: [*c]f64, dif: [*c]f64, mm: isize, m: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_ctgsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, vl: [*c]const cf32, ldvl: isize, vr: [*c]const cf32, ldvr: isize, s: [*c]f32, dif: [*c]f32, mm: isize, m: [*c]isize, work: [*c]cf32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_ztgsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, vl: [*c]const cf64, ldvl: isize, vr: [*c]const cf64, ldvr: isize, s: [*c]f64, dif: [*c]f64, mm: isize, m: [*c]isize, work: [*c]cf64, lwork: isize, iwork: [*c]isize) isize;
pub const stgsna_work = LAPACKE_stgsna_work;
pub const dtgsna_work = LAPACKE_dtgsna_work;
pub const ctgsna_work = LAPACKE_ctgsna_work;
pub const ztgsna_work = LAPACKE_ztgsna_work;

extern fn LAPACKE_stgsyl_work(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]f32, ldc: isize, d: [*c]const f32, ldd: isize, e: [*c]const f32, lde: isize, f: [*c]f32, ldf: isize, scale: [*c]f32, dif: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dtgsyl_work(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]f64, ldc: isize, d: [*c]const f64, ldd: isize, e: [*c]const f64, lde: isize, f: [*c]f64, ldf: isize, scale: [*c]f64, dif: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_ctgsyl_work(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize, d: [*c]const cf32, ldd: isize, e: [*c]const cf32, lde: isize, f: [*c]cf32, ldf: isize, scale: [*c]f32, dif: [*c]f32, work: [*c]cf32, lwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_ztgsyl_work(layout: c_int, trans: u8, ijob: isize, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, d: [*c]const cf64, ldd: isize, e: [*c]const cf64, lde: isize, f: [*c]cf64, ldf: isize, scale: [*c]f64, dif: [*c]f64, work: [*c]cf64, lwork: isize, iwork: [*c]isize) isize;
pub const stgsyl_work = LAPACKE_stgsyl_work;
pub const dtgsyl_work = LAPACKE_dtgsyl_work;
pub const ctgsyl_work = LAPACKE_ctgsyl_work;
pub const ztgsyl_work = LAPACKE_ztgsyl_work;

extern fn LAPACKE_stpcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dtpcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_ctpcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const cf32, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztpcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, ap: [*c]const cf64, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const stpcon_work = LAPACKE_stpcon_work;
pub const dtpcon_work = LAPACKE_dtpcon_work;
pub const ctpcon_work = LAPACKE_ctpcon_work;
pub const ztpcon_work = LAPACKE_ztpcon_work;

extern fn LAPACKE_stprfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f32, b: [*c]const f32, ldb: isize, x: [*c]const f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dtprfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f64, b: [*c]const f64, ldb: isize, x: [*c]const f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_ctprfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf32, b: [*c]const cf32, ldb: isize, x: [*c]const cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztprfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf64, b: [*c]const cf64, ldb: isize, x: [*c]const cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const stprfs_work = LAPACKE_stprfs_work;
pub const dtprfs_work = LAPACKE_dtprfs_work;
pub const ctprfs_work = LAPACKE_ctprfs_work;
pub const ztprfs_work = LAPACKE_ztprfs_work;

extern fn LAPACKE_stptri_work(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dtptri_work(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]f64) isize;
extern fn LAPACKE_ctptri_work(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_ztptri_work(layout: c_int, uplo: u8, diag: u8, n: isize, ap: [*c]cf64) isize;
pub const stptri_work = LAPACKE_stptri_work;
pub const dtptri_work = LAPACKE_dtptri_work;
pub const ctptri_work = LAPACKE_ctptri_work;
pub const ztptri_work = LAPACKE_ztptri_work;

extern fn LAPACKE_stptrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f32, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtptrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const f64, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctptrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf32, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztptrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, ap: [*c]const cf64, b: [*c]cf64, ldb: isize) isize;
pub const stptrs_work = LAPACKE_stptrs_work;
pub const dtptrs_work = LAPACKE_dtptrs_work;
pub const ctptrs_work = LAPACKE_ctptrs_work;
pub const ztptrs_work = LAPACKE_ztptrs_work;

extern fn LAPACKE_stpttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const f32, arf: [*c]f32) isize;
extern fn LAPACKE_dtpttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const f64, arf: [*c]f64) isize;
extern fn LAPACKE_ctpttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const cf32, arf: [*c]cf32) isize;
extern fn LAPACKE_ztpttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, ap: [*c]const cf64, arf: [*c]cf64) isize;
pub const stpttf_work = LAPACKE_stpttf_work;
pub const dtpttf_work = LAPACKE_dtpttf_work;
pub const ctpttf_work = LAPACKE_ctpttf_work;
pub const ztpttf_work = LAPACKE_ztpttf_work;

extern fn LAPACKE_stpttr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f32, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dtpttr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const f64, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_ctpttr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_ztpttr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, a: [*c]cf64, lda: isize) isize;
pub const stpttr_work = LAPACKE_stpttr_work;
pub const dtpttr_work = LAPACKE_dtpttr_work;
pub const ctpttr_work = LAPACKE_ctpttr_work;
pub const ztpttr_work = LAPACKE_ztpttr_work;

extern fn LAPACKE_strcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const f32, lda: isize, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dtrcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const f64, lda: isize, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_ctrcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const cf32, lda: isize, rcond: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztrcon_work(layout: c_int, norm: u8, uplo: u8, diag: u8, n: isize, a: [*c]const cf64, lda: isize, rcond: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const strcon_work = LAPACKE_strcon_work;
pub const dtrcon_work = LAPACKE_dtrcon_work;
pub const ctrcon_work = LAPACKE_ctrcon_work;
pub const ztrcon_work = LAPACKE_ztrcon_work;

extern fn LAPACKE_strevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]isize, n: isize, t: [*c]const f32, ldt: isize, vl: [*c]f32, ldvl: isize, vr: [*c]f32, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]f32) isize;
extern fn LAPACKE_dtrevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]isize, n: isize, t: [*c]const f64, ldt: isize, vl: [*c]f64, ldvl: isize, vr: [*c]f64, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]f64) isize;
extern fn LAPACKE_ctrevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]cf32, ldt: isize, vl: [*c]cf32, ldvl: isize, vr: [*c]cf32, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztrevc_work(layout: c_int, side: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]cf64, ldt: isize, vl: [*c]cf64, ldvl: isize, vr: [*c]cf64, ldvr: isize, mm: isize, m: [*c]isize, work: [*c]cf64, rwork: [*c]f64) isize;
pub const strevc_work = LAPACKE_strevc_work;
pub const dtrevc_work = LAPACKE_dtrevc_work;
pub const ctrevc_work = LAPACKE_ctrevc_work;
pub const ztrevc_work = LAPACKE_ztrevc_work;

extern fn LAPACKE_strexc_work(layout: c_int, compq: u8, n: isize, t: [*c]f32, ldt: isize, q: [*c]f32, ldq: isize, ifst: [*c]isize, ilst: [*c]isize, work: [*c]f32) isize;
extern fn LAPACKE_dtrexc_work(layout: c_int, compq: u8, n: isize, t: [*c]f64, ldt: isize, q: [*c]f64, ldq: isize, ifst: [*c]isize, ilst: [*c]isize, work: [*c]f64) isize;
extern fn LAPACKE_ctrexc_work(layout: c_int, compq: u8, n: isize, t: [*c]cf32, ldt: isize, q: [*c]cf32, ldq: isize, ifst: isize, ilst: isize) isize;
extern fn LAPACKE_ztrexc_work(layout: c_int, compq: u8, n: isize, t: [*c]cf64, ldt: isize, q: [*c]cf64, ldq: isize, ifst: isize, ilst: isize) isize;
pub const strexc_work = LAPACKE_strexc_work;
pub const dtrexc_work = LAPACKE_dtrexc_work;
pub const ctrexc_work = LAPACKE_ctrexc_work;
pub const ztrexc_work = LAPACKE_ztrexc_work;

extern fn LAPACKE_strrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, x: [*c]const f32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dtrrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, x: [*c]const f64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_ctrrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, x: [*c]const cf32, ldx: isize, ferr: [*c]f32, berr: [*c]f32, work: [*c]cf32, rwork: [*c]f32) isize;
extern fn LAPACKE_ztrrfs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, x: [*c]const cf64, ldx: isize, ferr: [*c]f64, berr: [*c]f64, work: [*c]cf64, rwork: [*c]f64) isize;
pub const strrfs_work = LAPACKE_strrfs_work;
pub const dtrrfs_work = LAPACKE_dtrrfs_work;
pub const ctrrfs_work = LAPACKE_ctrrfs_work;
pub const ztrrfs_work = LAPACKE_ztrrfs_work;

extern fn LAPACKE_strsen_work(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]f32, ldt: isize, q: [*c]f32, ldq: isize, wr: [*c]f32, wi: [*c]f32, m: [*c]isize, s: [*c]f32, sep: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dtrsen_work(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]f64, ldt: isize, q: [*c]f64, ldq: isize, wr: [*c]f64, wi: [*c]f64, m: [*c]isize, s: [*c]f64, sep: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_ctrsen_work(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]cf32, ldt: isize, q: [*c]cf32, ldq: isize, w: [*c]cf32, m: [*c]isize, s: [*c]f32, sep: [*c]f32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_ztrsen_work(layout: c_int, job: u8, compq: u8, select: [*c]const isize, n: isize, t: [*c]cf64, ldt: isize, q: [*c]cf64, ldq: isize, w: [*c]cf64, m: [*c]isize, s: [*c]f64, sep: [*c]f64, work: [*c]cf64, lwork: isize) isize;
pub const strsen_work = LAPACKE_strsen_work;
pub const dtrsen_work = LAPACKE_dtrsen_work;
pub const ctrsen_work = LAPACKE_ctrsen_work;
pub const ztrsen_work = LAPACKE_ztrsen_work;

extern fn LAPACKE_strsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const f32, ldt: isize, vl: [*c]const f32, ldvl: isize, vr: [*c]const f32, ldvr: isize, s: [*c]f32, sep: [*c]f32, mm: isize, m: [*c]isize, work: [*c]f32, ldwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_dtrsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const f64, ldt: isize, vl: [*c]const f64, ldvl: isize, vr: [*c]const f64, ldvr: isize, s: [*c]f64, sep: [*c]f64, mm: isize, m: [*c]isize, work: [*c]f64, ldwork: isize, iwork: [*c]isize) isize;
extern fn LAPACKE_ctrsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const cf32, ldt: isize, vl: [*c]const cf32, ldvl: isize, vr: [*c]const cf32, ldvr: isize, s: [*c]f32, sep: [*c]f32, mm: isize, m: [*c]isize, work: [*c]cf32, ldwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_ztrsna_work(layout: c_int, job: u8, howmny: u8, select: [*c]const isize, n: isize, t: [*c]const cf64, ldt: isize, vl: [*c]const cf64, ldvl: isize, vr: [*c]const cf64, ldvr: isize, s: [*c]f64, sep: [*c]f64, mm: isize, m: [*c]isize, work: [*c]cf64, ldwork: isize, rwork: [*c]f64) isize;
pub const strsna_work = LAPACKE_strsna_work;
pub const dtrsna_work = LAPACKE_dtrsna_work;
pub const ctrsna_work = LAPACKE_ctrsna_work;
pub const ztrsna_work = LAPACKE_ztrsna_work;

extern fn LAPACKE_strsyl_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]f32, ldc: isize, scale: [*c]f32) isize;
extern fn LAPACKE_dtrsyl_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]f64, ldc: isize, scale: [*c]f64) isize;
extern fn LAPACKE_ctrsyl_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize, scale: [*c]f32) isize;
extern fn LAPACKE_ztrsyl_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, scale: [*c]f64) isize;
pub const strsyl_work = LAPACKE_strsyl_work;
pub const dtrsyl_work = LAPACKE_dtrsyl_work;
pub const ctrsyl_work = LAPACKE_ctrsyl_work;
pub const ztrsyl_work = LAPACKE_ztrsyl_work;

extern fn LAPACKE_strsyl3_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f32, lda: isize, b: [*c]const f32, ldb: isize, c: [*c]f32, ldc: isize, scale: [*c]f32, iwork: [*c]isize, liwork: isize, swork: [*c]f32, ldswork: isize) isize;
extern fn LAPACKE_dtrsyl3_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const f64, lda: isize, b: [*c]const f64, ldb: isize, c: [*c]f64, ldc: isize, scale: [*c]f64, iwork: [*c]isize, liwork: isize, swork: [*c]f64, ldswork: isize) isize;
extern fn LAPACKE_ctrsyl3_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf32, lda: isize, b: [*c]const cf32, ldb: isize, c: [*c]cf32, ldc: isize, scale: [*c]f32, swork: [*c]f32, ldswork: isize) isize;
extern fn LAPACKE_ztrsyl3_work(layout: c_int, trana: u8, tranb: u8, isgn: isize, m: isize, n: isize, a: [*c]const cf64, lda: isize, b: [*c]const cf64, ldb: isize, c: [*c]cf64, ldc: isize, scale: [*c]f64, swork: [*c]f64, ldswork: isize) isize;
pub const strsyl3_work = LAPACKE_strsyl3_work;
pub const dtrsyl3_work = LAPACKE_dtrsyl3_work;
pub const ctrsyl3_work = LAPACKE_ctrsyl3_work;
pub const ztrsyl3_work = LAPACKE_ztrsyl3_work;

extern fn LAPACKE_strtri_work(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]f32, lda: isize) isize;
extern fn LAPACKE_dtrtri_work(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]f64, lda: isize) isize;
extern fn LAPACKE_ctrtri_work(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_ztrtri_work(layout: c_int, uplo: u8, diag: u8, n: isize, a: [*c]cf64, lda: isize) isize;
pub const strtri_work = LAPACKE_strtri_work;
pub const dtrtri_work = LAPACKE_dtrtri_work;
pub const ctrtri_work = LAPACKE_ctrtri_work;
pub const ztrtri_work = LAPACKE_ztrtri_work;

extern fn LAPACKE_strtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtrtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctrtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztrtrs_work(layout: c_int, uplo: u8, trans: u8, diag: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const strtrs_work = LAPACKE_strtrs_work;
pub const dtrtrs_work = LAPACKE_dtrtrs_work;
pub const ctrtrs_work = LAPACKE_ctrtrs_work;
pub const ztrtrs_work = LAPACKE_ztrtrs_work;

extern fn LAPACKE_strttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const f32, lda: isize, arf: [*c]f32) isize;
extern fn LAPACKE_dtrttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const f64, lda: isize, arf: [*c]f64) isize;
extern fn LAPACKE_ctrttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, arf: [*c]cf32) isize;
extern fn LAPACKE_ztrttf_work(layout: c_int, transr: u8, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, arf: [*c]cf64) isize;
pub const strttf_work = LAPACKE_strttf_work;
pub const dtrttf_work = LAPACKE_dtrttf_work;
pub const ctrttf_work = LAPACKE_ctrttf_work;
pub const ztrttf_work = LAPACKE_ztrttf_work;

extern fn LAPACKE_strttp_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, ap: [*c]f32) isize;
extern fn LAPACKE_dtrttp_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, ap: [*c]f64) isize;
extern fn LAPACKE_ctrttp_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, ap: [*c]cf32) isize;
extern fn LAPACKE_ztrttp_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, ap: [*c]cf64) isize;
pub const strttp_work = LAPACKE_strttp_work;
pub const dtrttp_work = LAPACKE_dtrttp_work;
pub const ctrttp_work = LAPACKE_ctrttp_work;
pub const ztrttp_work = LAPACKE_ztrttp_work;

extern fn LAPACKE_stzrzf_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, tau: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dtzrzf_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, tau: [*c]f64, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_ctzrzf_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, tau: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_ztzrzf_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, tau: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const stzrzf_work = LAPACKE_stzrzf_work;
pub const dtzrzf_work = LAPACKE_dtzrzf_work;
pub const ctzrzf_work = LAPACKE_ctzrzf_work;
pub const ztzrzf_work = LAPACKE_ztzrzf_work;

extern fn LAPACKE_cungbr_work(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zungbr_work(layout: c_int, vect: u8, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cungbr_work = LAPACKE_cungbr_work;
pub const zungbr_work = LAPACKE_zungbr_work;

extern fn LAPACKE_cunghr_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunghr_work(layout: c_int, n: isize, ilo: isize, ihi: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cunghr_work = LAPACKE_cunghr_work;
pub const zunghr_work = LAPACKE_zunghr_work;

extern fn LAPACKE_cunglq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunglq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cunglq_work = LAPACKE_cunglq_work;
pub const zunglq_work = LAPACKE_zunglq_work;

extern fn LAPACKE_cungql_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zungql_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cungql_work = LAPACKE_cungql_work;
pub const zungql_work = LAPACKE_zungql_work;

extern fn LAPACKE_cungqr_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zungqr_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cungqr_work = LAPACKE_cungqr_work;
pub const zungqr_work = LAPACKE_zungqr_work;

extern fn LAPACKE_cungrq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zungrq_work(layout: c_int, m: isize, n: isize, k: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cungrq_work = LAPACKE_cungrq_work;
pub const zungrq_work = LAPACKE_zungrq_work;

extern fn LAPACKE_cungtr_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, tau: [*c]const cf32, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zungtr_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, tau: [*c]const cf64, work: [*c]cf64, lwork: isize) isize;
pub const cungtr_work = LAPACKE_cungtr_work;
pub const zungtr_work = LAPACKE_zungtr_work;

extern fn LAPACKE_cungtsqr_row_work(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]cf32, lda: isize, t: [*c]const cf32, ldt: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zungtsqr_row_work(layout: c_int, m: isize, n: isize, mb: isize, nb: isize, a: [*c]cf64, lda: isize, t: [*c]const cf64, ldt: isize, work: [*c]cf64, lwork: isize) isize;
pub const cungtsqr_row_work = LAPACKE_cungtsqr_row_work;
pub const zungtsqr_row_work = LAPACKE_zungtsqr_row_work;

extern fn LAPACKE_cunmbr_work(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmbr_work(layout: c_int, vect: u8, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmbr_work = LAPACKE_cunmbr_work;
pub const zunmbr_work = LAPACKE_zunmbr_work;

extern fn LAPACKE_cunmhr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmhr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, ilo: isize, ihi: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmhr_work = LAPACKE_cunmhr_work;
pub const zunmhr_work = LAPACKE_zunmhr_work;

extern fn LAPACKE_cunmlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmlq_work = LAPACKE_cunmlq_work;
pub const zunmlq_work = LAPACKE_zunmlq_work;

extern fn LAPACKE_cunmql_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmql_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmql_work = LAPACKE_cunmql_work;
pub const zunmql_work = LAPACKE_zunmql_work;

extern fn LAPACKE_cunmqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmqr_work = LAPACKE_cunmqr_work;
pub const zunmqr_work = LAPACKE_zunmqr_work;

extern fn LAPACKE_cunmrq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmrq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmrq_work = LAPACKE_cunmrq_work;
pub const zunmrq_work = LAPACKE_zunmrq_work;

extern fn LAPACKE_cunmrz_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmrz_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmrz_work = LAPACKE_cunmrz_work;
pub const zunmrz_work = LAPACKE_zunmrz_work;

extern fn LAPACKE_cunmtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const cf32, lda: isize, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zunmtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, a: [*c]const cf64, lda: isize, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const cunmtr_work = LAPACKE_cunmtr_work;
pub const zunmtr_work = LAPACKE_zunmtr_work;

extern fn LAPACKE_cupgtr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf32, tau: [*c]const cf32, q: [*c]cf32, ldq: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zupgtr_work(layout: c_int, uplo: u8, n: isize, ap: [*c]const cf64, tau: [*c]const cf64, q: [*c]cf64, ldq: isize, work: [*c]cf64) isize;
pub const cupgtr_work = LAPACKE_cupgtr_work;
pub const zupgtr_work = LAPACKE_zupgtr_work;

extern fn LAPACKE_cupmtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const cf32, tau: [*c]const cf32, c: [*c]cf32, ldc: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zupmtr_work(layout: c_int, side: u8, uplo: u8, trans: u8, m: isize, n: isize, ap: [*c]const cf64, tau: [*c]const cf64, c: [*c]cf64, ldc: isize, work: [*c]cf64) isize;
pub const cupmtr_work = LAPACKE_cupmtr_work;
pub const zupmtr_work = LAPACKE_zupmtr_work;

extern fn LAPACKE_claghe(layout: c_int, n: isize, k: isize, d: [*c]const f32, a: [*c]cf32, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_zlaghe(layout: c_int, n: isize, k: isize, d: [*c]const f64, a: [*c]cf64, lda: isize, iseed: [*c]isize) isize;
pub const claghe = LAPACKE_claghe;
pub const zlaghe = LAPACKE_zlaghe;

extern fn LAPACKE_slagsy(layout: c_int, n: isize, k: isize, d: [*c]const f32, a: [*c]f32, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_dlagsy(layout: c_int, n: isize, k: isize, d: [*c]const f64, a: [*c]f64, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_clagsy(layout: c_int, n: isize, k: isize, d: [*c]const f32, a: [*c]cf32, lda: isize, iseed: [*c]isize) isize;
extern fn LAPACKE_zlagsy(layout: c_int, n: isize, k: isize, d: [*c]const f64, a: [*c]cf64, lda: isize, iseed: [*c]isize) isize;
pub const slagsy = LAPACKE_slagsy;
pub const dlagsy = LAPACKE_dlagsy;
pub const clagsy = LAPACKE_clagsy;
pub const zlagsy = LAPACKE_zlagsy;

extern fn LAPACKE_slapmr(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_dlapmr(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f64, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_clapmr(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_zlapmr(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf64, ldx: isize, k: [*c]isize) isize;
pub const slapmr = LAPACKE_slapmr;
pub const dlapmr = LAPACKE_dlapmr;
pub const clapmr = LAPACKE_clapmr;
pub const zlapmr = LAPACKE_zlapmr;

extern fn LAPACKE_slapmt(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_dlapmt(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]f64, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_clapmt(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf32, ldx: isize, k: [*c]isize) isize;
extern fn LAPACKE_zlapmt(layout: c_int, forwrd: isize, m: isize, n: isize, x: [*c]cf64, ldx: isize, k: [*c]isize) isize;
pub const slapmt = LAPACKE_slapmt;
pub const dlapmt = LAPACKE_dlapmt;
pub const clapmt = LAPACKE_clapmt;
pub const zlapmt = LAPACKE_zlapmt;

extern fn LAPACKE_slapy2(x: f32, y: f32) f32;
extern fn LAPACKE_dlapy2(x: f64, y: f64) f64;
pub const slapy2 = LAPACKE_slapy2;
pub const dlapy2 = LAPACKE_dlapy2;

extern fn LAPACKE_slapy3(x: f32, y: f32, z: f32) f32;
extern fn LAPACKE_dlapy3(x: f64, y: f64, z: f64) f64;
pub const slapy3 = LAPACKE_slapy3;
pub const dlapy3 = LAPACKE_dlapy3;

extern fn LAPACKE_slartgp(f: f32, g: f32, cs: [*c]f32, sn: [*c]f32, r: [*c]f32) isize;
extern fn LAPACKE_dlartgp(f: f64, g: f64, cs: [*c]f64, sn: [*c]f64, r: [*c]f64) isize;
pub const slartgp = LAPACKE_slartgp;
pub const dlartgp = LAPACKE_dlartgp;

extern fn LAPACKE_slartgs(x: f32, y: f32, sigma: f32, cs: [*c]f32, sn: [*c]f32) isize;
extern fn LAPACKE_dlartgs(x: f64, y: f64, sigma: f64, cs: [*c]f64, sn: [*c]f64) isize;
pub const slartgs = LAPACKE_slartgs;
pub const dlartgs = LAPACKE_dlartgs;

extern fn LAPACKE_cbbcsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f32, phi: [*c]f32, u1: [*c]cf32, ldu1: isize, u2: [*c]cf32, ldu2: isize, v1t: [*c]cf32, ldv1t: isize, v2t: [*c]cf32, ldv2t: isize, b11d: [*c]f32, b11e: [*c]f32, b12d: [*c]f32, b12e: [*c]f32, b21d: [*c]f32, b21e: [*c]f32, b22d: [*c]f32, b22e: [*c]f32) isize;
pub const cbbcsd = LAPACKE_cbbcsd;

extern fn LAPACKE_cbbcsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f32, phi: [*c]f32, u1: [*c]cf32, ldu1: isize, u2: [*c]cf32, ldu2: isize, v1t: [*c]cf32, ldv1t: isize, v2t: [*c]cf32, ldv2t: isize, b11d: [*c]f32, b11e: [*c]f32, b12d: [*c]f32, b12e: [*c]f32, b21d: [*c]f32, b21e: [*c]f32, b22d: [*c]f32, b22e: [*c]f32, rwork: [*c]f32, lrwork: isize) isize;
pub const cbbcsd_work = LAPACKE_cbbcsd_work;

extern fn LAPACKE_cheswapr(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, i1: isize, i2: isize) isize;
pub const cheswapr = LAPACKE_cheswapr;

extern fn LAPACKE_cheswapr_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, i1: isize, i2: isize) isize;
pub const cheswapr_work = LAPACKE_cheswapr_work;

extern fn LAPACKE_chetri2(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize) isize;
pub const chetri2 = LAPACKE_chetri2;

extern fn LAPACKE_chetri2_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32, lwork: isize) isize;
pub const chetri2_work = LAPACKE_chetri2_work;

extern fn LAPACKE_chetri2x(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, nb: isize) isize;
pub const chetri2x = LAPACKE_chetri2x;

extern fn LAPACKE_chetri2x_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32, nb: isize) isize;
pub const chetri2x_work = LAPACKE_chetri2x_work;

extern fn LAPACKE_chetrs2(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const chetrs2 = LAPACKE_chetrs2;

extern fn LAPACKE_chetrs2_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize, work: [*c]cf32) isize;
pub const chetrs2_work = LAPACKE_chetrs2_work;

extern fn LAPACKE_csyconv(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, e: [*c]cf32) isize;
pub const csyconv = LAPACKE_csyconv;

extern fn LAPACKE_csyconv_work(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, e: [*c]cf32) isize;
pub const csyconv_work = LAPACKE_csyconv_work;

extern fn LAPACKE_csyswapr(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, i1: isize, i2: isize) isize;
pub const csyswapr = LAPACKE_csyswapr;

extern fn LAPACKE_csyswapr_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, i1: isize, i2: isize) isize;
pub const csyswapr_work = LAPACKE_csyswapr_work;

extern fn LAPACKE_csytri2(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize) isize;
pub const csytri2 = LAPACKE_csytri2;

extern fn LAPACKE_csytri2_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32, lwork: isize) isize;
pub const csytri2_work = LAPACKE_csytri2_work;

extern fn LAPACKE_csytri2x(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, nb: isize) isize;
pub const csytri2x = LAPACKE_csytri2x;

extern fn LAPACKE_csytri2x_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]const isize, work: [*c]cf32, nb: isize) isize;
pub const csytri2x_work = LAPACKE_csytri2x_work;

extern fn LAPACKE_csytrs2(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const csytrs2 = LAPACKE_csytrs2;

extern fn LAPACKE_csytrs2_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize, work: [*c]cf32) isize;
pub const csytrs2_work = LAPACKE_csytrs2_work;

extern fn LAPACKE_cunbdb(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf32, ldx11: isize, x12: [*c]cf32, ldx12: isize, x21: [*c]cf32, ldx21: isize, x22: [*c]cf32, ldx22: isize, theta: [*c]f32, phi: [*c]f32, taup1: [*c]cf32, taup2: [*c]cf32, tauq1: [*c]cf32, tauq2: [*c]cf32) isize;
pub const cunbdb = LAPACKE_cunbdb;

extern fn LAPACKE_cunbdb_work(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf32, ldx11: isize, x12: [*c]cf32, ldx12: isize, x21: [*c]cf32, ldx21: isize, x22: [*c]cf32, ldx22: isize, theta: [*c]f32, phi: [*c]f32, taup1: [*c]cf32, taup2: [*c]cf32, tauq1: [*c]cf32, tauq2: [*c]cf32, work: [*c]cf32, lwork: isize) isize;
pub const cunbdb_work = LAPACKE_cunbdb_work;

extern fn LAPACKE_cuncsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf32, ldx11: isize, x12: [*c]cf32, ldx12: isize, x21: [*c]cf32, ldx21: isize, x22: [*c]cf32, ldx22: isize, theta: [*c]f32, u1: [*c]cf32, ldu1: isize, u2: [*c]cf32, ldu2: isize, v1t: [*c]cf32, ldv1t: isize, v2t: [*c]cf32, ldv2t: isize) isize;
pub const cuncsd = LAPACKE_cuncsd;

extern fn LAPACKE_cuncsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf32, ldx11: isize, x12: [*c]cf32, ldx12: isize, x21: [*c]cf32, ldx21: isize, x22: [*c]cf32, ldx22: isize, theta: [*c]f32, u1: [*c]cf32, ldu1: isize, u2: [*c]cf32, ldu2: isize, v1t: [*c]cf32, ldv1t: isize, v2t: [*c]cf32, ldv2t: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize) isize;
pub const cuncsd_work = LAPACKE_cuncsd_work;

extern fn LAPACKE_cuncsd2by1(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]cf32, ldx11: isize, x21: [*c]cf32, ldx21: isize, theta: [*c]f32, u1: [*c]cf32, ldu1: isize, u2: [*c]cf32, ldu2: isize, v1t: [*c]cf32, ldv1t: isize) isize;
pub const cuncsd2by1 = LAPACKE_cuncsd2by1;

extern fn LAPACKE_cuncsd2by1_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]cf32, ldx11: isize, x21: [*c]cf32, ldx21: isize, theta: [*c]f32, u1: [*c]cf32, ldu1: isize, u2: [*c]cf32, ldu2: isize, v1t: [*c]cf32, ldv1t: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize) isize;
pub const cuncsd2by1_work = LAPACKE_cuncsd2by1_work;

extern fn LAPACKE_dbbcsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f64, phi: [*c]f64, u1: [*c]f64, ldu1: isize, u2: [*c]f64, ldu2: isize, v1t: [*c]f64, ldv1t: isize, v2t: [*c]f64, ldv2t: isize, b11d: [*c]f64, b11e: [*c]f64, b12d: [*c]f64, b12e: [*c]f64, b21d: [*c]f64, b21e: [*c]f64, b22d: [*c]f64, b22e: [*c]f64) isize;
pub const dbbcsd = LAPACKE_dbbcsd;

extern fn LAPACKE_dbbcsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f64, phi: [*c]f64, u1: [*c]f64, ldu1: isize, u2: [*c]f64, ldu2: isize, v1t: [*c]f64, ldv1t: isize, v2t: [*c]f64, ldv2t: isize, b11d: [*c]f64, b11e: [*c]f64, b12d: [*c]f64, b12e: [*c]f64, b21d: [*c]f64, b21e: [*c]f64, b22d: [*c]f64, b22e: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const dbbcsd_work = LAPACKE_dbbcsd_work;

extern fn LAPACKE_dorbdb(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f64, ldx11: isize, x12: [*c]f64, ldx12: isize, x21: [*c]f64, ldx21: isize, x22: [*c]f64, ldx22: isize, theta: [*c]f64, phi: [*c]f64, taup1: [*c]f64, taup2: [*c]f64, tauq1: [*c]f64, tauq2: [*c]f64) isize;
pub const dorbdb = LAPACKE_dorbdb;

extern fn LAPACKE_dorbdb_work(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f64, ldx11: isize, x12: [*c]f64, ldx12: isize, x21: [*c]f64, ldx21: isize, x22: [*c]f64, ldx22: isize, theta: [*c]f64, phi: [*c]f64, taup1: [*c]f64, taup2: [*c]f64, tauq1: [*c]f64, tauq2: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const dorbdb_work = LAPACKE_dorbdb_work;

extern fn LAPACKE_dorcsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f64, ldx11: isize, x12: [*c]f64, ldx12: isize, x21: [*c]f64, ldx21: isize, x22: [*c]f64, ldx22: isize, theta: [*c]f64, u1: [*c]f64, ldu1: isize, u2: [*c]f64, ldu2: isize, v1t: [*c]f64, ldv1t: isize, v2t: [*c]f64, ldv2t: isize) isize;
pub const dorcsd = LAPACKE_dorcsd;

extern fn LAPACKE_dorcsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f64, ldx11: isize, x12: [*c]f64, ldx12: isize, x21: [*c]f64, ldx21: isize, x22: [*c]f64, ldx22: isize, theta: [*c]f64, u1: [*c]f64, ldu1: isize, u2: [*c]f64, ldu2: isize, v1t: [*c]f64, ldv1t: isize, v2t: [*c]f64, ldv2t: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
pub const dorcsd_work = LAPACKE_dorcsd_work;

extern fn LAPACKE_dorcsd2by1(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]f64, ldx11: isize, x21: [*c]f64, ldx21: isize, theta: [*c]f64, u1: [*c]f64, ldu1: isize, u2: [*c]f64, ldu2: isize, v1t: [*c]f64, ldv1t: isize) isize;
pub const dorcsd2by1 = LAPACKE_dorcsd2by1;

extern fn LAPACKE_dorcsd2by1_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]f64, ldx11: isize, x21: [*c]f64, ldx21: isize, theta: [*c]f64, u1: [*c]f64, ldu1: isize, u2: [*c]f64, ldu2: isize, v1t: [*c]f64, ldv1t: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize) isize;
pub const dorcsd2by1_work = LAPACKE_dorcsd2by1_work;

extern fn LAPACKE_dsyconv(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, e: [*c]f64) isize;
pub const dsyconv = LAPACKE_dsyconv;

extern fn LAPACKE_dsyconv_work(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, e: [*c]f64) isize;
pub const dsyconv_work = LAPACKE_dsyconv_work;

extern fn LAPACKE_dsyswapr(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, i1: isize, i2: isize) isize;
pub const dsyswapr = LAPACKE_dsyswapr;

extern fn LAPACKE_dsyswapr_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, i1: isize, i2: isize) isize;
pub const dsyswapr_work = LAPACKE_dsyswapr_work;

extern fn LAPACKE_dsytri2(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize) isize;
pub const dsytri2 = LAPACKE_dsytri2;

extern fn LAPACKE_dsytri2_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, work: [*c]f64, lwork: isize) isize;
pub const dsytri2_work = LAPACKE_dsytri2_work;

extern fn LAPACKE_dsytri2x(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, nb: isize) isize;
pub const dsytri2x = LAPACKE_dsytri2x;

extern fn LAPACKE_dsytri2x_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]const isize, work: [*c]f64, nb: isize) isize;
pub const dsytri2x_work = LAPACKE_dsytri2x_work;

extern fn LAPACKE_dsytrs2(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
pub const dsytrs2 = LAPACKE_dsytrs2;

extern fn LAPACKE_dsytrs2_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize, work: [*c]f64) isize;
pub const dsytrs2_work = LAPACKE_dsytrs2_work;

extern fn LAPACKE_sbbcsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f32, phi: [*c]f32, u1: [*c]f32, ldu1: isize, u2: [*c]f32, ldu2: isize, v1t: [*c]f32, ldv1t: isize, v2t: [*c]f32, ldv2t: isize, b11d: [*c]f32, b11e: [*c]f32, b12d: [*c]f32, b12e: [*c]f32, b21d: [*c]f32, b21e: [*c]f32, b22d: [*c]f32, b22e: [*c]f32) isize;
pub const sbbcsd = LAPACKE_sbbcsd;

extern fn LAPACKE_sbbcsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f32, phi: [*c]f32, u1: [*c]f32, ldu1: isize, u2: [*c]f32, ldu2: isize, v1t: [*c]f32, ldv1t: isize, v2t: [*c]f32, ldv2t: isize, b11d: [*c]f32, b11e: [*c]f32, b12d: [*c]f32, b12e: [*c]f32, b21d: [*c]f32, b21e: [*c]f32, b22d: [*c]f32, b22e: [*c]f32, work: [*c]f32, lwork: isize) isize;
pub const sbbcsd_work = LAPACKE_sbbcsd_work;

extern fn LAPACKE_sorbdb(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f32, ldx11: isize, x12: [*c]f32, ldx12: isize, x21: [*c]f32, ldx21: isize, x22: [*c]f32, ldx22: isize, theta: [*c]f32, phi: [*c]f32, taup1: [*c]f32, taup2: [*c]f32, tauq1: [*c]f32, tauq2: [*c]f32) isize;
pub const sorbdb = LAPACKE_sorbdb;

extern fn LAPACKE_sorbdb_work(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f32, ldx11: isize, x12: [*c]f32, ldx12: isize, x21: [*c]f32, ldx21: isize, x22: [*c]f32, ldx22: isize, theta: [*c]f32, phi: [*c]f32, taup1: [*c]f32, taup2: [*c]f32, tauq1: [*c]f32, tauq2: [*c]f32, work: [*c]f32, lwork: isize) isize;
pub const sorbdb_work = LAPACKE_sorbdb_work;

extern fn LAPACKE_sorcsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f32, ldx11: isize, x12: [*c]f32, ldx12: isize, x21: [*c]f32, ldx21: isize, x22: [*c]f32, ldx22: isize, theta: [*c]f32, u1: [*c]f32, ldu1: isize, u2: [*c]f32, ldu2: isize, v1t: [*c]f32, ldv1t: isize, v2t: [*c]f32, ldv2t: isize) isize;
pub const sorcsd = LAPACKE_sorcsd;

extern fn LAPACKE_sorcsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]f32, ldx11: isize, x12: [*c]f32, ldx12: isize, x21: [*c]f32, ldx21: isize, x22: [*c]f32, ldx22: isize, theta: [*c]f32, u1: [*c]f32, ldu1: isize, u2: [*c]f32, ldu2: isize, v1t: [*c]f32, ldv1t: isize, v2t: [*c]f32, ldv2t: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
pub const sorcsd_work = LAPACKE_sorcsd_work;

extern fn LAPACKE_sorcsd2by1(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]f32, ldx11: isize, x21: [*c]f32, ldx21: isize, theta: [*c]f32, u1: [*c]f32, ldu1: isize, u2: [*c]f32, ldu2: isize, v1t: [*c]f32, ldv1t: isize) isize;
pub const sorcsd2by1 = LAPACKE_sorcsd2by1;

extern fn LAPACKE_sorcsd2by1_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]f32, ldx11: isize, x21: [*c]f32, ldx21: isize, theta: [*c]f32, u1: [*c]f32, ldu1: isize, u2: [*c]f32, ldu2: isize, v1t: [*c]f32, ldv1t: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize) isize;
pub const sorcsd2by1_work = LAPACKE_sorcsd2by1_work;

extern fn LAPACKE_ssyconv(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, e: [*c]f32) isize;
pub const ssyconv = LAPACKE_ssyconv;

extern fn LAPACKE_ssyconv_work(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, e: [*c]f32) isize;
pub const ssyconv_work = LAPACKE_ssyconv_work;

extern fn LAPACKE_ssyswapr(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, i1: isize, i2: isize) isize;
pub const ssyswapr = LAPACKE_ssyswapr;

extern fn LAPACKE_ssyswapr_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, i1: isize, i2: isize) isize;
pub const ssyswapr_work = LAPACKE_ssyswapr_work;

extern fn LAPACKE_ssytri2(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize) isize;
pub const ssytri2 = LAPACKE_ssytri2;

extern fn LAPACKE_ssytri2_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, work: [*c]f32, lwork: isize) isize;
pub const ssytri2_work = LAPACKE_ssytri2_work;

extern fn LAPACKE_ssytri2x(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, nb: isize) isize;
pub const ssytri2x = LAPACKE_ssytri2x;

extern fn LAPACKE_ssytri2x_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]const isize, work: [*c]f32, nb: isize) isize;
pub const ssytri2x_work = LAPACKE_ssytri2x_work;

extern fn LAPACKE_ssytrs2(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
pub const ssytrs2 = LAPACKE_ssytrs2;

extern fn LAPACKE_ssytrs2_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize, work: [*c]f32) isize;
pub const ssytrs2_work = LAPACKE_ssytrs2_work;

extern fn LAPACKE_zbbcsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f64, phi: [*c]f64, u1: [*c]cf64, ldu1: isize, u2: [*c]cf64, ldu2: isize, v1t: [*c]cf64, ldv1t: isize, v2t: [*c]cf64, ldv2t: isize, b11d: [*c]f64, b11e: [*c]f64, b12d: [*c]f64, b12e: [*c]f64, b21d: [*c]f64, b21e: [*c]f64, b22d: [*c]f64, b22e: [*c]f64) isize;
pub const zbbcsd = LAPACKE_zbbcsd;

extern fn LAPACKE_zbbcsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, m: isize, p: isize, q: isize, theta: [*c]f64, phi: [*c]f64, u1: [*c]cf64, ldu1: isize, u2: [*c]cf64, ldu2: isize, v1t: [*c]cf64, ldv1t: isize, v2t: [*c]cf64, ldv2t: isize, b11d: [*c]f64, b11e: [*c]f64, b12d: [*c]f64, b12e: [*c]f64, b21d: [*c]f64, b21e: [*c]f64, b22d: [*c]f64, b22e: [*c]f64, rwork: [*c]f64, lrwork: isize) isize;
pub const zbbcsd_work = LAPACKE_zbbcsd_work;

extern fn LAPACKE_zheswapr(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, i1: isize, i2: isize) isize;
pub const zheswapr = LAPACKE_zheswapr;

extern fn LAPACKE_zheswapr_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, i1: isize, i2: isize) isize;
pub const zheswapr_work = LAPACKE_zheswapr_work;

extern fn LAPACKE_zhetri2(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize) isize;
pub const zhetri2 = LAPACKE_zhetri2;

extern fn LAPACKE_zhetri2_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64, lwork: isize) isize;
pub const zhetri2_work = LAPACKE_zhetri2_work;

extern fn LAPACKE_zhetri2x(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, nb: isize) isize;
pub const zhetri2x = LAPACKE_zhetri2x;

extern fn LAPACKE_zhetri2x_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64, nb: isize) isize;
pub const zhetri2x_work = LAPACKE_zhetri2x_work;

extern fn LAPACKE_zhetrs2(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zhetrs2 = LAPACKE_zhetrs2;

extern fn LAPACKE_zhetrs2_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize, work: [*c]cf64) isize;
pub const zhetrs2_work = LAPACKE_zhetrs2_work;

extern fn LAPACKE_zsyconv(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, e: [*c]cf64) isize;
pub const zsyconv = LAPACKE_zsyconv;

extern fn LAPACKE_zsyconv_work(layout: c_int, uplo: u8, way: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, e: [*c]cf64) isize;
pub const zsyconv_work = LAPACKE_zsyconv_work;

extern fn LAPACKE_zsyswapr(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, i1: isize, i2: isize) isize;
pub const zsyswapr = LAPACKE_zsyswapr;

extern fn LAPACKE_zsyswapr_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, i1: isize, i2: isize) isize;
pub const zsyswapr_work = LAPACKE_zsyswapr_work;

extern fn LAPACKE_zsytri2(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize) isize;
pub const zsytri2 = LAPACKE_zsytri2;

extern fn LAPACKE_zsytri2_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64, lwork: isize) isize;
pub const zsytri2_work = LAPACKE_zsytri2_work;

extern fn LAPACKE_zsytri2x(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, nb: isize) isize;
pub const zsytri2x = LAPACKE_zsytri2x;

extern fn LAPACKE_zsytri2x_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]const isize, work: [*c]cf64, nb: isize) isize;
pub const zsytri2x_work = LAPACKE_zsytri2x_work;

extern fn LAPACKE_zsytrs2(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zsytrs2 = LAPACKE_zsytrs2;

extern fn LAPACKE_zsytrs2_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize, work: [*c]cf64) isize;
pub const zsytrs2_work = LAPACKE_zsytrs2_work;

extern fn LAPACKE_zunbdb(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf64, ldx11: isize, x12: [*c]cf64, ldx12: isize, x21: [*c]cf64, ldx21: isize, x22: [*c]cf64, ldx22: isize, theta: [*c]f64, phi: [*c]f64, taup1: [*c]cf64, taup2: [*c]cf64, tauq1: [*c]cf64, tauq2: [*c]cf64) isize;
pub const zunbdb = LAPACKE_zunbdb;

extern fn LAPACKE_zunbdb_work(layout: c_int, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf64, ldx11: isize, x12: [*c]cf64, ldx12: isize, x21: [*c]cf64, ldx21: isize, x22: [*c]cf64, ldx22: isize, theta: [*c]f64, phi: [*c]f64, taup1: [*c]cf64, taup2: [*c]cf64, tauq1: [*c]cf64, tauq2: [*c]cf64, work: [*c]cf64, lwork: isize) isize;
pub const zunbdb_work = LAPACKE_zunbdb_work;

extern fn LAPACKE_zuncsd(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf64, ldx11: isize, x12: [*c]cf64, ldx12: isize, x21: [*c]cf64, ldx21: isize, x22: [*c]cf64, ldx22: isize, theta: [*c]f64, u1: [*c]cf64, ldu1: isize, u2: [*c]cf64, ldu2: isize, v1t: [*c]cf64, ldv1t: isize, v2t: [*c]cf64, ldv2t: isize) isize;
pub const zuncsd = LAPACKE_zuncsd;

extern fn LAPACKE_zuncsd_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, jobv2t: u8, trans: u8, signs: u8, m: isize, p: isize, q: isize, x11: [*c]cf64, ldx11: isize, x12: [*c]cf64, ldx12: isize, x21: [*c]cf64, ldx21: isize, x22: [*c]cf64, ldx22: isize, theta: [*c]f64, u1: [*c]cf64, ldu1: isize, u2: [*c]cf64, ldu2: isize, v1t: [*c]cf64, ldv1t: isize, v2t: [*c]cf64, ldv2t: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize) isize;
pub const zuncsd_work = LAPACKE_zuncsd_work;

extern fn LAPACKE_zuncsd2by1(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]cf64, ldx11: isize, x21: [*c]cf64, ldx21: isize, theta: [*c]f64, u1: [*c]cf64, ldu1: isize, u2: [*c]cf64, ldu2: isize, v1t: [*c]cf64, ldv1t: isize) isize;
pub const zuncsd2by1 = LAPACKE_zuncsd2by1;

extern fn LAPACKE_zuncsd2by1_work(layout: c_int, jobu1: u8, jobu2: u8, jobv1t: u8, m: isize, p: isize, q: isize, x11: [*c]cf64, ldx11: isize, x21: [*c]cf64, ldx21: isize, theta: [*c]f64, u1: [*c]cf64, ldu1: isize, u2: [*c]cf64, ldu2: isize, v1t: [*c]cf64, ldv1t: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize) isize;
pub const zuncsd2by1_work = LAPACKE_zuncsd2by1_work;

extern fn LAPACKE_sgemqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dgemqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, c: [*c]f64, ldc: isize) isize;
extern fn LAPACKE_cgemqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zgemqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, c: [*c]cf64, ldc: isize) isize;
pub const sgemqrt = LAPACKE_sgemqrt;
pub const dgemqrt = LAPACKE_dgemqrt;
pub const cgemqrt = LAPACKE_cgemqrt;
pub const zgemqrt = LAPACKE_zgemqrt;

extern fn LAPACKE_sgeqrt(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dgeqrt(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_cgeqrt(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zgeqrt(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize) isize;
pub const sgeqrt = LAPACKE_sgeqrt;
pub const dgeqrt = LAPACKE_dgeqrt;
pub const cgeqrt = LAPACKE_cgeqrt;
pub const zgeqrt = LAPACKE_zgeqrt;

extern fn LAPACKE_sgeqrt2(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dgeqrt2(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_cgeqrt2(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zgeqrt2(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize) isize;
pub const sgeqrt2 = LAPACKE_sgeqrt2;
pub const dgeqrt2 = LAPACKE_dgeqrt2;
pub const cgeqrt2 = LAPACKE_cgeqrt2;
pub const zgeqrt2 = LAPACKE_zgeqrt2;

extern fn LAPACKE_sgeqrt3(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dgeqrt3(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_cgeqrt3(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zgeqrt3(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize) isize;
pub const sgeqrt3 = LAPACKE_sgeqrt3;
pub const dgeqrt3 = LAPACKE_dgeqrt3;
pub const cgeqrt3 = LAPACKE_cgeqrt3;
pub const zgeqrt3 = LAPACKE_zgeqrt3;

extern fn LAPACKE_stpmqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtpmqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctpmqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztpmqrt(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const stpmqrt = LAPACKE_stpmqrt;
pub const dtpmqrt = LAPACKE_dtpmqrt;
pub const ctpmqrt = LAPACKE_ctpmqrt;
pub const ztpmqrt = LAPACKE_ztpmqrt;

extern fn LAPACKE_stpqrt(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dtpqrt(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_ctpqrt(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_ztpqrt(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, t: [*c]cf64, ldt: isize) isize;
pub const stpqrt = LAPACKE_stpqrt;
pub const dtpqrt = LAPACKE_dtpqrt;
pub const ctpqrt = LAPACKE_ctpqrt;
pub const ztpqrt = LAPACKE_ztpqrt;

extern fn LAPACKE_stpqrt2(layout: c_int, m: isize, n: isize, l: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dtpqrt2(layout: c_int, m: isize, n: isize, l: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_ctpqrt2(layout: c_int, m: isize, n: isize, l: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_ztpqrt2(layout: c_int, m: isize, n: isize, l: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, t: [*c]cf64, ldt: isize) isize;
pub const stpqrt2 = LAPACKE_stpqrt2;
pub const dtpqrt2 = LAPACKE_dtpqrt2;
pub const ctpqrt2 = LAPACKE_ctpqrt2;
pub const ztpqrt2 = LAPACKE_ztpqrt2;

extern fn LAPACKE_stprfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dtprfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_ctprfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_ztprfb(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const stprfb = LAPACKE_stprfb;
pub const dtprfb = LAPACKE_dtprfb;
pub const ctprfb = LAPACKE_ctprfb;
pub const ztprfb = LAPACKE_ztprfb;

extern fn LAPACKE_sgemqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, c: [*c]f32, ldc: isize, work: [*c]f32) isize;
extern fn LAPACKE_dgemqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, c: [*c]f64, ldc: isize, work: [*c]f64) isize;
extern fn LAPACKE_cgemqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, c: [*c]cf32, ldc: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zgemqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, nb: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, c: [*c]cf64, ldc: isize, work: [*c]cf64) isize;
pub const sgemqrt_work = LAPACKE_sgemqrt_work;
pub const dgemqrt_work = LAPACKE_dgemqrt_work;
pub const cgemqrt_work = LAPACKE_cgemqrt_work;
pub const zgemqrt_work = LAPACKE_zgemqrt_work;

extern fn LAPACKE_sgeqrt_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize, work: [*c]f32) isize;
extern fn LAPACKE_dgeqrt_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize, work: [*c]f64) isize;
extern fn LAPACKE_cgeqrt_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize, work: [*c]cf32) isize;
extern fn LAPACKE_zgeqrt_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize, work: [*c]cf64) isize;
pub const sgeqrt_work = LAPACKE_sgeqrt_work;
pub const dgeqrt_work = LAPACKE_dgeqrt_work;
pub const cgeqrt_work = LAPACKE_cgeqrt_work;
pub const zgeqrt_work = LAPACKE_zgeqrt_work;

extern fn LAPACKE_sgeqrt2_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dgeqrt2_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_cgeqrt2_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zgeqrt2_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize) isize;
pub const sgeqrt2_work = LAPACKE_sgeqrt2_work;
pub const dgeqrt2_work = LAPACKE_dgeqrt2_work;
pub const cgeqrt2_work = LAPACKE_cgeqrt2_work;
pub const zgeqrt2_work = LAPACKE_zgeqrt2_work;

extern fn LAPACKE_sgeqrt3_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dgeqrt3_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_cgeqrt3_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zgeqrt3_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize) isize;
pub const sgeqrt3_work = LAPACKE_sgeqrt3_work;
pub const dgeqrt3_work = LAPACKE_dgeqrt3_work;
pub const cgeqrt3_work = LAPACKE_cgeqrt3_work;
pub const zgeqrt3_work = LAPACKE_zgeqrt3_work;

extern fn LAPACKE_stpmqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, work: [*c]f32) isize;
extern fn LAPACKE_dtpmqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, work: [*c]f64) isize;
extern fn LAPACKE_ctpmqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, work: [*c]cf32) isize;
extern fn LAPACKE_ztpmqrt_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, l: isize, nb: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, work: [*c]cf64) isize;
pub const stpmqrt_work = LAPACKE_stpmqrt_work;
pub const dtpmqrt_work = LAPACKE_dtpmqrt_work;
pub const ctpmqrt_work = LAPACKE_ctpmqrt_work;
pub const ztpmqrt_work = LAPACKE_ztpmqrt_work;

extern fn LAPACKE_stpqrt_work(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, t: [*c]f32, ldt: isize, work: [*c]f32) isize;
extern fn LAPACKE_dtpqrt_work(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, t: [*c]f64, ldt: isize, work: [*c]f64) isize;
extern fn LAPACKE_ctpqrt_work(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, t: [*c]cf32, ldt: isize, work: [*c]cf32) isize;
extern fn LAPACKE_ztpqrt_work(layout: c_int, m: isize, n: isize, l: isize, nb: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, t: [*c]cf64, ldt: isize, work: [*c]cf64) isize;
pub const stpqrt_work = LAPACKE_stpqrt_work;
pub const dtpqrt_work = LAPACKE_dtpqrt_work;
pub const ctpqrt_work = LAPACKE_ctpqrt_work;
pub const ztpqrt_work = LAPACKE_ztpqrt_work;

extern fn LAPACKE_stpqrt2_work(layout: c_int, m: isize, n: isize, l: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dtpqrt2_work(layout: c_int, m: isize, n: isize, l: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_ctpqrt2_work(layout: c_int, m: isize, n: isize, l: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_ztpqrt2_work(layout: c_int, m: isize, n: isize, l: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, t: [*c]cf64, ldt: isize) isize;
pub const stpqrt2_work = LAPACKE_stpqrt2_work;
pub const dtpqrt2_work = LAPACKE_dtpqrt2_work;
pub const ctpqrt2_work = LAPACKE_ctpqrt2_work;
pub const ztpqrt2_work = LAPACKE_ztpqrt2_work;

extern fn LAPACKE_stprfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const f32, ldv: isize, t: [*c]const f32, ldt: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, work: [*c]f32, ldwork: isize) isize;
extern fn LAPACKE_dtprfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const f64, ldv: isize, t: [*c]const f64, ldt: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, work: [*c]f64, ldwork: isize) isize;
extern fn LAPACKE_ctprfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const cf32, ldv: isize, t: [*c]const cf32, ldt: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, ldwork: isize) isize;
extern fn LAPACKE_ztprfb_work(layout: c_int, side: u8, trans: u8, direct: u8, storev: u8, m: isize, n: isize, k: isize, l: isize, v: [*c]const cf64, ldv: isize, t: [*c]const cf64, ldt: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, ldwork: isize) isize;
pub const stprfb_work = LAPACKE_stprfb_work;
pub const dtprfb_work = LAPACKE_dtprfb_work;
pub const ctprfb_work = LAPACKE_ctprfb_work;
pub const ztprfb_work = LAPACKE_ztprfb_work;

extern fn LAPACKE_ssysv_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsysv_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csysv_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsysv_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const ssysv_rook = LAPACKE_ssysv_rook;
pub const dsysv_rook = LAPACKE_dsysv_rook;
pub const csysv_rook = LAPACKE_csysv_rook;
pub const zsysv_rook = LAPACKE_zsysv_rook;

extern fn LAPACKE_ssytrf_rook(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dsytrf_rook(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_csytrf_rook(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zsytrf_rook(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const ssytrf_rook = LAPACKE_ssytrf_rook;
pub const dsytrf_rook = LAPACKE_dsytrf_rook;
pub const csytrf_rook = LAPACKE_csytrf_rook;
pub const zsytrf_rook = LAPACKE_zsytrf_rook;

extern fn LAPACKE_ssytrs_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsytrs_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csytrs_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsytrs_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const ssytrs_rook = LAPACKE_ssytrs_rook;
pub const dsytrs_rook = LAPACKE_dsytrs_rook;
pub const csytrs_rook = LAPACKE_csytrs_rook;
pub const zsytrs_rook = LAPACKE_zsytrs_rook;

extern fn LAPACKE_chetrf_rook(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zhetrf_rook(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const chetrf_rook = LAPACKE_chetrf_rook;
pub const zhetrf_rook = LAPACKE_zhetrf_rook;

extern fn LAPACKE_chetrs_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhetrs_rook(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const chetrs_rook = LAPACKE_chetrs_rook;
pub const zhetrs_rook = LAPACKE_zhetrs_rook;

extern fn LAPACKE_csyr(layout: c_int, uplo: u8, n: isize, alpha: cf32, x: [*c]const cf32, incx: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zsyr(layout: c_int, uplo: u8, n: isize, alpha: cf64, x: [*c]const cf64, incx: isize, a: [*c]cf64, lda: isize) isize;
pub const csyr = LAPACKE_csyr;
pub const zsyr = LAPACKE_zsyr;

extern fn LAPACKE_ssysv_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsysv_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csysv_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsysv_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const ssysv_rook_work = LAPACKE_ssysv_rook_work;
pub const dsysv_rook_work = LAPACKE_dsysv_rook_work;
pub const csysv_rook_work = LAPACKE_csysv_rook_work;
pub const zsysv_rook_work = LAPACKE_zsysv_rook_work;

extern fn LAPACKE_ssytrf_rook_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsytrf_rook_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csytrf_rook_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsytrf_rook_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const ssytrf_rook_work = LAPACKE_ssytrf_rook_work;
pub const dsytrf_rook_work = LAPACKE_dsytrf_rook_work;
pub const csytrf_rook_work = LAPACKE_csytrf_rook_work;
pub const zsytrf_rook_work = LAPACKE_zsytrf_rook_work;

extern fn LAPACKE_ssytrs_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dsytrs_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_csytrs_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zsytrs_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const ssytrs_rook_work = LAPACKE_ssytrs_rook_work;
pub const dsytrs_rook_work = LAPACKE_dsytrs_rook_work;
pub const csytrs_rook_work = LAPACKE_csytrs_rook_work;
pub const zsytrs_rook_work = LAPACKE_zsytrs_rook_work;

extern fn LAPACKE_chetrf_rook_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhetrf_rook_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const chetrf_rook_work = LAPACKE_chetrf_rook_work;
pub const zhetrf_rook_work = LAPACKE_zhetrf_rook_work;

extern fn LAPACKE_chetrs_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zhetrs_rook_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const chetrs_rook_work = LAPACKE_chetrs_rook_work;
pub const zhetrs_rook_work = LAPACKE_zhetrs_rook_work;

extern fn LAPACKE_csyr_work(layout: c_int, uplo: u8, n: isize, alpha: cf32, x: [*c]const cf32, incx: isize, a: [*c]cf32, lda: isize) isize;
extern fn LAPACKE_zsyr_work(layout: c_int, uplo: u8, n: isize, alpha: cf64, x: [*c]const cf64, incx: isize, a: [*c]cf64, lda: isize) isize;
pub const csyr_work = LAPACKE_csyr_work;
pub const zsyr_work = LAPACKE_zsyr_work;

extern fn LAPACKE_ilaver(vers_major: [*c]isize, vers_minor: [*c]isize, vers_patch: [*c]isize) void;
pub const ilaver = LAPACKE_ilaver;

extern fn LAPACKE_ssysv_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
pub const ssysv_aa = LAPACKE_ssysv_aa;

extern fn LAPACKE_ssysv_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
pub const ssysv_aa_work = LAPACKE_ssysv_aa_work;

extern fn LAPACKE_dsysv_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
pub const dsysv_aa = LAPACKE_dsysv_aa;

extern fn LAPACKE_dsysv_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
pub const dsysv_aa_work = LAPACKE_dsysv_aa_work;

extern fn LAPACKE_csysv_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const csysv_aa = LAPACKE_csysv_aa;

extern fn LAPACKE_csysv_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const csysv_aa_work = LAPACKE_csysv_aa_work;

extern fn LAPACKE_zsysv_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zsysv_aa = LAPACKE_zsysv_aa;

extern fn LAPACKE_zsysv_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zsysv_aa_work = LAPACKE_zsysv_aa_work;

extern fn LAPACKE_chesv_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const chesv_aa = LAPACKE_chesv_aa;

extern fn LAPACKE_chesv_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const chesv_aa_work = LAPACKE_chesv_aa_work;

extern fn LAPACKE_zhesv_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zhesv_aa = LAPACKE_zhesv_aa;

extern fn LAPACKE_zhesv_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zhesv_aa_work = LAPACKE_zhesv_aa_work;

extern fn LAPACKE_ssytrf_aa(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_dsytrf_aa(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_csytrf_aa(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zsytrf_aa(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const ssytrf_aa = LAPACKE_ssytrf_aa;
pub const dsytrf_aa = LAPACKE_dsytrf_aa;
pub const csytrf_aa = LAPACKE_csytrf_aa;
pub const zsytrf_aa = LAPACKE_zsytrf_aa;

extern fn LAPACKE_chetrf_aa(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize) isize;
extern fn LAPACKE_zhetrf_aa(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize) isize;
pub const chetrf_aa = LAPACKE_chetrf_aa;
pub const zhetrf_aa = LAPACKE_zhetrf_aa;

extern fn LAPACKE_ssytrf_aa_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, ipiv: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsytrf_aa_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, ipiv: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csytrf_aa_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsytrf_aa_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const ssytrf_aa_work = LAPACKE_ssytrf_aa_work;
pub const dsytrf_aa_work = LAPACKE_dsytrf_aa_work;
pub const csytrf_aa_work = LAPACKE_csytrf_aa_work;
pub const zsytrf_aa_work = LAPACKE_zsytrf_aa_work;

extern fn LAPACKE_chetrf_aa_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhetrf_aa_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const chetrf_aa_work = LAPACKE_chetrf_aa_work;
pub const zhetrf_aa_work = LAPACKE_zhetrf_aa_work;

extern fn LAPACKE_csytrs_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const csytrs_aa = LAPACKE_csytrs_aa;

extern fn LAPACKE_csytrs_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const csytrs_aa_work = LAPACKE_csytrs_aa_work;

extern fn LAPACKE_chetrs_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const chetrs_aa = LAPACKE_chetrs_aa;

extern fn LAPACKE_chetrs_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const chetrs_aa_work = LAPACKE_chetrs_aa_work;

extern fn LAPACKE_dsytrs_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
pub const dsytrs_aa = LAPACKE_dsytrs_aa;

extern fn LAPACKE_dsytrs_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, ipiv: [*c]const isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
pub const dsytrs_aa_work = LAPACKE_dsytrs_aa_work;

extern fn LAPACKE_ssytrs_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
pub const ssytrs_aa = LAPACKE_ssytrs_aa;

extern fn LAPACKE_ssytrs_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, ipiv: [*c]const isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
pub const ssytrs_aa_work = LAPACKE_ssytrs_aa_work;

extern fn LAPACKE_zsytrs_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zsytrs_aa = LAPACKE_zsytrs_aa;

extern fn LAPACKE_zsytrs_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zsytrs_aa_work = LAPACKE_zsytrs_aa_work;

extern fn LAPACKE_zhetrs_aa(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zhetrs_aa = LAPACKE_zhetrs_aa;

extern fn LAPACKE_zhetrs_aa_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zhetrs_aa_work = LAPACKE_zhetrs_aa_work;

extern fn LAPACKE_ssysv_rk(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, e: [*c]f32, ipiv: [*c]isize, b: [*c]f32, ldb: isize) isize;
pub const ssysv_rk = LAPACKE_ssysv_rk;

extern fn LAPACKE_ssysv_rk_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, e: [*c]f32, ipiv: [*c]isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
pub const ssysv_rk_work = LAPACKE_ssysv_rk_work;

extern fn LAPACKE_dsysv_rk(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, e: [*c]f64, ipiv: [*c]isize, b: [*c]f64, ldb: isize) isize;
pub const dsysv_rk = LAPACKE_dsysv_rk;

extern fn LAPACKE_dsysv_rk_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, e: [*c]f64, ipiv: [*c]isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
pub const dsysv_rk_work = LAPACKE_dsysv_rk_work;

extern fn LAPACKE_csysv_rk(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const csysv_rk = LAPACKE_csysv_rk;

extern fn LAPACKE_csysv_rk_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const csysv_rk_work = LAPACKE_csysv_rk_work;

extern fn LAPACKE_zsysv_rk(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zsysv_rk = LAPACKE_zsysv_rk;

extern fn LAPACKE_zsysv_rk_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zsysv_rk_work = LAPACKE_zsysv_rk_work;

extern fn LAPACKE_chesv_rk(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const chesv_rk = LAPACKE_chesv_rk;

extern fn LAPACKE_chesv_rk_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const chesv_rk_work = LAPACKE_chesv_rk_work;

extern fn LAPACKE_zhesv_rk(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zhesv_rk = LAPACKE_zhesv_rk;

extern fn LAPACKE_zhesv_rk_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zhesv_rk_work = LAPACKE_zhesv_rk_work;

extern fn LAPACKE_ssytrf_rk(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, e: [*c]f32, ipiv: [*c]isize) isize;
extern fn LAPACKE_dsytrf_rk(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, e: [*c]f64, ipiv: [*c]isize) isize;
extern fn LAPACKE_csytrf_rk(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zsytrf_rk(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize) isize;
pub const ssytrf_rk = LAPACKE_ssytrf_rk;
pub const dsytrf_rk = LAPACKE_dsytrf_rk;
pub const csytrf_rk = LAPACKE_csytrf_rk;
pub const zsytrf_rk = LAPACKE_zsytrf_rk;

extern fn LAPACKE_chetrf_rk(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize) isize;
extern fn LAPACKE_zhetrf_rk(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize) isize;
pub const chetrf_rk = LAPACKE_chetrf_rk;
pub const zhetrf_rk = LAPACKE_zhetrf_rk;

extern fn LAPACKE_ssytrf_rk_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, e: [*c]f32, ipiv: [*c]isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsytrf_rk_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, e: [*c]f64, ipiv: [*c]isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csytrf_rk_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsytrf_rk_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const ssytrf_rk_work = LAPACKE_ssytrf_rk_work;
pub const dsytrf_rk_work = LAPACKE_dsytrf_rk_work;
pub const csytrf_rk_work = LAPACKE_csytrf_rk_work;
pub const zsytrf_rk_work = LAPACKE_zsytrf_rk_work;

extern fn LAPACKE_chetrf_rk_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]cf32, ipiv: [*c]isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhetrf_rk_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]cf64, ipiv: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const chetrf_rk_work = LAPACKE_chetrf_rk_work;
pub const zhetrf_rk_work = LAPACKE_zhetrf_rk_work;

extern fn LAPACKE_csytrs_3(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const csytrs_3 = LAPACKE_csytrs_3;

extern fn LAPACKE_csytrs_3_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const csytrs_3_work = LAPACKE_csytrs_3_work;

extern fn LAPACKE_chetrs_3(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const chetrs_3 = LAPACKE_chetrs_3;

extern fn LAPACKE_chetrs_3_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, b: [*c]cf32, ldb: isize) isize;
pub const chetrs_3_work = LAPACKE_chetrs_3_work;

extern fn LAPACKE_dsytrs_3(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, e: [*c]const f64, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
pub const dsytrs_3 = LAPACKE_dsytrs_3;

extern fn LAPACKE_dsytrs_3_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f64, lda: isize, e: [*c]const f64, ipiv: [*c]const isize, b: [*c]f64, ldb: isize) isize;
pub const dsytrs_3_work = LAPACKE_dsytrs_3_work;

extern fn LAPACKE_ssytrs_3(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, e: [*c]const f32, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
pub const ssytrs_3 = LAPACKE_ssytrs_3;

extern fn LAPACKE_ssytrs_3_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const f32, lda: isize, e: [*c]const f32, ipiv: [*c]const isize, b: [*c]f32, ldb: isize) isize;
pub const ssytrs_3_work = LAPACKE_ssytrs_3_work;

extern fn LAPACKE_zsytrs_3(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zsytrs_3 = LAPACKE_zsytrs_3;

extern fn LAPACKE_zsytrs_3_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zsytrs_3_work = LAPACKE_zsytrs_3_work;

extern fn LAPACKE_zhetrs_3(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zhetrs_3 = LAPACKE_zhetrs_3;

extern fn LAPACKE_zhetrs_3_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, b: [*c]cf64, ldb: isize) isize;
pub const zhetrs_3_work = LAPACKE_zhetrs_3_work;

extern fn LAPACKE_ssytri_3(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, e: [*c]const f32, ipiv: [*c]const isize) isize;
extern fn LAPACKE_dsytri_3(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, e: [*c]const f64, ipiv: [*c]const isize) isize;
extern fn LAPACKE_csytri_3(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zsytri_3(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize) isize;
pub const ssytri_3 = LAPACKE_ssytri_3;
pub const dsytri_3 = LAPACKE_dsytri_3;
pub const csytri_3 = LAPACKE_csytri_3;
pub const zsytri_3 = LAPACKE_zsytri_3;

extern fn LAPACKE_chetri_3(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize) isize;
extern fn LAPACKE_zhetri_3(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize) isize;
pub const chetri_3 = LAPACKE_chetri_3;
pub const zhetri_3 = LAPACKE_zhetri_3;

extern fn LAPACKE_ssytri_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, e: [*c]const f32, ipiv: [*c]const isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsytri_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, e: [*c]const f64, ipiv: [*c]const isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_csytri_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zsytri_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, work: [*c]cf64, lwork: isize) isize;
pub const ssytri_3_work = LAPACKE_ssytri_3_work;
pub const dsytri_3_work = LAPACKE_dsytri_3_work;
pub const csytri_3_work = LAPACKE_csytri_3_work;
pub const zsytri_3_work = LAPACKE_zsytri_3_work;

extern fn LAPACKE_chetri_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zhetri_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, work: [*c]cf64, lwork: isize) isize;
pub const chetri_3_work = LAPACKE_chetri_3_work;
pub const zhetri_3_work = LAPACKE_zhetri_3_work;

extern fn LAPACKE_ssycon_3(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, e: [*c]const f32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_dsycon_3(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, e: [*c]const f64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
extern fn LAPACKE_csycon_3(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zsycon_3(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const ssycon_3 = LAPACKE_ssycon_3;
pub const dsycon_3 = LAPACKE_dsycon_3;
pub const csycon_3 = LAPACKE_csycon_3;
pub const zsycon_3 = LAPACKE_zsycon_3;

extern fn LAPACKE_checon_3(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32) isize;
extern fn LAPACKE_zhecon_3(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64) isize;
pub const checon_3 = LAPACKE_checon_3;
pub const zhecon_3 = LAPACKE_zhecon_3;

extern fn LAPACKE_ssycon_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f32, lda: isize, e: [*c]const f32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]f32, iwork: [*c]isize) isize;
extern fn LAPACKE_dsycon_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]const f64, lda: isize, e: [*c]const f64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]f64, iwork: [*c]isize) isize;
extern fn LAPACKE_csycon_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zsycon_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const ssycon_3_work = LAPACKE_ssycon_3_work;
pub const dsycon_3_work = LAPACKE_dsycon_3_work;
pub const csycon_3_work = LAPACKE_csycon_3_work;
pub const zsycon_3_work = LAPACKE_zsycon_3_work;

extern fn LAPACKE_checon_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf32, lda: isize, e: [*c]const cf32, ipiv: [*c]const isize, anorm: f32, rcond: [*c]f32, work: [*c]cf32) isize;
extern fn LAPACKE_zhecon_3_work(layout: c_int, uplo: u8, n: isize, a: [*c]const cf64, lda: isize, e: [*c]const cf64, ipiv: [*c]const isize, anorm: f64, rcond: [*c]f64, work: [*c]cf64) isize;
pub const checon_3_work = LAPACKE_checon_3_work;
pub const zhecon_3_work = LAPACKE_zhecon_3_work;

extern fn LAPACKE_sgelq(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, tsize: isize) isize;
extern fn LAPACKE_dgelq(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, tsize: isize) isize;
extern fn LAPACKE_cgelq(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, tsize: isize) isize;
extern fn LAPACKE_zgelq(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, tsize: isize) isize;
pub const sgelq = LAPACKE_sgelq;
pub const dgelq = LAPACKE_dgelq;
pub const cgelq = LAPACKE_cgelq;
pub const zgelq = LAPACKE_zgelq;

extern fn LAPACKE_sgelq_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, tsize: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgelq_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, tsize: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgelq_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, tsize: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgelq_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, tsize: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgelq_work = LAPACKE_sgelq_work;
pub const dgelq_work = LAPACKE_dgelq_work;
pub const cgelq_work = LAPACKE_cgelq_work;
pub const zgelq_work = LAPACKE_zgelq_work;

extern fn LAPACKE_sgemlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, t: [*c]const f32, tsize: isize, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dgemlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, t: [*c]const f64, tsize: isize, c: [*c]f64, ldc: isize) isize;
extern fn LAPACKE_cgemlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, t: [*c]const cf32, tsize: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zgemlq(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, t: [*c]const cf64, tsize: isize, c: [*c]cf64, ldc: isize) isize;
pub const sgemlq = LAPACKE_sgemlq;
pub const dgemlq = LAPACKE_dgemlq;
pub const cgemlq = LAPACKE_cgemlq;
pub const zgemlq = LAPACKE_zgemlq;

extern fn LAPACKE_sgemlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, t: [*c]const f32, tsize: isize, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgemlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, t: [*c]const f64, tsize: isize, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgemlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, t: [*c]const cf32, tsize: isize, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgemlq_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, t: [*c]const cf64, tsize: isize, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgemlq_work = LAPACKE_sgemlq_work;
pub const dgemlq_work = LAPACKE_dgemlq_work;
pub const cgemlq_work = LAPACKE_cgemlq_work;
pub const zgemlq_work = LAPACKE_zgemlq_work;

extern fn LAPACKE_sgeqr(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, tsize: isize) isize;
extern fn LAPACKE_dgeqr(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, tsize: isize) isize;
extern fn LAPACKE_cgeqr(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, tsize: isize) isize;
extern fn LAPACKE_zgeqr(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, tsize: isize) isize;
pub const sgeqr = LAPACKE_sgeqr;
pub const dgeqr = LAPACKE_dgeqr;
pub const cgeqr = LAPACKE_cgeqr;
pub const zgeqr = LAPACKE_zgeqr;

extern fn LAPACKE_sgeqr_work(layout: c_int, m: isize, n: isize, a: [*c]f32, lda: isize, t: [*c]f32, tsize: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgeqr_work(layout: c_int, m: isize, n: isize, a: [*c]f64, lda: isize, t: [*c]f64, tsize: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgeqr_work(layout: c_int, m: isize, n: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, tsize: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgeqr_work(layout: c_int, m: isize, n: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, tsize: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgeqr_work = LAPACKE_sgeqr_work;
pub const dgeqr_work = LAPACKE_dgeqr_work;
pub const cgeqr_work = LAPACKE_cgeqr_work;
pub const zgeqr_work = LAPACKE_zgeqr_work;

extern fn LAPACKE_sgemqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, t: [*c]const f32, tsize: isize, c: [*c]f32, ldc: isize) isize;
extern fn LAPACKE_dgemqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, t: [*c]const f64, tsize: isize, c: [*c]f64, ldc: isize) isize;
extern fn LAPACKE_cgemqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, t: [*c]const cf32, tsize: isize, c: [*c]cf32, ldc: isize) isize;
extern fn LAPACKE_zgemqr(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, t: [*c]const cf64, tsize: isize, c: [*c]cf64, ldc: isize) isize;
pub const sgemqr = LAPACKE_sgemqr;
pub const dgemqr = LAPACKE_dgemqr;
pub const cgemqr = LAPACKE_cgemqr;
pub const zgemqr = LAPACKE_zgemqr;

extern fn LAPACKE_sgemqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f32, lda: isize, t: [*c]const f32, tsize: isize, c: [*c]f32, ldc: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgemqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const f64, lda: isize, t: [*c]const f64, tsize: isize, c: [*c]f64, ldc: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgemqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf32, lda: isize, t: [*c]const cf32, tsize: isize, c: [*c]cf32, ldc: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgemqr_work(layout: c_int, side: u8, trans: u8, m: isize, n: isize, k: isize, a: [*c]const cf64, lda: isize, t: [*c]const cf64, tsize: isize, c: [*c]cf64, ldc: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgemqr_work = LAPACKE_sgemqr_work;
pub const dgemqr_work = LAPACKE_dgemqr_work;
pub const cgemqr_work = LAPACKE_cgemqr_work;
pub const zgemqr_work = LAPACKE_zgemqr_work;

extern fn LAPACKE_sgetsls(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize) isize;
extern fn LAPACKE_dgetsls(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize) isize;
extern fn LAPACKE_cgetsls(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize) isize;
extern fn LAPACKE_zgetsls(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize) isize;
pub const sgetsls = LAPACKE_sgetsls;
pub const dgetsls = LAPACKE_dgetsls;
pub const cgetsls = LAPACKE_cgetsls;
pub const zgetsls = LAPACKE_zgetsls;

extern fn LAPACKE_sgetsls_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgetsls_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgetsls_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgetsls_work(layout: c_int, trans: u8, m: isize, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgetsls_work = LAPACKE_sgetsls_work;
pub const dgetsls_work = LAPACKE_dgetsls_work;
pub const cgetsls_work = LAPACKE_cgetsls_work;
pub const zgetsls_work = LAPACKE_zgetsls_work;

extern fn LAPACKE_sgetsqrhrt(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize) isize;
extern fn LAPACKE_dgetsqrhrt(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize) isize;
extern fn LAPACKE_cgetsqrhrt(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize) isize;
extern fn LAPACKE_zgetsqrhrt(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize) isize;
pub const sgetsqrhrt = LAPACKE_sgetsqrhrt;
pub const dgetsqrhrt = LAPACKE_dgetsqrhrt;
pub const cgetsqrhrt = LAPACKE_cgetsqrhrt;
pub const zgetsqrhrt = LAPACKE_zgetsqrhrt;

extern fn LAPACKE_sgetsqrhrt_work(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dgetsqrhrt_work(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize, work: [*c]f64, lwork: isize) isize;
extern fn LAPACKE_cgetsqrhrt_work(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize, work: [*c]cf32, lwork: isize) isize;
extern fn LAPACKE_zgetsqrhrt_work(layout: c_int, m: isize, n: isize, mb1: isize, nb1: isize, nb2: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize, work: [*c]cf64, lwork: isize) isize;
pub const sgetsqrhrt_work = LAPACKE_sgetsqrhrt_work;
pub const dgetsqrhrt_work = LAPACKE_dgetsqrhrt_work;
pub const cgetsqrhrt_work = LAPACKE_cgetsqrhrt_work;
pub const zgetsqrhrt_work = LAPACKE_zgetsqrhrt_work;

extern fn LAPACKE_ssyev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsyev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64) isize;
pub const ssyev_2stage = LAPACKE_ssyev_2stage;
pub const dsyev_2stage = LAPACKE_dsyev_2stage;

extern fn LAPACKE_ssyevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsyevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64) isize;
pub const ssyevd_2stage = LAPACKE_ssyevd_2stage;
pub const dsyevd_2stage = LAPACKE_dsyevd_2stage;

extern fn LAPACKE_ssyevr_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_dsyevr_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize) isize;
pub const ssyevr_2stage = LAPACKE_ssyevr_2stage;
pub const dsyevr_2stage = LAPACKE_dsyevr_2stage;

extern fn LAPACKE_ssyevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsyevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const ssyevx_2stage = LAPACKE_ssyevx_2stage;
pub const dsyevx_2stage = LAPACKE_dsyevx_2stage;

extern fn LAPACKE_ssyev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsyev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const ssyev_2stage_work = LAPACKE_ssyev_2stage_work;
pub const dsyev_2stage_work = LAPACKE_dsyev_2stage_work;

extern fn LAPACKE_ssyevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, w: [*c]f32, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsyevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, w: [*c]f64, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssyevd_2stage_work = LAPACKE_ssyevd_2stage_work;
pub const dsyevd_2stage_work = LAPACKE_dsyevd_2stage_work;

extern fn LAPACKE_ssyevr_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, isuppz: [*c]isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsyevr_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, isuppz: [*c]isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssyevr_2stage_work = LAPACKE_ssyevr_2stage_work;
pub const dsyevr_2stage_work = LAPACKE_dsyevr_2stage_work;

extern fn LAPACKE_ssyevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsyevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const ssyevx_2stage_work = LAPACKE_ssyevx_2stage_work;
pub const dsyevx_2stage_work = LAPACKE_dsyevx_2stage_work;

extern fn LAPACKE_cheev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_zheev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64) isize;
pub const cheev_2stage = LAPACKE_cheev_2stage;
pub const zheev_2stage = LAPACKE_zheev_2stage;

extern fn LAPACKE_cheevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32) isize;
extern fn LAPACKE_zheevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64) isize;
pub const cheevd_2stage = LAPACKE_cheevd_2stage;
pub const zheevd_2stage = LAPACKE_zheevd_2stage;

extern fn LAPACKE_cheevr_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, isuppz: [*c]isize) isize;
extern fn LAPACKE_zheevr_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, isuppz: [*c]isize) isize;
pub const cheevr_2stage = LAPACKE_cheevr_2stage;
pub const zheevr_2stage = LAPACKE_zheevr_2stage;

extern fn LAPACKE_cheevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zheevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const cheevx_2stage = LAPACKE_cheevx_2stage;
pub const zheevx_2stage = LAPACKE_zheevx_2stage;

extern fn LAPACKE_cheev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zheev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const cheev_2stage_work = LAPACKE_cheev_2stage_work;
pub const zheev_2stage_work = LAPACKE_zheev_2stage_work;

extern fn LAPACKE_cheevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zheevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const cheevd_2stage_work = LAPACKE_cheevd_2stage_work;
pub const zheevd_2stage_work = LAPACKE_zheevd_2stage_work;

extern fn LAPACKE_cheevr_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, isuppz: [*c]isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zheevr_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, isuppz: [*c]isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const cheevr_2stage_work = LAPACKE_cheevr_2stage_work;
pub const zheevr_2stage_work = LAPACKE_zheevr_2stage_work;

extern fn LAPACKE_cheevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zheevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const cheevx_2stage_work = LAPACKE_cheevx_2stage_work;
pub const zheevx_2stage_work = LAPACKE_zheevx_2stage_work;

extern fn LAPACKE_ssbev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsbev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const ssbev_2stage = LAPACKE_ssbev_2stage;
pub const dsbev_2stage = LAPACKE_dsbev_2stage;

extern fn LAPACKE_ssbevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize) isize;
extern fn LAPACKE_dsbevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize) isize;
pub const ssbevd_2stage = LAPACKE_ssbevd_2stage;
pub const dsbevd_2stage = LAPACKE_dsbevd_2stage;

extern fn LAPACKE_ssbevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, q: [*c]f32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsbevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, q: [*c]f64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, ifail: [*c]isize) isize;
pub const ssbevx_2stage = LAPACKE_ssbevx_2stage;
pub const dsbevx_2stage = LAPACKE_dsbevx_2stage;

extern fn LAPACKE_ssbev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsbev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize) isize;
pub const ssbev_2stage_work = LAPACKE_ssbev_2stage_work;
pub const dsbev_2stage_work = LAPACKE_dsbev_2stage_work;

extern fn LAPACKE_ssbevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_dsbevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const ssbevd_2stage_work = LAPACKE_ssbevd_2stage_work;
pub const dsbevd_2stage_work = LAPACKE_dsbevd_2stage_work;

extern fn LAPACKE_ssbevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f32, ldab: isize, q: [*c]f32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]f32, ldz: isize, work: [*c]f32, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_dsbevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]f64, ldab: isize, q: [*c]f64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]f64, ldz: isize, work: [*c]f64, lwork: isize, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const ssbevx_2stage_work = LAPACKE_ssbevx_2stage_work;
pub const dsbevx_2stage_work = LAPACKE_dsbevx_2stage_work;

extern fn LAPACKE_chbev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhbev_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chbev_2stage = LAPACKE_chbev_2stage;
pub const zhbev_2stage = LAPACKE_zhbev_2stage;

extern fn LAPACKE_chbevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize) isize;
extern fn LAPACKE_zhbevd_2stage(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize) isize;
pub const chbevd_2stage = LAPACKE_chbevd_2stage;
pub const zhbevd_2stage = LAPACKE_zhbevd_2stage;

extern fn LAPACKE_chbevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, q: [*c]cf32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhbevx_2stage(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, q: [*c]cf64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, ifail: [*c]isize) isize;
pub const chbevx_2stage = LAPACKE_chbevx_2stage;
pub const zhbevx_2stage = LAPACKE_zhbevx_2stage;

extern fn LAPACKE_chbev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zhbev_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const chbev_2stage_work = LAPACKE_chbev_2stage_work;
pub const zhbev_2stage_work = LAPACKE_zhbev_2stage_work;

extern fn LAPACKE_chbevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
extern fn LAPACKE_zhbevd_2stage_work(layout: c_int, jobz: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, lrwork: isize, iwork: [*c]isize, liwork: isize) isize;
pub const chbevd_2stage_work = LAPACKE_chbevd_2stage_work;
pub const zhbevd_2stage_work = LAPACKE_zhbevd_2stage_work;

extern fn LAPACKE_chbevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf32, ldab: isize, q: [*c]cf32, ldq: isize, vl: f32, vu: f32, il: isize, iu: isize, abstol: f32, m: [*c]isize, w: [*c]f32, z: [*c]cf32, ldz: isize, work: [*c]cf32, lwork: isize, rwork: [*c]f32, iwork: [*c]isize, ifail: [*c]isize) isize;
extern fn LAPACKE_zhbevx_2stage_work(layout: c_int, jobz: u8, range: u8, uplo: u8, n: isize, kd: isize, ab: [*c]cf64, ldab: isize, q: [*c]cf64, ldq: isize, vl: f64, vu: f64, il: isize, iu: isize, abstol: f64, m: [*c]isize, w: [*c]f64, z: [*c]cf64, ldz: isize, work: [*c]cf64, lwork: isize, rwork: [*c]f64, iwork: [*c]isize, ifail: [*c]isize) isize;
pub const chbevx_2stage_work = LAPACKE_chbevx_2stage_work;
pub const zhbevx_2stage_work = LAPACKE_zhbevx_2stage_work;

extern fn LAPACKE_ssygv_2stage(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, w: [*c]f32) isize;
extern fn LAPACKE_dsygv_2stage(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, w: [*c]f64) isize;
pub const ssygv_2stage = LAPACKE_ssygv_2stage;
pub const dsygv_2stage = LAPACKE_dsygv_2stage;

extern fn LAPACKE_ssygv_2stage_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f32, lda: isize, b: [*c]f32, ldb: isize, w: [*c]f32, work: [*c]f32, lwork: isize) isize;
extern fn LAPACKE_dsygv_2stage_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]f64, lda: isize, b: [*c]f64, ldb: isize, w: [*c]f64, work: [*c]f64, lwork: isize) isize;
pub const ssygv_2stage_work = LAPACKE_ssygv_2stage_work;
pub const dsygv_2stage_work = LAPACKE_dsygv_2stage_work;

extern fn LAPACKE_chegv_2stage(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, w: [*c]f32) isize;
extern fn LAPACKE_zhegv_2stage(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, w: [*c]f64) isize;
pub const chegv_2stage = LAPACKE_chegv_2stage;
pub const zhegv_2stage = LAPACKE_zhegv_2stage;

extern fn LAPACKE_chegv_2stage_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf32, lda: isize, b: [*c]cf32, ldb: isize, w: [*c]f32, work: [*c]cf32, lwork: isize, rwork: [*c]f32) isize;
extern fn LAPACKE_zhegv_2stage_work(layout: c_int, itype: isize, jobz: u8, uplo: u8, n: isize, a: [*c]cf64, lda: isize, b: [*c]cf64, ldb: isize, w: [*c]f64, work: [*c]cf64, lwork: isize, rwork: [*c]f64) isize;
pub const chegv_2stage_work = LAPACKE_chegv_2stage_work;
pub const zhegv_2stage_work = LAPACKE_zhegv_2stage_work;

extern fn LAPACKE_ssysv_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, tb: [*c]f32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f32, ldb: isize) isize;
pub const ssysv_aa_2stage = LAPACKE_ssysv_aa_2stage;

extern fn LAPACKE_ssysv_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, tb: [*c]f32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f32, ldb: isize, work: [*c]f32, lwork: isize) isize;
pub const ssysv_aa_2stage_work = LAPACKE_ssysv_aa_2stage_work;

extern fn LAPACKE_dsysv_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, tb: [*c]f64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f64, ldb: isize) isize;
pub const dsysv_aa_2stage = LAPACKE_dsysv_aa_2stage;

extern fn LAPACKE_dsysv_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, tb: [*c]f64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f64, ldb: isize, work: [*c]f64, lwork: isize) isize;
pub const dsysv_aa_2stage_work = LAPACKE_dsysv_aa_2stage_work;

extern fn LAPACKE_csysv_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const csysv_aa_2stage = LAPACKE_csysv_aa_2stage;

extern fn LAPACKE_csysv_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const csysv_aa_2stage_work = LAPACKE_csysv_aa_2stage_work;

extern fn LAPACKE_zsysv_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zsysv_aa_2stage = LAPACKE_zsysv_aa_2stage;

extern fn LAPACKE_zsysv_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zsysv_aa_2stage_work = LAPACKE_zsysv_aa_2stage_work;

extern fn LAPACKE_chesv_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const chesv_aa_2stage = LAPACKE_chesv_aa_2stage;

extern fn LAPACKE_chesv_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize, work: [*c]cf32, lwork: isize) isize;
pub const chesv_aa_2stage_work = LAPACKE_chesv_aa_2stage_work;

extern fn LAPACKE_zhesv_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zhesv_aa_2stage = LAPACKE_zhesv_aa_2stage;

extern fn LAPACKE_zhesv_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize, work: [*c]cf64, lwork: isize) isize;
pub const zhesv_aa_2stage_work = LAPACKE_zhesv_aa_2stage_work;

extern fn LAPACKE_ssytrf_aa_2stage(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, tb: [*c]f32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize) isize;
pub const ssytrf_aa_2stage = LAPACKE_ssytrf_aa_2stage;

extern fn LAPACKE_ssytrf_aa_2stage_work(layout: c_int, uplo: u8, n: isize, a: [*c]f32, lda: isize, tb: [*c]f32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, work: [*c]f32, lwork: isize) isize;
pub const ssytrf_aa_2stage_work = LAPACKE_ssytrf_aa_2stage_work;

extern fn LAPACKE_dsytrf_aa_2stage(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, tb: [*c]f64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize) isize;
pub const dsytrf_aa_2stage = LAPACKE_dsytrf_aa_2stage;

extern fn LAPACKE_dsytrf_aa_2stage_work(layout: c_int, uplo: u8, n: isize, a: [*c]f64, lda: isize, tb: [*c]f64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, work: [*c]f64, lwork: isize) isize;
pub const dsytrf_aa_2stage_work = LAPACKE_dsytrf_aa_2stage_work;

extern fn LAPACKE_csytrf_aa_2stage(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize) isize;
pub const csytrf_aa_2stage = LAPACKE_csytrf_aa_2stage;

extern fn LAPACKE_csytrf_aa_2stage_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, work: [*c]cf32, lwork: isize) isize;
pub const csytrf_aa_2stage_work = LAPACKE_csytrf_aa_2stage_work;

extern fn LAPACKE_zsytrf_aa_2stage(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize) isize;
pub const zsytrf_aa_2stage = LAPACKE_zsytrf_aa_2stage;

extern fn LAPACKE_zsytrf_aa_2stage_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const zsytrf_aa_2stage_work = LAPACKE_zsytrf_aa_2stage_work;

extern fn LAPACKE_chetrf_aa_2stage(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize) isize;
pub const chetrf_aa_2stage = LAPACKE_chetrf_aa_2stage;

extern fn LAPACKE_chetrf_aa_2stage_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, work: [*c]cf32, lwork: isize) isize;
pub const chetrf_aa_2stage_work = LAPACKE_chetrf_aa_2stage_work;

extern fn LAPACKE_zhetrf_aa_2stage(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize) isize;
pub const zhetrf_aa_2stage = LAPACKE_zhetrf_aa_2stage;

extern fn LAPACKE_zhetrf_aa_2stage_work(layout: c_int, uplo: u8, n: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, work: [*c]cf64, lwork: isize) isize;
pub const zhetrf_aa_2stage_work = LAPACKE_zhetrf_aa_2stage_work;

extern fn LAPACKE_ssytrs_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, tb: [*c]f32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f32, ldb: isize) isize;
pub const ssytrs_aa_2stage = LAPACKE_ssytrs_aa_2stage;

extern fn LAPACKE_ssytrs_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f32, lda: isize, tb: [*c]f32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f32, ldb: isize) isize;
pub const ssytrs_aa_2stage_work = LAPACKE_ssytrs_aa_2stage_work;

extern fn LAPACKE_dsytrs_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, tb: [*c]f64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f64, ldb: isize) isize;
pub const dsytrs_aa_2stage = LAPACKE_dsytrs_aa_2stage;

extern fn LAPACKE_dsytrs_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]f64, lda: isize, tb: [*c]f64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]f64, ldb: isize) isize;
pub const dsytrs_aa_2stage_work = LAPACKE_dsytrs_aa_2stage_work;

extern fn LAPACKE_csytrs_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const csytrs_aa_2stage = LAPACKE_csytrs_aa_2stage;

extern fn LAPACKE_csytrs_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const csytrs_aa_2stage_work = LAPACKE_csytrs_aa_2stage_work;

extern fn LAPACKE_zsytrs_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zsytrs_aa_2stage = LAPACKE_zsytrs_aa_2stage;

extern fn LAPACKE_zsytrs_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zsytrs_aa_2stage_work = LAPACKE_zsytrs_aa_2stage_work;

extern fn LAPACKE_chetrs_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const chetrs_aa_2stage = LAPACKE_chetrs_aa_2stage;

extern fn LAPACKE_chetrs_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf32, lda: isize, tb: [*c]cf32, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf32, ldb: isize) isize;
pub const chetrs_aa_2stage_work = LAPACKE_chetrs_aa_2stage_work;

extern fn LAPACKE_zhetrs_aa_2stage(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zhetrs_aa_2stage = LAPACKE_zhetrs_aa_2stage;

extern fn LAPACKE_zhetrs_aa_2stage_work(layout: c_int, uplo: u8, n: isize, nrhs: isize, a: [*c]cf64, lda: isize, tb: [*c]cf64, ltb: isize, ipiv: [*c]isize, ipiv2: [*c]isize, b: [*c]cf64, ldb: isize) isize;
pub const zhetrs_aa_2stage_work = LAPACKE_zhetrs_aa_2stage_work;

extern fn LAPACKE_sorhr_col(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize, d: [*c]f32) isize;
pub const sorhr_col = LAPACKE_sorhr_col;

extern fn LAPACKE_sorhr_col_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f32, lda: isize, t: [*c]f32, ldt: isize, d: [*c]f32) isize;
pub const sorhr_col_work = LAPACKE_sorhr_col_work;

extern fn LAPACKE_dorhr_col(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize, d: [*c]f64) isize;
pub const dorhr_col = LAPACKE_dorhr_col;

extern fn LAPACKE_dorhr_col_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]f64, lda: isize, t: [*c]f64, ldt: isize, d: [*c]f64) isize;
pub const dorhr_col_work = LAPACKE_dorhr_col_work;

extern fn LAPACKE_cunhr_col(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize, d: [*c]cf32) isize;
pub const cunhr_col = LAPACKE_cunhr_col;

extern fn LAPACKE_cunhr_col_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf32, lda: isize, t: [*c]cf32, ldt: isize, d: [*c]cf32) isize;
pub const cunhr_col_work = LAPACKE_cunhr_col_work;

extern fn LAPACKE_zunhr_col(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize, d: [*c]cf64) isize;
pub const zunhr_col = LAPACKE_zunhr_col;

extern fn LAPACKE_zunhr_col_work(layout: c_int, m: isize, n: isize, nb: isize, a: [*c]cf64, lda: isize, t: [*c]cf64, ldt: isize, d: [*c]cf64) isize;
pub const zunhr_col_work = LAPACKE_zunhr_col_work;
