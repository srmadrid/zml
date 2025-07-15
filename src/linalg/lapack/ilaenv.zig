const std = @import("std");

const types = @import("../../types.zig");
const int = @import("../../int.zig");
const float = @import("../../float.zig");

const lapack = @import("../lapack.zig");

pub fn ilaenv(ispec: isize, comptime name: []const u8, comptime opts: []const u8, n1: isize, n2: isize, n3: isize, n4: isize) isize {
    // Handle computed GOTO equivalent
    if (ispec == 1 or ispec == 2 or ispec == 3) {
        // Continue to label 10
    } else if (ispec == 4) {
        // ispec = 4: number of shifts (used by xHSEQR)
        return 6;
    } else if (ispec == 5) {
        // ispec = 5: minimum column dimension (not used)
        return 2;
    } else if (ispec == 6) {
        // ispec = 6: crossover point for SVD (used by xGELSS and xGESVD)
        return types.scast(isize, types.scast(f32, if (n1 < n2) n1 else n2) * 1.6);
    } else if (ispec == 7) {
        // ispec = 7: number of processors (not used)
        return 1;
    } else if (ispec == 8) {
        // ispec = 8: crossover point for multishift (used by xHSEQR)
        return 50;
    } else if (ispec == 9) {
        // ispec = 9: maximum size of the subproblems at the bottom of the
        //           computation tree in the divide-and-conquer algorithm
        //           (used by xGELSD and xGESDD)
        return 25;
    } else if (ispec == 10) {
        // ispec = 10: ieee and infinity NaN arithmetic can be trusted not to trap
        var ILAENV: isize = 1;
        if (ILAENV == 1) {
            ILAENV = ieeeck(1, 0, 1);
        }
        return ILAENV;
    } else if (ispec == 11) {
        // ispec = 11: ieee infinity arithmetic can be trusted not to trap
        var ILAENV: isize = 1;
        if (ILAENV == 1) {
            ILAENV = ieeeck(0, 0.0, 1.0);
        }
        return ILAENV;
    } else if (ispec >= 12 and ispec <= 17) {
        // 12 <= ispec <= 17: xHSEQR or related subroutines.
        return iparmq(ispec, name, opts, n1, n2, n3, n4);
    } else {
        // Invalid value for ispec
        return -1;
    }

    // Convert name to upper case if the first character is lower case.
    var SUBNAM: [16]u8 = .{0} ** 16;
    @memcpy(SUBNAM[0..name.len], name);

    var IC: isize = types.scast(isize, SUBNAM[0]);
    const IZ: isize = types.scast(isize, 'Z');

    if (IZ == 90 or IZ == 122) {
        // ASCII character set
        if (IC >= 97 and IC <= 122) {
            SUBNAM[0] = types.scast(u8, IC - 32);
            var I: usize = 1;
            while (I < 6 and I < name.len) {
                IC = types.scast(isize, SUBNAM[I]);
                if (IC >= 97 and IC <= 122) {
                    SUBNAM[I] = types.scast(u8, IC - 32);
                }
                I += 1;
            }
        }
    } else if (IZ == 233 or IZ == 169) {
        // EBCDIC character set
        if ((IC >= 129 and IC <= 137) or (IC >= 145 and IC <= 153) or (IC >= 162 and IC <= 169)) {
            SUBNAM[0] = types.scast(u8, IC + 64);
            var I: usize = 1;
            while (I < 6 and I < name.len) {
                IC = types.scast(isize, SUBNAM[I]);
                if ((IC >= 129 and IC <= 137) or (IC >= 145 and IC <= 153) or (IC >= 162 and IC <= 169)) {
                    SUBNAM[I] = types.scast(u8, IC + 64);
                }
                I += 1;
            }
        }
    } else if (IZ == 218 or IZ == 250) {
        // Prime machines: ASCII+128
        if (IC >= 225 and IC <= 250) {
            SUBNAM[0] = types.scast(u8, IC - 32);
            var I: usize = 1;
            while (I < 6 and I < name.len) {
                IC = types.scast(isize, SUBNAM[I]);
                if (IC >= 225 and IC <= 250) {
                    SUBNAM[I] = types.scast(u8, IC - 32);
                }
                I += 1;
            }
        }
    }

    const C1: u8 = SUBNAM[0];
    const Sname: bool = (C1 == 'S' or C1 == 'D');
    const Cname: bool = (C1 == 'C' or C1 == 'Z');

    if (!(Cname or Sname)) {
        return 1;
    }

    var C2: [2]u8 = .{0} ** 2;
    if (name.len >= 2) {
        C2[0] = SUBNAM[1];
        C2[1] = if (name.len >= 3) SUBNAM[2] else 0;
    }

    var C3: [3]u8 = .{0} ** 3;
    var C4: [2]u8 = .{0} ** 2;
    if (name.len >= 6) {
        C3[0] = SUBNAM[3];
        C3[1] = SUBNAM[4];
        C3[2] = SUBNAM[5];
        C4[0] = C3[1];
        C4[1] = C3[2];
    }

    const TWOSTAGE: bool = (name.len >= 11 and SUBNAM[10] == '2');

    if (ispec == 1) {
        // ispec = 1: block size
        // In these examples, separate code is provided for setting NB for
        // real and complex. We assume that NB will take the same value in
        // single or double precision.
        var NB: isize = 1;

        if (std.mem.eql(u8, SUBNAM[1..6], "LAORH")) {
            // This is for *LAORHR_GETRFNP routine
            if (Sname) {
                NB = 32;
            } else {
                NB = 32;
            }
        } else if (std.mem.eql(u8, &C2, "GE")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (Sname) {
                    NB = 64;
                } else {
                    NB = 64;
                }
            } else if (std.mem.eql(u8, &C3, "QRF") or std.mem.eql(u8, &C3, "RQF") or std.mem.eql(u8, &C3, "LQF") or std.mem.eql(u8, &C3, "QLF")) {
                if (Sname) {
                    NB = 32;
                } else {
                    NB = 32;
                }
            } else if (std.mem.eql(u8, &C3, "QR ")) {
                if (n3 == 1) {
                    if (Sname) {
                        // M*N
                        if ((n1 * n2 <= 131072) or (n1 <= 8192)) {
                            NB = n1;
                        } else {
                            NB = int.div(32768, n2);
                        }
                    } else {
                        if ((n1 * n2 <= 131072) or (n1 <= 8192)) {
                            NB = n1;
                        } else {
                            NB = int.div(32768, n2);
                        }
                    }
                } else {
                    if (Sname) {
                        NB = 1;
                    } else {
                        NB = 1;
                    }
                }
            } else if (std.mem.eql(u8, &C3, "LQ ")) {
                if (n3 == 2) {
                    if (Sname) {
                        // M*N
                        if ((n1 * n2 <= 131072) or (n1 <= 8192)) {
                            NB = n1;
                        } else {
                            NB = int.div(3276, n2);
                        }
                    } else {
                        if ((n1 * n2 <= 131072) or (n1 <= 8192)) {
                            NB = n1;
                        } else {
                            NB = int.div(32768, n2);
                        }
                    }
                } else {
                    if (Sname) {
                        NB = 1;
                    } else {
                        NB = 1;
                    }
                }
            } else if (std.mem.eql(u8, &C3, "HRD")) {
                if (Sname) {
                    NB = 32;
                } else {
                    NB = 32;
                }
            } else if (std.mem.eql(u8, &C3, "BRD")) {
                if (Sname) {
                    NB = 32;
                } else {
                    NB = 32;
                }
            } else if (std.mem.eql(u8, &C3, "TRI")) {
                if (Sname) {
                    NB = 64;
                } else {
                    NB = 64;
                }
            } else if (name.len >= 7 and std.mem.eql(u8, SUBNAM[3..8], "QP3RK")) {
                if (Sname) {
                    NB = 32;
                } else {
                    NB = 32;
                }
            }
        } else if (std.mem.eql(u8, &C2, "PO")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (Sname) {
                    NB = 64;
                } else {
                    NB = 64;
                }
            }
        } else if (std.mem.eql(u8, &C2, "SY")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (Sname) {
                    if (TWOSTAGE) {
                        NB = 192;
                    } else {
                        NB = 64;
                    }
                } else {
                    if (TWOSTAGE) {
                        NB = 192;
                    } else {
                        NB = 64;
                    }
                }
            } else if (Sname and std.mem.eql(u8, &C3, "TRD")) {
                NB = 32;
            } else if (Sname and std.mem.eql(u8, &C3, "GST")) {
                NB = 64;
            }
        } else if (Cname and std.mem.eql(u8, &C2, "HE")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (TWOSTAGE) {
                    NB = 192;
                } else {
                    NB = 64;
                }
            } else if (std.mem.eql(u8, &C3, "TRD")) {
                NB = 32;
            } else if (std.mem.eql(u8, &C3, "GST")) {
                NB = 64;
            }
        } else if (Sname and std.mem.eql(u8, &C2, "OR")) {
            if (C3[0] == 'G') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NB = 32;
                }
            } else if (C3[0] == 'M') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NB = 32;
                }
            }
        } else if (Cname and std.mem.eql(u8, &C2, "UN")) {
            if (C3[0] == 'G') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NB = 32;
                }
            } else if (C3[0] == 'M') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NB = 32;
                }
            }
        } else if (std.mem.eql(u8, &C2, "GB")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (Sname) {
                    if (n4 <= 64) {
                        NB = 1;
                    } else {
                        NB = 32;
                    }
                } else {
                    if (n4 <= 64) {
                        NB = 1;
                    } else {
                        NB = 32;
                    }
                }
            }
        } else if (std.mem.eql(u8, &C2, "PB")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (Sname) {
                    if (n2 <= 64) {
                        NB = 1;
                    } else {
                        NB = 32;
                    }
                } else {
                    if (n2 <= 64) {
                        NB = 1;
                    } else {
                        NB = 32;
                    }
                }
            }
        } else if (std.mem.eql(u8, &C2, "TR")) {
            if (std.mem.eql(u8, &C3, "TRI")) {
                if (Sname) {
                    NB = 64;
                } else {
                    NB = 64;
                }
            } else if (std.mem.eql(u8, &C3, "EVC")) {
                if (Sname) {
                    NB = 64;
                } else {
                    NB = 64;
                }
            } else if (std.mem.eql(u8, &C3, "SYL")) {
                // The upper bound is to prevent overly aggressive scaling.
                if (Sname) {
                    const temp: isize = types.scast(isize, types.scast(f32, if (n1 < n2) n1 else n2) * 16 / 100);
                    NB = if ((if (temp > 48) temp else 48) < 240) (if (temp > 48) temp else 48) else 240;
                } else {
                    const temp: isize = types.scast(isize, types.scast(f32, if (n1 < n2) n1 else n2) * 8 / 100);
                    NB = if ((if (temp > 24) temp else 24) < 80) (if (temp > 24) temp else 24) else 80;
                }
            }
        } else if (std.mem.eql(u8, &C2, "LA")) {
            if (std.mem.eql(u8, &C3, "UUM")) {
                if (Sname) {
                    NB = 64;
                } else {
                    NB = 64;
                }
            } else if (std.mem.eql(u8, &C3, "TRS")) {
                if (Sname) {
                    NB = 32;
                } else {
                    NB = 32;
                }
            }
        } else if (Sname and std.mem.eql(u8, &C2, "ST")) {
            if (std.mem.eql(u8, &C3, "EBZ")) {
                NB = 1;
            }
        } else if (std.mem.eql(u8, &C2, "GG")) {
            NB = 32;
            if (std.mem.eql(u8, &C3, "HD3")) {
                if (Sname) {
                    NB = 32;
                } else {
                    NB = 32;
                }
            }
        }
        return NB;
    } else if (ispec == 2) {
        // ispec = 2: minimum block size
        var NBMIN: isize = 2;
        if (std.mem.eql(u8, &C2, "GE")) {
            if (std.mem.eql(u8, &C3, "QRF") or std.mem.eql(u8, &C3, "RQF") or std.mem.eql(u8, &C3, "LQF") or std.mem.eql(u8, &C3, "QLF")) {
                if (Sname) {
                    NBMIN = 2;
                } else {
                    NBMIN = 2;
                }
            } else if (std.mem.eql(u8, &C3, "HRD")) {
                if (Sname) {
                    NBMIN = 2;
                } else {
                    NBMIN = 2;
                }
            } else if (std.mem.eql(u8, &C3, "BRD")) {
                if (Sname) {
                    NBMIN = 2;
                } else {
                    NBMIN = 2;
                }
            } else if (std.mem.eql(u8, &C3, "TRI")) {
                if (Sname) {
                    NBMIN = 2;
                } else {
                    NBMIN = 2;
                }
            } else if (name.len >= 7 and std.mem.eql(u8, SUBNAM[3..8], "QP3RK")) {
                if (Sname) {
                    NBMIN = 2;
                } else {
                    NBMIN = 2;
                }
            }
        } else if (std.mem.eql(u8, &C2, "SY")) {
            if (std.mem.eql(u8, &C3, "TRF")) {
                if (Sname) {
                    NBMIN = 8;
                } else {
                    NBMIN = 8;
                }
            } else if (Sname and std.mem.eql(u8, &C3, "TRD")) {
                NBMIN = 2;
            }
        } else if (Cname and std.mem.eql(u8, &C2, "HE")) {
            if (std.mem.eql(u8, &C3, "TRD")) {
                NBMIN = 2;
            }
        } else if (Sname and std.mem.eql(u8, &C2, "OR")) {
            if (C3[0] == 'G') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NBMIN = 2;
                }
            } else if (C3[0] == 'M') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NBMIN = 2;
                }
            }
        } else if (Cname and std.mem.eql(u8, &C2, "UN")) {
            if (C3[0] == 'G') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NBMIN = 2;
                }
            } else if (C3[0] == 'M') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NBMIN = 2;
                }
            }
        } else if (std.mem.eql(u8, &C2, "GG")) {
            NBMIN = 2;
            if (std.mem.eql(u8, &C3, "HD3")) {
                NBMIN = 2;
            }
        }
        return NBMIN;
    } else if (ispec == 3) {
        // ispec = 3: crossover point
        var NX: isize = 0;
        if (std.mem.eql(u8, &C2, "GE")) {
            if (std.mem.eql(u8, &C3, "QRF") or std.mem.eql(u8, &C3, "RQF") or std.mem.eql(u8, &C3, "LQF") or std.mem.eql(u8, &C3, "QLF")) {
                if (Sname) {
                    NX = 128;
                } else {
                    NX = 128;
                }
            } else if (std.mem.eql(u8, &C3, "HRD")) {
                if (Sname) {
                    NX = 128;
                } else {
                    NX = 128;
                }
            } else if (std.mem.eql(u8, &C3, "BRD")) {
                if (Sname) {
                    NX = 128;
                } else {
                    NX = 128;
                }
            } else if (name.len >= 7 and std.mem.eql(SUBNAM[3..8], "QP3RK")) {
                if (Sname) {
                    NX = 128;
                } else {
                    NX = 128;
                }
            }
        } else if (std.mem.eql(u8, &C2, "SY")) {
            if (Sname and std.mem.eql(u8, &C3, "TRD")) {
                NX = 32;
            }
        } else if (Cname and std.mem.eql(u8, &C2, "HE")) {
            if (std.mem.eql(u8, &C3, "TRD")) {
                NX = 32;
            }
        } else if (Sname and std.mem.eql(u8, &C2, "OR")) {
            if (C3[0] == 'G') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NX = 128;
                }
            }
        } else if (Cname and std.mem.eql(u8, &C2, "UN")) {
            if (C3[0] == 'G') {
                if (std.mem.eql(u8, &C4, "QR") or std.mem.eql(u8, &C4, "RQ") or std.mem.eql(u8, &C4, "LQ") or std.mem.eql(u8, &C4, "QL") or std.mem.eql(u8, &C4, "HR") or std.mem.eql(u8, &C4, "TR") or std.mem.eql(u8, &C4, "BR")) {
                    NX = 128;
                }
            }
        } else if (std.mem.eql(u8, &C2, "GG")) {
            NX = 128;
            if (std.mem.eql(u8, &C3, "HD3")) {
                NX = 128;
            }
        }
        return NX;
    }

    // This should never be reached
    return 1;
}

fn ieeeck(ispec: isize, zero: f32, one: f32) isize {
    var posinf: f32 = one / zero;
    if (posinf <= one) {
        return 0;
    }

    var neginf: f32 = -one / zero;
    if (neginf >= zero) {
        return 0;
    }

    const negzro: f32 = one / (neginf + one);
    if (negzro != zero) {
        return 0;
    }

    neginf = one / negzro;
    if (neginf >= zero) {
        return 0;
    }

    const newzro: f32 = negzro + zero;
    if (newzro != zero) {
        return 0;
    }

    posinf = one / newzro;
    if (posinf <= one) {
        return 0;
    }

    neginf = neginf * posinf;
    if (neginf >= zero) {
        return 0;
    }

    posinf = posinf * posinf;
    if (posinf <= one) {
        return 0;
    }

    //
    //     Return if we were only asked to check infinity arithmetic
    //
    if (ispec != 0) {
        return 1;
    }

    const nan1: f32 = posinf + neginf;
    const nan2: f32 = posinf / neginf;
    const nan3: f32 = posinf / posinf;
    const nan4: f32 = posinf * zero;
    const nan5: f32 = neginf * negzro;
    const nan6: f32 = nan5 * zero;

    if (nan1 == nan1) {
        return 0;
    }

    if (nan2 == nan2) {
        return 0;
    }

    if (nan3 == nan3) {
        return 0;
    }

    if (nan4 == nan4) {
        return 0;
    }

    if (nan5 == nan5) {
        return 0;
    }

    if (nan6 == nan6) {
        return 0;
    }

    return 1;
}

fn iparmq(ispec: isize, comptime name: []const u8, comptime opts: []const u8, n: isize, ilo: isize, ihi: isize, lwork: isize) isize {
    _ = opts; // opts is not used
    _ = n; // n is not used
    _ = lwork; // lwork is not used

    // Parameters
    const INMIN: isize = 12;
    const INWIN: isize = 13;
    const INIBL: isize = 14;
    const ISHFTS: isize = 15;
    const IACC22: isize = 16;
    const ICOST: isize = 17;
    const NMIN: isize = 75;
    const K22MIN: isize = 14;
    const KACMIN: isize = 14;
    const NIBBLE: isize = 14;
    const KNWSWP: isize = 500;
    const RCOST: isize = 10;

    // Initialize
    var iparmq_result: isize = -1;

    var nh: isize = undefined;
    var ns: isize = undefined;
    if ((ispec == ISHFTS) or (ispec == INWIN) or (ispec == IACC22)) {
        // Set the number simultaneous shifts
        nh = ihi - ilo + 1;
        ns = 2;
        if (nh >= 30)
            ns = 4;
        if (nh >= 60)
            ns = 10;
        if (nh >= 150)
            ns = int.max(10, int.div(nh, types.scast(isize, (float.log(nh) / float.log(2)) + 0.5)));
        if (nh >= 590)
            ns = 64;
        if (nh >= 3000)
            ns = 128;
        if (nh >= 6000)
            ns = 256;
        ns = int.max(2, ns - @mod(ns, 2));
    }

    if (ispec == INMIN) {
        // Matrices of order smaller than NMIN get sent
        // to xLAHQR, the classic double shift algorithm.
        // This must be at least 11.
        iparmq_result = NMIN;
    } else if (ispec == INIBL) {
        // INIBL: skip a multi-shift qr iteration and
        // whenever aggressive early deflation finds
        // at least (NIBBLE*(window size)/100) deflations.
        iparmq_result = NIBBLE;
    } else if (ispec == ISHFTS) {
        // NSHFTS: The number of simultaneous shifts
        iparmq_result = ns;
    } else if (ispec == INWIN) {
        // NW: deflation window size.
        if (nh <= KNWSWP) {
            iparmq_result = ns;
        } else {
            iparmq_result = int.div(3 * ns, 2);
        }
    } else if (ispec == IACC22) {
        // IACC22: Whether to accumulate reflections
        // before updating the far-from-diagonal elements
        // and whether to use 2-by-2 block structure while
        // doing it. A small amount of work could be saved
        // by making this choice dependent also upon the
        // NH=IHI-ILO+1.

        // Convert NAME to upper case if the first character is lower case.
        iparmq_result = 0;
        var subnam: [6]u8 = .{0} ** 6;
        @memcpy(subnam[0..name.len], name);

        var ic: isize = types.scast(isize, subnam[0]);
        const iz: isize = types.scast(isize, 'Z');

        if (iz == 90 or iz == 122) {
            // ASCII character set
            if (ic >= 97 and ic <= 122) {
                subnam[0] = types.scast(u8, ic - 32);
                var i: usize = 1;
                while (i < 6) {
                    ic = types.scast(isize, subnam[i]);
                    if (ic >= 97 and ic <= 122)
                        subnam[i] = types.scast(u8, ic - 32);

                    i += 1;
                }
            }
        } else if (iz == 233 or iz == 169) {
            // EBCDIC character set
            if ((ic >= 129 and ic <= 137) or
                (ic >= 145 and ic <= 153) or
                (ic >= 162 and ic <= 169))
            {
                subnam[0] = types.scast(u8, ic + 64);
                var i: usize = 1;
                while (i < 6) {
                    ic = types.scast(isize, subnam[i]);
                    if ((ic >= 129 and ic <= 137) or
                        (ic >= 145 and ic <= 153) or
                        (ic >= 162 and ic <= 169))
                        subnam[i] = types.scast(u8, ic + 64);

                    i += 1;
                }
            }
        } else if (iz == 218 or iz == 250) {
            // Prime machines: ASCII+128
            if (ic >= 225 and ic <= 250) {
                subnam[0] = types.scast(u8, ic - 32);
                var i: usize = 1;
                while (i < 6) {
                    ic = types.scast(isize, subnam[i]);
                    if (ic >= 225 and ic <= 250)
                        subnam[i] = types.scast(u8, ic - 32);

                    i += 1;
                }
            }
        }

        if (std.mem.eql(u8, subnam[1..6], "GGHRD") or
            std.mem.eql(u8, subnam[1..6], "GGHD3"))
        {
            iparmq_result = 1;
            if (nh >= K22MIN)
                iparmq_result = 2;
        } else if (std.mem.eql(u8, subnam[3..6], "EXC")) {
            if (nh >= KACMIN)
                iparmq_result = 1;
            if (nh >= K22MIN)
                iparmq_result = 2;
        } else if (std.mem.eql(u8, subnam[1..6], "HSEQR") or
            std.mem.eql(u8, subnam[1..5], "LAQR"))
        {
            if (ns >= KACMIN)
                iparmq_result = 1;
            if (ns >= K22MIN)
                iparmq_result = 2;
        }
    } else if (ispec == ICOST) {
        // Relative cost of near-the-diagonal chase vs
        // BLAS updates
        iparmq_result = RCOST;
    } else {
        // invalid value of ispec
        iparmq_result = -1;
    }

    return iparmq_result;
}
