const std = @import("std");
const zml = @import("zml");
const cf64 = zml.cf64;
const geru = zml.linalg.blas.geru;

test geru {
    const a = std.testing.allocator;

    const m = 4;
    const n = 5;
    const alpha: cf64 = cf64.init(1, 1);

    const A = try a.alloc(cf64, m * n);
    defer a.free(A);
    const x1 = try a.alloc(cf64, 2 * m);
    defer a.free(x1);
    const y1 = try a.alloc(cf64, 2 * n);
    defer a.free(y1);

    @memcpy(A.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(2, 2),
        cf64.init(3, 3),
        cf64.init(4, 4),
        cf64.init(5, 5),
        cf64.init(6, 6),
        cf64.init(7, 7),
        cf64.init(8, 8),
        cf64.init(9, 9),
        cf64.init(10, 10),
        cf64.init(11, 11),
        cf64.init(12, 12),
        cf64.init(13, 13),
        cf64.init(14, 14),
        cf64.init(15, 15),
        cf64.init(16, 16),
        cf64.init(17, 17),
        cf64.init(18, 18),
        cf64.init(19, 19),
        cf64.init(20, 20),
    });
    @memcpy(x1.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
    });
    @memcpy(y1.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(5, 5),
        cf64.init(0, 0),
    });

    geru(.row_major, m, n, alpha, x1.ptr, 2, y1.ptr, 2, A.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-3, A[0].re);
    try std.testing.expectEqual(3, A[0].im);
    try std.testing.expectEqual(-6, A[1].re);
    try std.testing.expectEqual(6, A[1].im);
    try std.testing.expectEqual(-9, A[2].re);
    try std.testing.expectEqual(9, A[2].im);
    try std.testing.expectEqual(-12, A[3].re);
    try std.testing.expectEqual(12, A[3].im);
    try std.testing.expectEqual(-15, A[4].re);
    try std.testing.expectEqual(15, A[4].im);
    try std.testing.expectEqual(-2, A[5].re);
    try std.testing.expectEqual(12, A[5].im);
    try std.testing.expectEqual(-9, A[6].re);
    try std.testing.expectEqual(19, A[6].im);
    try std.testing.expectEqual(-16, A[7].re);
    try std.testing.expectEqual(26, A[7].im);
    try std.testing.expectEqual(-23, A[8].re);
    try std.testing.expectEqual(33, A[8].im);
    try std.testing.expectEqual(-30, A[9].re);
    try std.testing.expectEqual(40, A[9].im);
    try std.testing.expectEqual(-1, A[10].re);
    try std.testing.expectEqual(21, A[10].im);
    try std.testing.expectEqual(-12, A[11].re);
    try std.testing.expectEqual(32, A[11].im);
    try std.testing.expectEqual(-23, A[12].re);
    try std.testing.expectEqual(43, A[12].im);
    try std.testing.expectEqual(-34, A[13].re);
    try std.testing.expectEqual(54, A[13].im);
    try std.testing.expectEqual(-45, A[14].re);
    try std.testing.expectEqual(65, A[14].im);
    try std.testing.expectEqual(0, A[15].re);
    try std.testing.expectEqual(30, A[15].im);
    try std.testing.expectEqual(-15, A[16].re);
    try std.testing.expectEqual(45, A[16].im);
    try std.testing.expectEqual(-30, A[17].re);
    try std.testing.expectEqual(60, A[17].im);
    try std.testing.expectEqual(-45, A[18].re);
    try std.testing.expectEqual(75, A[18].im);
    try std.testing.expectEqual(-60, A[19].re);
    try std.testing.expectEqual(90, A[19].im);

    const x2 = try a.alloc(cf64, 2 * m);
    defer a.free(x2);
    const y2 = try a.alloc(cf64, 2 * n);
    defer a.free(y2);

    @memcpy(x2.ptr, &[_]cf64{
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y2.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(9, 10),
        cf64.init(0, 0),
    });

    geru(.row_major, m, n, alpha, x2.ptr, -2, y2.ptr, -2, A.ptr, n, .{}) catch unreachable;

    try std.testing.expectEqual(-42, A[0].re);
    try std.testing.expectEqual(20, A[0].im);
    try std.testing.expectEqual(-37, A[1].re);
    try std.testing.expectEqual(19, A[1].im);
    try std.testing.expectEqual(-32, A[2].re);
    try std.testing.expectEqual(18, A[2].im);
    try std.testing.expectEqual(-27, A[3].re);
    try std.testing.expectEqual(17, A[3].im);
    try std.testing.expectEqual(-22, A[4].re);
    try std.testing.expectEqual(16, A[4].im);
    try std.testing.expectEqual(-81, A[5].re);
    try std.testing.expectEqual(65, A[5].im);
    try std.testing.expectEqual(-72, A[6].re);
    try std.testing.expectEqual(60, A[6].im);
    try std.testing.expectEqual(-63, A[7].re);
    try std.testing.expectEqual(55, A[7].im);
    try std.testing.expectEqual(-54, A[8].re);
    try std.testing.expectEqual(50, A[8].im);
    try std.testing.expectEqual(-45, A[9].re);
    try std.testing.expectEqual(45, A[9].im);
    try std.testing.expectEqual(-120, A[10].re);
    try std.testing.expectEqual(110, A[10].im);
    try std.testing.expectEqual(-107, A[11].re);
    try std.testing.expectEqual(101, A[11].im);
    try std.testing.expectEqual(-94, A[12].re);
    try std.testing.expectEqual(92, A[12].im);
    try std.testing.expectEqual(-81, A[13].re);
    try std.testing.expectEqual(83, A[13].im);
    try std.testing.expectEqual(-68, A[14].re);
    try std.testing.expectEqual(74, A[14].im);
    try std.testing.expectEqual(-159, A[15].re);
    try std.testing.expectEqual(155, A[15].im);
    try std.testing.expectEqual(-142, A[16].re);
    try std.testing.expectEqual(142, A[16].im);
    try std.testing.expectEqual(-125, A[17].re);
    try std.testing.expectEqual(129, A[17].im);
    try std.testing.expectEqual(-108, A[18].re);
    try std.testing.expectEqual(116, A[18].im);
    try std.testing.expectEqual(-91, A[19].re);
    try std.testing.expectEqual(103, A[19].im);

    const x3 = try a.alloc(cf64, 2 * m);
    defer a.free(x3);
    const y3 = try a.alloc(cf64, 2 * n);
    defer a.free(y3);

    @memcpy(x3.ptr, &[_]cf64{
        cf64.init(1, 2),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(7, 8),
        cf64.init(0, 0),
    });
    @memcpy(y3.ptr, &[_]cf64{
        cf64.init(1, 1),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(5, 5),
        cf64.init(0, 0),
    });

    geru(.col_major, m, n, alpha, x3.ptr, 2, y3.ptr, 2, A.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-46, A[0].re);
    try std.testing.expectEqual(22, A[0].im);
    try std.testing.expectEqual(-45, A[1].re);
    try std.testing.expectEqual(25, A[1].im);
    try std.testing.expectEqual(-44, A[2].re);
    try std.testing.expectEqual(28, A[2].im);
    try std.testing.expectEqual(-43, A[3].re);
    try std.testing.expectEqual(31, A[3].im);
    try std.testing.expectEqual(-30, A[4].re);
    try std.testing.expectEqual(20, A[4].im);
    try std.testing.expectEqual(-97, A[5].re);
    try std.testing.expectEqual(77, A[5].im);
    try std.testing.expectEqual(-96, A[6].re);
    try std.testing.expectEqual(80, A[6].im);
    try std.testing.expectEqual(-95, A[7].re);
    try std.testing.expectEqual(83, A[7].im);
    try std.testing.expectEqual(-66, A[8].re);
    try std.testing.expectEqual(56, A[8].im);
    try std.testing.expectEqual(-69, A[9].re);
    try std.testing.expectEqual(63, A[9].im);
    try std.testing.expectEqual(-156, A[10].re);
    try std.testing.expectEqual(140, A[10].im);
    try std.testing.expectEqual(-155, A[11].re);
    try std.testing.expectEqual(143, A[11].im);
    try std.testing.expectEqual(-110, A[12].re);
    try std.testing.expectEqual(100, A[12].im);
    try std.testing.expectEqual(-113, A[13].re);
    try std.testing.expectEqual(107, A[13].im);
    try std.testing.expectEqual(-116, A[14].re);
    try std.testing.expectEqual(114, A[14].im);
    try std.testing.expectEqual(-223, A[15].re);
    try std.testing.expectEqual(211, A[15].im);
    try std.testing.expectEqual(-162, A[16].re);
    try std.testing.expectEqual(152, A[16].im);
    try std.testing.expectEqual(-165, A[17].re);
    try std.testing.expectEqual(159, A[17].im);
    try std.testing.expectEqual(-168, A[18].re);
    try std.testing.expectEqual(166, A[18].im);
    try std.testing.expectEqual(-171, A[19].re);
    try std.testing.expectEqual(173, A[19].im);

    const x4 = try a.alloc(cf64, 2 * m);
    defer a.free(x4);
    const y4 = try a.alloc(cf64, 2 * n);
    defer a.free(y4);

    @memcpy(x4.ptr, &[_]cf64{
        cf64.init(7, 8),
        cf64.init(0, 0),
        cf64.init(5, 6),
        cf64.init(0, 0),
        cf64.init(3, 4),
        cf64.init(0, 0),
        cf64.init(1, 2),
        cf64.init(0, 0),
    });
    @memcpy(y4.ptr, &[_]cf64{
        cf64.init(5, 5),
        cf64.init(0, 0),
        cf64.init(4, 4),
        cf64.init(0, 0),
        cf64.init(3, 3),
        cf64.init(0, 0),
        cf64.init(2, 2),
        cf64.init(0, 0),
        cf64.init(1, 1),
        cf64.init(0, 0),
    });

    geru(.col_major, m, n, alpha, x4.ptr, -2, y4.ptr, -2, A.ptr, m, .{}) catch unreachable;

    try std.testing.expectEqual(-50, A[0].re);
    try std.testing.expectEqual(24, A[0].im);
    try std.testing.expectEqual(-53, A[1].re);
    try std.testing.expectEqual(31, A[1].im);
    try std.testing.expectEqual(-56, A[2].re);
    try std.testing.expectEqual(38, A[2].im);
    try std.testing.expectEqual(-59, A[3].re);
    try std.testing.expectEqual(45, A[3].im);
    try std.testing.expectEqual(-38, A[4].re);
    try std.testing.expectEqual(24, A[4].im);
    try std.testing.expectEqual(-113, A[5].re);
    try std.testing.expectEqual(89, A[5].im);
    try std.testing.expectEqual(-120, A[6].re);
    try std.testing.expectEqual(100, A[6].im);
    try std.testing.expectEqual(-127, A[7].re);
    try std.testing.expectEqual(111, A[7].im);
    try std.testing.expectEqual(-78, A[8].re);
    try std.testing.expectEqual(62, A[8].im);
    try std.testing.expectEqual(-93, A[9].re);
    try std.testing.expectEqual(81, A[9].im);
    try std.testing.expectEqual(-192, A[10].re);
    try std.testing.expectEqual(170, A[10].im);
    try std.testing.expectEqual(-203, A[11].re);
    try std.testing.expectEqual(185, A[11].im);
    try std.testing.expectEqual(-126, A[12].re);
    try std.testing.expectEqual(108, A[12].im);
    try std.testing.expectEqual(-145, A[13].re);
    try std.testing.expectEqual(131, A[13].im);
    try std.testing.expectEqual(-164, A[14].re);
    try std.testing.expectEqual(154, A[14].im);
    try std.testing.expectEqual(-287, A[15].re);
    try std.testing.expectEqual(267, A[15].im);
    try std.testing.expectEqual(-182, A[16].re);
    try std.testing.expectEqual(162, A[16].im);
    try std.testing.expectEqual(-205, A[17].re);
    try std.testing.expectEqual(189, A[17].im);
    try std.testing.expectEqual(-228, A[18].re);
    try std.testing.expectEqual(216, A[18].im);
    try std.testing.expectEqual(-251, A[19].re);
    try std.testing.expectEqual(243, A[19].im);
}
