const types = @import("../types.zig");
const math = @import("../math.zig");
const atnat = @import("atnat.zig");
const cast = types.cast;
const endian = @import("builtin").cpu.arch.endian();

const t576: [2]u32 = if (endian == .big) .{ 0x63f00000, 0x00000000 } else .{ 0x00000000, 0x63f00000 }; // 2 ^ 576
const tm600: [2]u32 = if (endian == .big) .{ 0x1a700000, 0x00000000 } else .{ 0x00000000, 0x1a700000 }; // 2 ^- 600
const tm24: [2]u32 = if (endian == .big) .{ 0x3e700000, 0x00000000 } else .{ 0x00000000, 0x3e700000 }; // 2 ^- 24
const big: [2]u32 = if (endian == .big) .{ 0x43380000, 0x00000000 } else .{ 0x00000000, 0x43380000 }; //  6755399441055744
const big1: [2]u32 = if (endian == .big) .{ 0x43580000, 0x00000000 } else .{ 0x00000000, 0x43580000 }; // 27021597764222976
const hp0: [2]u32 = if (endian == .big) .{ 0x3ff921fb, 0x54442d18 } else .{ 0x54442d18, 0x3ff921fb }; // 1.5707963267948966
const hp1: [2]u32 = if (endian == .big) .{ 0x3c91a626, 0x33145c07 } else .{ 0x33145c07, 0x3c91a626 }; // 6.123233995736766e-17
const mp1: [2]u32 = if (endian == .big) .{ 0x3ff921fb, 0x58000000 } else .{ 0x58000000, 0x3ff921fb }; // 1.5707963407039642
const mp2: [2]u32 = if (endian == .big) .{ 0xbe4dde97, 0x40000000 } else .{ 0x40000000, 0xbe4dde97 }; //-1.3909067675399456e-08

const toverp: [75]f64 = .{ // 2/PI base 24
    10680707.0, 7228996.0,  1387004.0,  2578385.0,  16069853.0,
    12639074.0, 9804092.0,  4427841.0,  16666979.0, 11263675.0,
    12935607.0, 2387514.0,  4345298.0,  14681673.0, 3074569.0,
    13734428.0, 16653803.0, 1880361.0,  10960616.0, 8533493.0,
    3062596.0,  8710556.0,  7349940.0,  6258241.0,  3772886.0,
    3769171.0,  3798172.0,  8675211.0,  12450088.0, 3874808.0,
    9961438.0,  366607.0,   15675153.0, 9132554.0,  7151469.0,
    3571407.0,  2607881.0,  12013382.0, 4155038.0,  6285869.0,
    7677882.0,  13102053.0, 15825725.0, 473591.0,   9065106.0,
    15363067.0, 6271263.0,  9264392.0,  5636912.0,  4652155.0,
    7056368.0,  13614112.0, 10155062.0, 1944035.0,  9527646.0,
    15080200.0, 6658437.0,  6231200.0,  6832269.0,  16767104.0,
    5075751.0,  3212806.0,  1398474.0,  7579849.0,  6349435.0,
    12618859.0, 4703257.0,  12806093.0, 14477321.0, 2786137.0,
    12875403.0, 9837734.0,  14528324.0, 13719321.0, 343717.0,
};

const split: f64 = 134217729.0; // 2^27 + 1

//******************************************************************/
// Routine  branred() performs range  reduction of a double number */
// x into Double length number a+aa,such that                      */
// x=n*pi/2+(a+aa), abs(a+aa)<pi/4, n=0,+-1,+-2,....               */
// Routine return integer (n mod 4)                                */
//******************************************************************/
pub fn branred(x: f64, a: *f64, aa: *f64) i32 {
    const xx: f64 = x * @as(f64, @bitCast(tm600));
    var t: f64 = xx * split; // split x to two numbers
    const x1: f64 = t - (t - xx);
    const x2: f64 = xx - x1;
    var sum: f64 = 0;
    var u: [2]i32 = @bitCast(x1);
    var k: i32 = (u[atnat.HIGH_HALF] >> 20) & 2047;
    k = @divTrunc(k - 450, 24);
    if (k < 0)
        k = 0;
    var gor: [2]i32 = @bitCast(t576);
    gor[atnat.HIGH_HALF] -= ((k * 24) << 20);
    var r: [6]f64 = undefined;
    var i: i32 = 0;
    while (i < 6) {
        r[@intCast(i)] = x1 * toverp[@intCast(k + i)] * @as(f64, @bitCast(gor));
        gor = @bitCast(@as(f64, @bitCast(gor)) * @as(f64, @bitCast(tm24)));
        i += 1;
    }
    i = 0;
    while (i < 3) {
        const s: f64 = (r[@intCast(i)] + @as(f64, @bitCast(big))) - @as(f64, @bitCast(big));
        sum += s;
        r[@intCast(i)] -= s;
        i += 1;
    }
    t = 0;
    i = 0;
    while (i < 6) {
        t += r[@intCast(5 - i)];
        i += 1;
    }
    var bb: f64 = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
    var s: f64 = (t + @as(f64, @bitCast(big))) - @as(f64, @bitCast(big));
    sum += s;
    t -= s;
    var b: f64 = t + bb;
    bb = (t - b) + bb;
    s = (sum + @as(f64, @bitCast(big1))) - @as(f64, @bitCast(big1));
    sum -= s;
    const b1: f64 = b;
    const bb1: f64 = bb;
    const sum1: f64 = sum;
    sum = 0;

    u = @bitCast(x2);
    k = (u[atnat.HIGH_HALF] >> 20) & 2047;
    k = @divTrunc(k - 450, 24);
    if (k < 0)
        k = 0;
    gor = @bitCast(t576);
    gor[atnat.HIGH_HALF] -= ((k * 24) << 20);
    i = 0;
    while (i < 6) {
        r[@intCast(i)] = x2 * toverp[@intCast(k + i)] * @as(f64, @bitCast(gor));
        gor = @bitCast(@as(f64, @bitCast(gor)) * @as(f64, @bitCast(tm24)));
        i += 1;
    }
    i = 0;
    while (i < 3) {
        s = (r[@intCast(i)] + @as(f64, @bitCast(big))) - @as(f64, @bitCast(big));
        sum += s;
        r[@intCast(i)] -= s;
        i += 1;
    }
    t = 0;
    i = 0;
    while (i < 6) {
        t += r[@intCast(5 - i)];
        i += 1;
    }
    bb = (((((r[0] - t) + r[1]) + r[2]) + r[3]) + r[4]) + r[5];
    s = (t + @as(f64, @bitCast(big))) - @as(f64, @bitCast(big));
    sum += s;
    t -= s;
    b = t + bb;
    bb = (t - b) + bb;
    s = (sum + @as(f64, @bitCast(big1))) - @as(f64, @bitCast(big1));
    sum -= s;

    const b2: f64 = b;
    const bb2: f64 = bb;
    const sum2: f64 = sum;

    sum = sum1 + sum2;
    b = b1 + b2;
    bb = if (math.abs(b1) > math.abs(b2)) (b1 - b) + b2 else (b2 - b) + b1;
    if (b > 0.5) {
        b -= 1.0;
        sum += 1.0;
    } else if (b < -0.5) {
        b += 1.0;
        sum -= 1.0;
    }
    s = b + (bb + bb1 + bb2);
    t = ((b - s) + bb) + (bb1 + bb2);
    b = s * split;
    const t1: f64 = b - (b - s);
    const t2: f64 = s - t1;
    b = s * @as(f64, @bitCast(hp0));
    bb = (((t1 * @as(f64, @bitCast(mp1)) - b) + t1 * @as(f64, @bitCast(mp2))) + t2 * @as(f64, @bitCast(mp1))) + (t2 * @as(f64, @bitCast(mp2)) + s * @as(f64, @bitCast(hp1)) + t * @as(f64, @bitCast(hp0)));
    s = b + bb;
    t = (b - s) + bb;
    a.* = s;
    aa.* = t;
    return (cast(i32, sum, .{})) & 3; // return quarter of unit circle
}
