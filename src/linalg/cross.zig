const meta = @import("../meta.zig");
const numeric = @import("../numeric.zig");

const vector = @import("../vector.zig");

const linalg = @import("../linalg.zig");

pub fn Cross(comptime X: type, comptime Y: type) type {
    comptime if (!meta.isStaticVector(X) or !meta.isStaticVector(Y) or
        (X.len != 3 and X.len != 7) or
        (Y.len != 3 and Y.len != 7) or
        X.len != Y.len)
        @compileError("zsl.linalg.Cross: X and Y must be static vectors of length 3 or 7, got\n\tX = " ++
            @typeName(X) ++ "\n\tY = " ++ @typeName(Y) ++ "\n");

    return if (comptime X.len == 3)
        vector.Static(
            X.len,
            numeric.Sub(
                numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
                numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
            ),
        )
    else
        vector.Static(
            X.len,
            numeric.Add(
                numeric.Add(
                    numeric.Sub(
                        numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
                        numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
                    ),
                    numeric.Sub(
                        numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
                        numeric.Mul(meta.Numeric(X), meta.Numeric(Y)),
                    ),
                ),
                numeric.Sub(numeric.Mul(meta.Numeric(X), meta.Numeric(Y)), numeric.Mul(meta.Numeric(X), meta.Numeric(Y))),
            ),
        );
}

/// Computes the cross product of two static vectors, x × y.
///
/// ## Signature
/// ```zig
/// linalg.cross(x: X, y: Y) linalg.Cross(X, Y)
/// ```
///
/// ## Arguments
/// * `x` (`anytype`): The left vector.
/// * `y` (`anytype`): The right vector.
///
/// ## Returns
/// `linalg.Cross(@TypeOf(x), @TypeOf(y))`: The cross product x × y.
pub fn cross(x: anytype, y: anytype) linalg.Cross(@TypeOf(x), @TypeOf(y)) {
    const X = @TypeOf(x);
    const Y = @TypeOf(y);

    var result: linalg.Cross(X, Y) = .init;

    if (comptime X.len == 3) {
        // result[0] = x[1] * y[2] - x[2] * y[1]
        numeric.subInto(
            &result.data[0],
            numeric.mul(x.data[1], y.data[2]),
            numeric.mul(x.data[2], y.data[1]),
        );

        // result[1] = x[2] * y[0] - x[0] * y[2]
        numeric.subInto(
            &result.data[1],
            numeric.mul(x.data[2], y.data[0]),
            numeric.mul(x.data[0], y.data[2]),
        );

        // result[2] = x[0] * y[1] - x[1] * y[0]
        numeric.subInto(
            &result.data[2],
            numeric.mul(x.data[0], y.data[1]),
            numeric.mul(x.data[1], y.data[0]),
        );
    } else { // X.len == 7
        const pair = struct {
            inline fn calc(a1: meta.Numeric(X), b1: meta.Numeric(Y), a2: meta.Numeric(X), b2: meta.Numeric(Y)) numeric.Sub(numeric.Mul(meta.Numeric(X), meta.Numeric(Y)), numeric.Mul(meta.Numeric(X), meta.Numeric(Y))) {
                // a1 * b1 - a2 * b2
                return numeric.sub(
                    numeric.mul(a1, b1),
                    numeric.mul(a2, b2),
                );
            }
        }.calc;

        // result[0] = x[1] * y[2] - x[2] * y[1] + x[3] * y[4] - x[4] * y[3] + x[6] * y[5] - x[5] * y[6]
        numeric.addInto(
            &result.data[0],
            numeric.add(
                pair(x.data[1], y.data[2], x.data[2], y.data[1]),
                pair(x.data[3], y.data[4], x.data[4], y.data[3]),
            ),
            pair(x.data[6], y.data[5], x.data[5], y.data[6]),
        );

        // result[1] = x[2] * y[0] - x[0] * y[2] + x[3] * y[5] - x[5] * y[3] + x[4] * y[6] - x[6] * y[4]
        numeric.addInto(
            &result.data[1],
            numeric.add(
                pair(x.data[2], y.data[0], x.data[0], y.data[2]),
                pair(x.data[3], y.data[5], x.data[5], y.data[3]),
            ),
            pair(x.data[4], y.data[6], x.data[6], y.data[4]),
        );

        // result[2] = x[0] * y[1] - x[1] * y[0] + x[3] * y[6] - x[6] * y[3] + x[5] * y[4] - x[4] * y[5]
        numeric.addInto(
            &result.data[2],
            numeric.add(
                pair(x.data[0], y.data[1], x.data[1], y.data[0]),
                pair(x.data[3], y.data[6], x.data[6], y.data[3]),
            ),
            pair(x.data[5], y.data[4], x.data[4], y.data[5]),
        );

        // result[3] = x[4] * y[0] - x[0] * y[4] + x[5] * y[1] - x[1] * y[5] + x[6] * y[2] - x[2] * y[6]
        numeric.addInto(
            &result.data[3],
            numeric.add(
                pair(x.data[4], y.data[0], x.data[0], y.data[4]),
                pair(x.data[5], y.data[1], x.data[1], y.data[5]),
            ),
            pair(x.data[6], y.data[2], x.data[2], y.data[6]),
        );

        // result[4] = x[0] * y[3] - x[3] * y[0] + x[6] * y[1] - x[1] * y[6] + x[2] * y[5] - x[5] * y[2]
        numeric.addInto(
            &result.data[4],
            numeric.add(
                pair(x.data[0], y.data[3], x.data[3], y.data[0]),
                pair(x.data[6], y.data[1], x.data[1], y.data[6]),
            ),
            pair(x.data[2], y.data[5], x.data[5], y.data[2]),
        );

        // result[5] = x[6] * y[0] - x[0] * y[6] + x[1] * y[3] - x[3] * y[1] + x[4] * y[2] - x[2] * y[4]
        numeric.addInto(
            &result.data[5],
            numeric.add(
                pair(x.data[6], y.data[0], x.data[0], y.data[6]),
                pair(x.data[1], y.data[3], x.data[3], y.data[1]),
            ),
            pair(x.data[4], y.data[2], x.data[2], y.data[4]),
        );

        // result[6] = x[5] * y[0] - x[0] * y[5] + x[1] * y[4] - x[4] * y[1] + x[2] * y[3] - x[3] * y[2]
        numeric.addInto(
            &result.data[6],
            numeric.add(
                pair(x.data[5], y.data[0], x.data[0], y.data[5]),
                pair(x.data[1], y.data[4], x.data[4], y.data[1]),
            ),
            pair(x.data[2], y.data[3], x.data[3], y.data[2]),
        );
    }

    return result;
}
