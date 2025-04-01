/// Returns a value with the magnitude of `x` and the sign of `y`.
pub fn copysign(x: anytype, y: anytype) @TypeOf(x) {
    return if ((x < 0 and y > 0) or (x > 0 and y < 0)) -x else x;
}
