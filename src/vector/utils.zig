pub fn hasNZ(v: anytype, row: usize, col: usize) bool {
    if (row > 0 and col > 0) return false;

    const target = if (row > 0) row else col;

    var lo: usize = 0;
    var hi: usize = v.nnz;
    while (lo < hi) {
        const mid = lo + (hi - lo) / 2;
        if (v.idx[mid] == target) return true;
        if (v.idx[mid] < target) lo = mid + 1 else hi = mid;
    }
    return false;
}
