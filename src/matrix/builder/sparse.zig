const std = @import("std");

const meta = @import("../../meta.zig");
const Layout = meta.Layout;
const Uplo = meta.Uplo;
const Diag = meta.Diag;

const numeric = @import("../../numeric.zig");

const matrix = @import("../../matrix.zig");

/// Sparse builder matrix type, represented in COO format. Three arrays are
/// used to store the row indices, column indices, and values of the non-zero
/// elements. Indices are not sorted and duplicate entries are allowed, getting
/// summed at compilation. This type cannot be used for matrix computations
/// directly; it must first be compiled into a standard sparse matrix.
pub fn Sparse(N: type) type {
    if (!meta.isNumeric(N))
        @compileError("zsl.matrix.builder.Sparse: N must be a numeric type, got \n\tN = " ++ @typeName(N) ++ "\n");

    return struct {
        data: [*]N,
        ridx: [*]usize,
        cidx: [*]usize,
        nnz: usize,
        rows: usize,
        cols: usize,
        _dlen: usize,
        _rlen: usize,
        _clen: usize,
        flags: matrix.Flags,

        // Type signatures
        pub const is_matrix = true;
        pub const is_builder = true;
        pub const is_sparse = true;
        pub const storage_layout = meta.default_layout;
        pub const storage_uplo = meta.default_uplo;
        pub const storage_diag = meta.default_diag;

        // Numeric type
        pub const Numeric = N;

        pub const empty = matrix.builder.Sparse(N){
            .data = &.{},
            .ridx = &.{},
            .cidx = &.{},
            .nnz = 0,
            .rows = 0,
            .cols = 0,
            ._dlen = 0,
            ._rlen = 0,
            ._clen = 0,
            .flags = .{ .owns_data = false },
        };

        /// Initializes a new `matrix.builder.Sparse(N)` with the specified rows
        /// and columns, and an initial capacity for non-zero elements.
        ///
        /// ## Arguments
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        /// * `rows` (`usize`): The rows of the builder matrix.
        /// * `cols` (`usize`): The columns of the builder matrix.
        /// * `nnz` (`usize`): The initial capacity for non-zero elements.
        ///
        /// ## Returns
        /// `matrix.builder.Sparse(N)`: The newly initialized builder matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.ZeroDimension`: If either `rows` or `cols` is zero.
        /// * `matrix.Error.DimensionMismatch`: If `nnz` is zero or greater than
        ///   `rows * cols`.
        pub fn init(allocator: std.mem.Allocator, rows: usize, cols: usize, nnz: usize) !matrix.builder.Sparse(N) {
            if (rows == 0 or cols == 0)
                return matrix.Error.ZeroDimension;

            if (nnz == 0 or nnz > rows * cols)
                return matrix.Error.DimensionMismatch;

            const data: []N = try allocator.alloc(N, nnz);
            errdefer allocator.free(data);

            const ridx: []usize = try allocator.alloc(usize, nnz);
            errdefer allocator.free(ridx);

            return .{
                .data = data.ptr,
                .ridx = ridx.ptr,
                .cidx = (try allocator.alloc(usize, nnz)).ptr,
                .nnz = 0,
                .rows = rows,
                .cols = cols,
                ._dlen = nnz,
                ._rlen = nnz,
                ._clen = nnz,
                .flags = .{ .owns_data = true },
            };
        }

        // pub fn initBuffer

        /// Deinitializes the builder matrix, freeing any allocated memory and
        /// invalidating it.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix to deinitialize.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   deallocation. Must be the same allocator used to initialize
        ///   `self`.
        ///
        /// ## Returns
        /// `void`
        pub fn deinit(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator) void {
            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.ridx[0..self._rlen]);
                allocator.free(self.cidx[0..self._clen]);
            }

            self.* = undefined;
        }

        /// Reserves space for at least `new_nnz` non-zero elements.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix to reserve space for.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `new_nnz` (`usize`): The new capacity for non-zero elements.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn reserve(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, new_nnz: usize) !void {
            if (!self.flags.owns_data)
                return;

            if (new_nnz <= self._dlen and new_nnz <= self._rlen and new_nnz <= self._clen)
                return;

            if (new_nnz > self._dlen) {
                self.data = (try allocator.realloc(self.data[0..self._dlen], new_nnz)).ptr;
                self._dlen = new_nnz;
            }

            if (new_nnz > self._rlen) {
                self.ridx = (try allocator.realloc(self.ridx[0..self._rlen], new_nnz)).ptr;
                self._rlen = new_nnz;
            }

            if (new_nnz > self._clen) {
                self.cidx = (try allocator.realloc(self.cidx[0..self._clen], new_nnz)).ptr;
                self._clen = new_nnz;
            }
        }

        /// Appends a new non-zero element to the builder matrix.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `r` (`usize`): The row index of the element.
        /// * `c` (`usize`): The column index of the element.
        /// * `value` (`N`): The value to insert.
        ///
        /// ## Returns
        /// `void`
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails
        ///   when resizing the internal buffers.
        /// * `matrix.Error.PositionOutOfBounds`: If `r` or `c` is out of bounds.
        /// * `matrix.Error.DataNotOwned`: If the builder matrix does not own
        ///   its data and a resize is required.
        pub fn append(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, r: usize, c: usize, value: N) !void {
            if (r >= self.rows or c >= self.cols)
                return matrix.Error.PositionOutOfBounds;

            if (self.nnz == self._dlen or self.nnz == self._rlen or self.nnz == self._clen) {
                if (!self.flags.owns_data)
                    return matrix.Error.DataNotOwned;

                var new_nnz = self.nnz * 2;
                if (new_nnz == 0)
                    new_nnz = 2;

                try self.reserve(allocator, new_nnz);
            }

            self.data[self.nnz] = value;
            self.ridx[self.nnz] = r;
            self.cidx[self.nnz] = c;
            self.nnz += 1;
        }

        /// Appends a new non-zero element without performing bounds checks
        /// or verifying capacity.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix.
        /// * `r` (`usize`): The row index of the element. Assumed to be within
        ///   bounds.
        /// * `c` (`usize`): The column index of the element. Assumed to be
        ///   within bounds.
        /// * `value` (`N`): The value to insert.
        ///
        /// ## Returns
        /// `void`
        pub fn appendAssumeCapacity(self: *matrix.builder.Sparse(N), r: usize, c: usize, value: N) void {
            self.data[self.nnz] = value;
            self.ridx[self.nnz] = r;
            self.cidx[self.nnz] = c;
            self.nnz += 1;
        }

        // pub fn transpose (swap ridx and cidx)

        /// Creates a copy of the builder matrix.
        ///
        /// ## Arguments
        /// * `self` (`matrix.builder.Sparse(N)`): The builder matrix to copy.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations.
        ///
        /// ## Returns
        /// `matrix.builder.Sparse(N)`: The copied builder matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn copy(self: matrix.builder.Sparse(N), allocator: std.mem.Allocator) !Sparse(N) {
            var data: []N = try allocator.alloc(N, self.nnz);
            errdefer allocator.free(data);
            var ridx: []usize = try allocator.alloc(usize, self.nnz);
            errdefer allocator.free(ridx);
            var cidx: []usize = try allocator.alloc(usize, self.nnz);
            errdefer allocator.free(cidx);

            var i: usize = 0;

            while (i < self.nnz) : (i += 1) {
                data[i] = self.data[i];
                ridx[i] = self.ridx[i];
                cidx[i] = self.cidx[i];
            }

            return .{
                .data = data.ptr,
                .ridx = ridx.ptr,
                .cidx = cidx.ptr,
                .nnz = self.nnz,
                .rows = self.rows,
                .cols = self.cols,
                ._dlen = self.nnz,
                ._rlen = self.nnz,
                ._clen = self.nnz,
                .flags = .{ .owns_data = true },
            };
        }

        fn compileInternal(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime M: type) !M {
            const perm = try allocator.alloc(usize, self.nnz);
            defer allocator.free(perm);

            var i: usize = 0;
            while (i < perm.len) : (i += 1) {
                perm[i] = i;
            }

            const Context = struct {
                r: [*]usize,
                c: [*]usize,

                pub fn lessThan(ctx: @This(), a: usize, b: usize) bool {
                    if (comptime meta.layoutOf(M) == .col_major) {
                        if (ctx.c[a] != ctx.c[b])
                            return ctx.c[a] < ctx.c[b];

                        return ctx.r[a] < ctx.r[b];
                    } else {
                        if (ctx.r[a] != ctx.r[b])
                            return ctx.r[a] < ctx.r[b];

                        return ctx.c[a] < ctx.c[b];
                    }
                }
            };

            std.mem.sortUnstable(usize, perm, Context{ .r = self.ridx, .c = self.cidx }, Context.lessThan);

            var unique_nnz: usize = 0;
            if (self.nnz > 0) {
                unique_nnz = 1;
                var last_r = self.ridx[perm[0]];
                var last_c = self.cidx[perm[0]];

                i = 1;
                while (i < self.nnz) : (i += 1) {
                    const r = self.ridx[perm[i]];
                    const c = self.cidx[perm[i]];
                    if (r != last_r or c != last_c) {
                        unique_nnz += 1;
                        last_r = r;
                        last_c = c;
                    }
                }
            }

            var ptr = try allocator.alloc(usize, if (comptime meta.layoutOf(M) == .col_major) self.cols + 1 else self.rows + 1);
            errdefer allocator.free(ptr);
            @memset(ptr, 0);

            var data = try allocator.alloc(N, unique_nnz);
            errdefer allocator.free(data);

            var idx = try allocator.alloc(usize, unique_nnz);
            errdefer allocator.free(idx);

            if (self.nnz > 0) {
                var write_idx: usize = 0;
                var current_r = self.ridx[perm[0]];
                var current_c = self.cidx[perm[0]];
                var current_val = self.data[perm[0]];

                i = 1;
                while (i < self.nnz) : (i += 1) {
                    const p = perm[i];
                    const r = self.ridx[p];
                    const c = self.cidx[p];
                    const val = self.data[p];

                    if (r == current_r and c == current_c) {
                        numeric.addInto(&current_val, current_val, val);
                    } else {
                        data[write_idx] = current_val;
                        idx[write_idx] = if (comptime meta.layoutOf(M) == .col_major) current_r else current_c;
                        const major_idx = if (comptime meta.layoutOf(M) == .col_major) current_c else current_r;
                        ptr[major_idx + 1] += 1;

                        write_idx += 1;
                        current_r = r;
                        current_c = c;
                        current_val = val;
                    }
                }

                data[write_idx] = current_val;
                idx[write_idx] = if (comptime meta.layoutOf(M) == .col_major) current_r else current_c;
                const major_idx = if (comptime meta.layoutOf(M) == .col_major) current_c else current_r;
                ptr[major_idx + 1] += 1;
            }

            i = 0;
            while (i < ptr.len - 1) : (i += 1) {
                ptr[i + 1] += ptr[i];
            }

            if (self.flags.owns_data) {
                allocator.free(self.data[0..self._dlen]);
                allocator.free(self.ridx[0..self._rlen]);
                allocator.free(self.cidx[0..self._clen]);
            }

            var result: M = undefined;
            result.data = data.ptr;
            result._dlen = unique_nnz;
            result.idx = idx.ptr;
            result._ilen = unique_nnz;
            result.ptr = ptr.ptr;
            result.nnz = unique_nnz;
            result.flags = .{ .owns_data = true };

            result.rows = self.rows;
            result.cols = self.cols;

            self.* = undefined;

            return result;
        }

        /// Compiles the builder matrix into a general sparse matrix,
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix to compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.general.Sparse(N, layout)`: The compiled general sparse
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn compile(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime layout: Layout) !matrix.general.Sparse(N, layout) {
            return self.compileInternal(allocator, matrix.general.Sparse(N, layout));
        }

        /// Compiles the builder matrix into a general sparse matrix by copying
        /// the data, leaving the builder intact.
        ///
        /// ## Arguments
        /// * `self` (`matrix.builder.Sparse(N)`): The builder matrix to
        ///   compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.general.Sparse(N)`: The compiled general sparse matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn compileCopy(self: matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime layout: Layout) !matrix.general.Sparse(N, layout) {
            var cpy = try self.copy(allocator);
            errdefer cpy.deinit(allocator);
            return cpy.compileInternal(allocator, matrix.general.Sparse(N, layout));
        }

        fn removeTriangle(self: *matrix.builder.Sparse(N), comptime uplo: Uplo, comptime diagonal: bool) void {
            var i: usize = 0;
            var j: usize = 0;
            while (i < self.nnz) : (i += 1) {
                const r = self.ridx[i];
                const c = self.cidx[i];

                const keep = switch (comptime uplo) {
                    .upper => if (diagonal) r > c else r >= c,
                    .lower => if (diagonal) r < c else r <= c,
                };

                if (keep) {
                    if (i != j) {
                        self.data[j] = self.data[i];
                        self.ridx[j] = self.ridx[i];
                        self.cidx[j] = self.cidx[i];
                    }

                    j += 1;
                }
            }

            self.nnz = j;
        }

        /// Compiles the builder matrix into a symmetric sparse matrix, keeping
        /// only the specified triangle part and discarding the other, and
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix to compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `uplo` (`comptime Uplo`): Specifies which triangle part to keep.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.symmetric.Sparse(N, uplo)`: The compiled symmetric sparse
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.NotSquare`: If the builder matrix is not square.
        pub fn compileSymmetric(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime layout: Layout) !matrix.symmetric.Sparse(N, uplo, layout) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            self.removeTriangle(comptime uplo.invert(), false);
            return self.compileInternal(allocator, matrix.symmetric.Sparse(N, uplo, layout));
        }

        /// Compiles the builder matrix into a symmetric sparse matrix by
        /// copying the data, leaving the builder intact, keeping only the
        /// specified triangle part and discarding the other.
        ///
        /// ## Arguments
        /// * `self` (`matrix.builder.Sparse(N)`): The builder matrix to
        ///   compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `uplo` (`comptime Uplo`): Specifies which triangle part to keep.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.symmetric.Sparse(N, uplo)`: The compiled symmetric sparse
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.NotSquare`: If the builder matrix is not square.
        pub fn compileSymmetricCopy(self: matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime layout: Layout) !matrix.symmetric.Sparse(N, uplo, layout) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            var cpy = try self.copy(allocator);
            errdefer cpy.deinit(allocator);
            cpy.removeTriangle(comptime uplo.invert(), false);
            return cpy.compileInternal(allocator, matrix.symmetric.Sparse(N, uplo, layout));
        }

        /// Compiles the builder matrix into a Hermitian sparse matrix, keeping
        /// only the specified triangle part and discarding the other, and
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix to compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `uplo` (`comptime Uplo`): Specifies which triangle part to keep.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.hermitian.Sparse(N, uplo)`: The compiled Hermitian sparse
        /// matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.NotSquare`: If the builder matrix is not square.
        pub fn compileHermitian(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime layout: Layout) !matrix.hermitian.Sparse(N, uplo, layout) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            self.removeTriangle(comptime uplo.invert(), false);
            return self.compileInternal(allocator, matrix.hermitian.Sparse(N, uplo, layout));
        }

        /// Compiles the builder matrix into a hermitian sparse matrix by
        /// copying the data, leaving the builder intact, keeping only the
        /// specified triangle part and discarding the other.
        ///
        /// ## Arguments
        /// * `self` (`matrix.builder.Sparse(N)`): The builder matrix to
        ///   compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `uplo` (`comptime Uplo`): Specifies which triangle part to keep.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.hermitian.Sparse(N, uplo)`: The compiled hermitian
        /// sparse matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        /// * `matrix.Error.NotSquare`: If the builder matrix is not square.
        pub fn compileHermitianCopy(self: matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime layout: Layout) !matrix.hermitian.Sparse(N, uplo, layout) {
            if (self.rows != self.cols)
                return matrix.Error.NotSquare;

            var cpy = try self.copy(allocator);
            errdefer cpy.deinit(allocator);
            cpy.removeTriangle(comptime uplo.invert(), false);
            return cpy.compileInternal(allocator, matrix.hermitian.Sparse(N, uplo, layout));
        }

        /// Compiles the builder matrix into a triangular sparse matrix, keeping
        /// only the specified triangle part and discarding the other, and
        /// transferring ownership of the data and invalidating the builder.
        ///
        /// ## Arguments
        /// * `self` (`*matrix.builder.Sparse(N)`): A pointer to the builder
        ///   matrix to compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `uplo` (`comptime Uplo`): Specifies which triangle part to keep.
        /// * `diag` (`comptime Diag`): Specifies whether the diagonal is unit
        ///   or non-unit.
        /// * `layout` (`comptime Layout`): Specifies the memory layout of the
        ///   compiled matrix: `.col_major` for CSC and `.row_major` for CSR.
        ///
        /// ## Returns
        /// `matrix.triangular.Sparse(N, uplo, diag, laayout)`: The compiled
        /// triangular sparse matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn compileTriangular(self: *matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime diag: Diag, comptime layout: Layout) !matrix.triangular.Sparse(N, uplo, diag, layout) {
            self.removeTriangle(comptime uplo.invert(), comptime diag == .unit);
            return self.compileInternal(allocator, matrix.triangular.Sparse(N, uplo, diag, layout));
        }

        /// Compiles the builder matrix into a triangular sparse matrix by
        /// copying the data, leaving the builder intact, keeping only the
        /// specified triangle part and discarding the other.
        ///
        /// ## Arguments
        /// * `self` (`matrix.builder.Sparse(N)`): The builder matrix to
        ///   compile.
        /// * `allocator` (`std.mem.Allocator`): The allocator to use for memory
        ///   allocations. Must be the same allocator used to initialize `self`.
        /// * `uplo` (`comptime Uplo`): Specifies which triangle part to keep.
        /// * `diag` (`comptime Diag`): Specifies whether the diagonal is unit
        ///   or non-unit.
        ///
        /// ## Returns
        /// `matrix.triangular.Sparse(N, uplo, diag, layout)`: The compiled
        /// triangular sparse matrix.
        ///
        /// ## Errors
        /// * `std.mem.Allocator.Error.OutOfMemory`: If memory allocation fails.
        pub fn compileTriangularCopy(self: matrix.builder.Sparse(N), allocator: std.mem.Allocator, comptime uplo: Uplo, comptime diag: Diag, comptime layout: Layout) !matrix.triangular.Sparse(N, uplo, diag, layout) {
            var cpy = try self.copy(allocator);
            errdefer cpy.deinit(allocator);
            cpy.removeTriangle(comptime uplo.invert(), comptime diag == .unit);
            return cpy.compileInternal(allocator, matrix.triangular.Sparse(N, uplo, diag, layout));
        }
    };
}
