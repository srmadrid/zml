const std = @import("std");
const options = @import("options");

const types = @import("../types.zig");
const Order = types.Order;
const Numeric = types.Numeric;
const Coerce = types.Coerce;
const orderOf = types.orderOf;
const int = @import("../int.zig");

const vector = @import("../vector.zig");
const matrix = @import("../matrix.zig");

const linalg = @import("../linalg.zig");

pub const Mode = enum {
    full,
    economical,
};

pub fn QR(T: type, order: Order) type {
    return struct {
        _rows: u32,
        _cols: u32,
        _qr: [*]T,
        _tau: [*]T,

        pub fn init(rows: u32, cols: u32, a: [*]T, tau: [*]T) QR(T, order) {
            return .{
                ._rows = rows,
                ._cols = cols,
                ._qr = a,
                ._tau = tau,
            };
        }

        pub fn deinit(self: *QR(T, order), allocator: std.mem.Allocator) void {
            allocator.free(self._qr[0 .. self._rows * self._cols]);
            allocator.free(self._tau[0..int.min(self._rows, self._cols)]);

            self.* = undefined;
        }

        pub fn q(self: *const QR(T, order), allocator: std.mem.Allocator, mode: Mode, ctx: anytype) !matrix.General(T, order) {
            var Q: matrix.General(T, order) = try .full(
                allocator,
                self._rows,
                if (mode == .full) self._rows else int.min(self._rows, self._cols),
                0,
                ctx,
            );
            errdefer Q.deinit(allocator);

            try linalg.lapack.lacpy(
                order,
                .{ .full = {} },
                types.scast(i32, Q.rows),
                types.scast(i32, Q.cols),
                self._qr,
                types.scast(i32, if (order == .col_major) self._rows else self._cols),
                Q.data,
                types.scast(i32, Q.ld),
                ctx,
            );

            if (comptime options.link_lapacke != null) {
                if (comptime !types.isComplex(T)) {
                    try linalg.lapack.orgqr(
                        order,
                        types.scast(i32, Q.rows),
                        types.scast(i32, if (mode == .full) Q.rows else int.min(Q.rows, Q.cols)),
                        types.scast(i32, int.min(Q.rows, Q.cols)),
                        Q.data,
                        types.scast(i32, Q.ld),
                        self._tau,
                        undefined,
                        undefined,
                        .{},
                    );
                } else {
                    try linalg.lapack.ungqr(
                        order,
                        types.scast(i32, Q.rows),
                        types.scast(i32, if (mode == .full) Q.rows else int.min(Q.rows, Q.cols)),
                        types.scast(i32, int.min(Q.rows, Q.cols)),
                        Q.data,
                        types.scast(i32, Q.ld),
                        self._tau,
                        undefined,
                        undefined,
                        .{},
                    );
                }
            } else {
                if (comptime !types.isComplex(T)) {
                    var lwork: i32 = -1;
                    try linalg.lapack.orgqr(
                        order,
                        types.scast(i32, Q.rows),
                        types.scast(i32, if (mode == .full) Q.rows else int.min(Q.rows, Q.cols)),
                        types.scast(i32, int.min(Q.rows, Q.cols)),
                        @as([*]f64, undefined),
                        types.scast(i32, Q.ld),
                        @as([*]f64, undefined),
                        @as([*]i32, @ptrCast(&lwork)),
                        -1,
                        .{},
                    );

                    var work: vector.Vector(T) = try .init(allocator, types.scast(u32, lwork));
                    defer work.deinit(allocator);

                    try linalg.lapack.orgqr(
                        order,
                        types.scast(i32, Q.rows),
                        types.scast(i32, if (mode == .full) Q.rows else int.min(Q.rows, Q.cols)),
                        types.scast(i32, int.min(Q.rows, Q.cols)),
                        Q.data,
                        types.scast(i32, Q.ld),
                        self._tau,
                        work.data,
                        lwork,
                        .{},
                    );
                } else {
                    var lwork: i32 = -1;
                    try linalg.lapack.ungqr(
                        order,
                        types.scast(i32, Q.rows),
                        types.scast(i32, if (mode == .full) Q.rows else int.min(Q.rows, Q.cols)),
                        types.scast(i32, int.min(Q.rows, Q.cols)),
                        @as([*]f64, undefined),
                        types.scast(i32, Q.ld),
                        @as([*]f64, undefined),
                        @as([*]i32, @ptrCast(&lwork)),
                        -1,
                        .{},
                    );

                    var work: vector.Vector(T) = try .init(allocator, types.scast(u32, lwork));
                    defer work.deinit(allocator);

                    try linalg.lapack.ungqr(
                        order,
                        types.scast(i32, Q.rows),
                        types.scast(i32, if (mode == .full) Q.rows else int.min(Q.rows, Q.cols)),
                        types.scast(i32, int.min(Q.rows, Q.cols)),
                        Q.data,
                        types.scast(i32, Q.ld),
                        self._tau,
                        work.data,
                        lwork,
                        .{},
                    );
                }
            }

            return Q;
        }

        pub fn r(self: *const QR(T, order), mode: Mode) matrix.Triangular(T, .upper, .non_unit, order) {
            return .{
                .data = self._qr,
                .rows = if (mode == .full)
                    self._rows
                else
                    int.min(self._rows, self._cols),
                .cols = self._cols,
                .ld = if (order == .col_major)
                    self._rows
                else
                    self._cols,
                .flags = .{ .owns_data = false },
            };
        }

        pub fn reconstruct(self: *const QR(T, order), allocator: std.mem.Allocator) !matrix.General(T, order) {
            _ = self;
            _ = allocator;
            // A = Q * R
            @compileError("QR.reconstruct not implemented yet");
        }
    };
}

pub fn QRP(T: type, order: Order) type {
    _ = T;
    _ = order;
    return void;
}

pub fn qr(allocator: std.mem.Allocator, a: anytype, ctx: anytype) !QR(Numeric(@TypeOf(a)), orderOf(@TypeOf(a))) {
    const A: type = @TypeOf(a);

    comptime if (!types.isMatrix(A) or types.matrixType(A) != .general)
        @compileError("qr: a must be a general matrix");

    comptime if (types.isArbitraryPrecision(Numeric(A))) {
        // When implemented, expand if
        @compileError("zml.linalg.qr not implemented for arbitrary precision types yet");
    } else {
        types.validateContext(@TypeOf(ctx), .{});
    };

    const m: u32 = a.rows;
    const n: u32 = a.cols;

    if (comptime options.link_lapacke != null) {
        var _qr = try a.copy(allocator, ctx);
        errdefer _qr.deinit(allocator);

        var tau: vector.Vector(Numeric(A)) = try .init(allocator, int.min(m, n));
        errdefer tau.deinit(allocator);

        try linalg.lapack.geqrf(
            types.orderOf(A),
            types.scast(i32, m),
            types.scast(i32, n),
            _qr.data,
            types.scast(i32, _qr.ld),
            tau.data,
            undefined,
            undefined,
            ctx,
        );

        return .init(m, n, _qr.data, tau.data);
    } else {
        var lwork: i32 = -1;
        try linalg.lapack.geqrf(
            types.orderOf(A),
            types.scast(i32, m),
            types.scast(i32, n),
            @as([*]f64, undefined),
            types.scast(i32, a.ld),
            @as([*]f64, undefined),
            @as([*]i32, @ptrCast(&lwork)),
            -1,
            ctx,
        );

        var _qr = try a.copy(allocator, ctx);
        errdefer _qr.deinit(allocator);

        var tau: vector.Vector(Numeric(A)) = try .init(allocator, int.min(m, n));
        errdefer tau.deinit(allocator);

        var work: vector.Vector(Numeric(A)) = try .init(allocator, types.scast(u32, lwork));
        defer work.deinit(allocator);

        try linalg.lapack.geqrf(
            types.orderOf(A),
            types.scast(i32, m),
            types.scast(i32, n),
            _qr.data,
            types.scast(i32, _qr.ld),
            tau.data,
            work.data,
            lwork,
            ctx,
        );

        return .init(m, n, _qr.data, tau.data);
    }
}

pub fn qrp() void {
    // Make sure to convert from column permutation to row permutation
    // (Permutation matrix type assumes row permutations and, in case whe
    // multiplicate from the right, like A * P, matmul already takes care of
    // converting it to column permutation internally, i.e., doing the inverse
    // permutation). Also from 1 to 0-based indexing.
}
