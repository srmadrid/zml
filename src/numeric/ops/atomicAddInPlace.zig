const std = @import("std");
const meta = @import("../../meta.zig");
const numeric = @import("../../numeric.zig");

const cache_line: usize = switch (@import("builtin").cpu.arch) {
    .aarch64, .aarch64_be, .powerpc64, .powerpc64le => 128,
    else => 64,
};

pub const SpinLock = struct {
    state: std.atomic.Value(bool) align(cache_line) = std.atomic.Value(bool).init(false),
    _pad: [cache_line - @sizeOf(std.atomic.Value(bool))]u8 = undefined,

    pub inline fn lock(self: *SpinLock) void {
        if (self.state.cmpxchgWeak(false, true, .acquire, .monotonic) != null) {
            @branchHint(.unlikely);

            self.lockSlow();
        }
    }

    fn lockSlow(self: *SpinLock) void {
        @branchHint(.cold);

        while (true) {
            while (self.state.load(.monotonic))
                std.atomic.spinLoopHint();

            if (self.state.cmpxchgWeak(false, true, .acquire, .monotonic) == null)
                return;
        }
    }

    pub inline fn unlock(self: *SpinLock) void {
        self.state.store(false, .release);
    }
};

fn TypeLocks(comptime T: type) type {
    return struct {
        const pool_size: usize = if (@sizeOf(T) > 64) 256 else 1024;
        const golden: usize = if (@sizeOf(usize) == 8) 0x9E3779B97F4A7C15 else 0x9E3779B9;
        const shift_amt = @bitSizeOf(usize) - @as(comptime_int, std.math.log2_int(usize, pool_size));

        var locks: [pool_size]SpinLock align(cache_line) = .{SpinLock{}} ** pool_size;

        inline fn get(ptr: *T) *SpinLock {
            // Fibonacci hash: multiply by 2^N/φ, keep the high log2(pool_size) bits.
            const addr = @intFromPtr(ptr);
            return &locks[(addr *% golden) >> shift_amt];
        }
    };
}

pub inline fn atomicAddInPlace(ptr: anytype, val: anytype) void {
    comptime var P = @TypeOf(ptr);
    const V = @TypeOf(val);

    comptime if (!meta.isPointer(P) or meta.isConstPointer(P) or !meta.isNumeric(meta.Child(P)) or
        !meta.isNumeric(V))
        @compileError("zsl.numeric.atomicAddInPlace: ptr must be a mutable one-item pointer to a numeric, and val must be a numeric, got \n\tptr: " ++ @typeName(P) ++ "\n\tval: " ++ @typeName(V) ++ "\n");

    P = meta.Child(P);

    switch (comptime meta.numericType(P)) {
        .int, .float => _ = @atomicRmw(P, ptr, .Add, val, .monotonic),
        else => {
            const lock = TypeLocks(P).get(ptr);

            lock.lock();
            defer lock.unlock();

            numeric.addInto(ptr, ptr.*, val);
        },
    }
}
