const std = @import("std");
const builtin = @import("builtin");

/// Defines the maximum number of concurrent chunks allowed, capping the size of
/// the stack-allocated job array used during parallel iterations.
pub const max_workers: usize = 128;

pub const Task = struct {
    /// An optional pointer linking to the next task in the pool's intrusive
    /// LIFO queue.
    next: ?*Task = null,
    /// The function pointer that executes the actual work payload, receiving
    /// the task pointer and the executing worker's id.
    runFn: *const fn (*Task, worker_id: usize) void,
};

const Futex = struct {
    /// Suspends the current thread until the atomic value changes from the
    /// expect parameter.
    pub fn wait(ptr: *const std.atomic.Value(u32), expect: u32) void {
        if (builtin.os.tag == .windows)
            _ = std.os.windows.ntdll.RtlWaitOnAddress(
                @as(?*const anyopaque, @ptrCast(ptr)),
                @as(?*const anyopaque, @ptrCast(&expect)),
                4,
                null,
            )
        else if (builtin.os.tag.isDarwin())
            _ = std.c.__ulock_wait(
                .{ .op = .COMPARE_AND_WAIT, .NO_ERRNO = true },
                @ptrCast(ptr),
                expect,
                0,
            )
        else if (builtin.os.tag == .linux)
            _ = std.os.linux.futex_4arg(
                @ptrCast(&ptr.raw),
                .{ .cmd = .WAIT, .private = true },
                expect,
                null,
            )
        else
            @compileError("zsl.Futex.wait: not implemented for " ++ @tagName(builtin.os.tag));
    }

    /// Wakes up to `max_waiters` threads that are currently sleeping on the
    /// given atomic value pointer.
    pub fn wake(ptr: *const std.atomic.Value(u32), max_waiters: u32) void {
        if (max_waiters == 0)
            return;

        if (builtin.os.tag == .windows) {
            if (max_waiters == 1)
                std.os.windows.ntdll.RtlWakeAddressSingle(@as(?*const anyopaque, @ptrCast(ptr)))
            else
                std.os.windows.ntdll.RtlWakeAddressAll(@as(?*const anyopaque, @ptrCast(ptr)));
        } else if (builtin.os.tag.isDarwin()) {
            _ = std.c.__ulock_wake(
                .{ .op = .COMPARE_AND_WAIT, .NO_ERRNO = true, .WAKE_ALL = max_waiters > 1 },
                @as(*const anyopaque, @ptrCast(ptr)),
                0,
            );
        } else if (builtin.os.tag == .linux) {
            _ = std.os.linux.futex_3arg(
                @ptrCast(&ptr.raw),
                .{ .cmd = .WAKE, .private = true },
                max_waiters,
            );
        } else {
            @compileError("zsl.Futex.wake: not implemented for " ++ @tagName(builtin.os.tag));
        }
    }
};

pub const Semaphore = struct {
    /// An atomic counter tracking the number of available permits for the
    /// semaphore.
    state: std.atomic.Value(u32) = std.atomic.Value(u32).init(0),

    /// Increments the semaphore's internal state and wakes up a single waiting
    /// thread via a futex.
    pub fn post(self: *Semaphore) void {
        _ = self.state.fetchAdd(1, .release);
        Futex.wake(&self.state, 1);
    }

    /// Loops to acquire a permit by decrementing the state, putting the thread
    /// to sleep via a futex if the count is currently zero.
    pub fn wait(self: *Semaphore) void {
        var current = self.state.load(.monotonic);
        while (true) {
            if (current == 0) {
                Futex.wait(&self.state, 0);
                current = self.state.load(.monotonic);
                continue;
            }

            if (self.state.cmpxchgWeak(current, current - 1, .acquire, .monotonic)) |updated|
                current = updated
            else
                return;
        }
    }
};

pub const WaitGroup = struct {
    /// An atomic integer tracking the number of pending operations that must
    /// complete before the wait group is satisfied.
    counter: std.atomic.Value(usize) = std.atomic.Value(usize).init(0),
    /// An atomic flag used as a synchronization point to block and wake waiting
    /// threads.
    event_state: std.atomic.Value(u32) = std.atomic.Value(u32).init(0),

    /// Adds `n` pending operations to the wait group's counter and resets the
    /// event state.
    pub fn start(self: *WaitGroup, n: usize) void {
        self.event_state.store(0, .monotonic);
        _ = self.counter.fetchAdd(n, .monotonic);
    }

    /// Decrements the operation counter, waking all threads blocked on the wait
    /// group's event state if the counter reaches zero.
    pub fn finish(self: *WaitGroup) void {
        if (self.counter.fetchSub(1, .release) == 1) {
            self.event_state.store(1, .release);
            Futex.wake(&self.event_state, std.math.maxInt(u32));
        }
    }

    /// Returns `true` if all registered operations have finished and the
    /// counter is zero.
    pub fn isDone(self: *WaitGroup) bool {
        return self.counter.load(.acquire) == 0;
    }

    /// Blocks the calling thread until the wait group's event state indicates
    /// all operations have completed.
    pub fn wait(self: *WaitGroup) void {
        while (self.event_state.load(.acquire) == 0) {
            Futex.wait(&self.event_state, 0);
        }
    }
};

pub const Pool = struct {
    /// An atomic pointer representing the top of the intrusive LIFO stack for
    /// scheduled tasks.
    head: std.atomic.Value(?*Task) align(std.atomic.cache_line) = .init(null),
    /// An atomic counter tracking how many worker threads are currently
    /// sleeping and waiting for work.
    idle: std.atomic.Value(u32) align(std.atomic.cache_line) = .init(0),
    /// An atomic boolean flag that dictates whether worker threads should
    /// continue executing or shut down.
    running: std.atomic.Value(bool) = .init(true),
    /// A synchronization primitive used to wake idle worker threads when new
    /// tasks are scheduled.
    semaphore: Semaphore = .{},
    /// A dynamically allocated slice holding the active thread handles for the
    /// pool.
    threads: []std.Thread,

    pub const Opts = struct {
        /// An optional configuration field specifying the number of threads to
        /// use, defaulting to the system's cpu count if left null.
        n_jobs: ?usize = null,
    };

    /// Initializes the pool state and spawns the specified number of worker
    /// threads minus one, reserving one slot for the calling thread.
    pub fn init(pool: *Pool, allocator: std.mem.Allocator, opts: Opts) !void {
        pool.* = .{ .threads = &.{} };

        if (builtin.single_threaded)
            return;

        const n_jobs = opts.n_jobs orelse (std.Thread.getCpuCount() catch 1);
        const n_workers = if (n_jobs > 1) n_jobs - 1 else 0;

        pool.threads = try allocator.alloc(std.Thread, n_workers);
        var spawned: usize = 0;
        errdefer pool.shutdown(spawned);

        for (pool.threads, 0..) |*t, i| {
            t.* = try std.Thread.spawn(.{}, workerLoop, .{ pool, i + 1 });
            spawned += 1;
        }
    }

    /// Signals the pool to shut down, joins all worker threads, and frees the
    /// allocated thread slice.
    pub fn deinit(pool: *Pool, allocator: std.mem.Allocator) void {
        pool.shutdown(pool.threads.len);
        allocator.free(pool.threads);
        pool.* = undefined;
    }

    /// Flags the pool as inactive and posts to the semaphore enough times to
    /// wake and terminate all spawned threads.
    fn shutdown(pool: *Pool, spawned: usize) void {
        if (builtin.single_threaded)
            return;

        pool.running.store(false, .release);
        for (0..spawned) |_|
            pool.semaphore.post();

        for (pool.threads[0..spawned]) |t|
            t.join();
    }

    /// Returns the total number of participant threads, which includes the
    /// calling thread.
    pub fn idCount(pool: *Pool) usize {
        return pool.threads.len + 1;
    }

    /// Pushes a new task onto the atomic LIFO queue and wakes an idle thread if
    /// one is available.
    pub fn schedule(pool: *Pool, task: *Task) void {
        var head = pool.head.load(.monotonic);
        while (true) {
            task.next = head;
            head = pool.head.cmpxchgWeak(head, task, .release, .monotonic) orelse break;
        }

        if (pool.idle.load(.seq_cst) > 0)
            pool.semaphore.post();
    }

    /// Atomically attempts to remove and return the most recently scheduled
    /// task from the LIFO queue.
    fn pop(pool: *Pool) ?*Task {
        var head = pool.head.load(.acquire);
        while (head) |task| {
            head = pool.head.cmpxchgWeak(head, task.next, .acquire, .acquire) orelse return task;
        }

        return null;
    }

    /// Continuously pops and executes tasks from the pool's queue until the
    /// specified `WaitGroup` completes, sleeping if no tasks are available.
    pub fn wait(pool: *Pool, wg: *WaitGroup) void {
        while (!wg.isDone()) {
            if (pool.pop()) |task| {
                task.runFn(task, 0);
            } else {
                wg.wait();
                return;
            }
        }
    }

    /// The main execution loop for spawned threads, which briefly spin-waits
    /// for incoming tasks before sleeping on the semaphore to save cpu cycles.
    fn workerLoop(pool: *Pool, id: usize) void {
        const spin_limit = 40;
        var spins: u32 = 0;

        while (pool.running.load(.acquire)) {
            if (pool.pop()) |task| {
                spins = 0;
                task.runFn(task, id);
                continue;
            }

            if (spins < spin_limit) {
                spins += 1;
                std.atomic.spinLoopHint();
                continue;
            }

            _ = pool.idle.fetchAdd(1, .seq_cst);
            if (pool.pop()) |task| {
                _ = pool.idle.fetchSub(1, .monotonic);
                spins = 0;
                task.runFn(task, id);
                continue;
            }

            pool.semaphore.wait();
            _ = pool.idle.fetchSub(1, .monotonic);
        }
    }

    /// Divides a range of items into evenly sized chunks and schedules them
    /// across the pool, actively blocking and helping execute tasks until all
    /// chunks finish.
    pub fn parallelFor(
        pool: *Pool,
        len: usize,
        ctx: anytype,
        comptime kernel: fn (ctx: @TypeOf(ctx), start: usize, end: usize, worker_id: usize) void,
    ) void {
        if (len == 0)
            return;

        if (builtin.single_threaded or pool.threads.len == 0) {
            kernel(ctx, 0, len, 0);
            return;
        }

        const Ctx = @TypeOf(ctx);
        const Job = struct {
            task: Task = .{ .runFn = run },
            wg: *WaitGroup,
            ctx: Ctx,
            start: usize,
            end: usize,

            fn run(task: *Task, worker_id: usize) void {
                const self: *@This() = @alignCast(@fieldParentPtr("task", task));
                kernel(self.ctx, self.start, self.end, worker_id);
                self.wg.finish();
            }
        };

        const n_chunks = @min(@min(pool.idCount(), max_workers), len);
        const chunk = (len + n_chunks - 1) / n_chunks;

        var jobs: [max_workers]Job = undefined;
        var wg: WaitGroup = .{};
        wg.start(n_chunks);

        var start: usize = 0;
        var i: usize = 0;
        while (i < n_chunks) : (i += 1) {
            const end = @min(start + chunk, len);
            jobs[i] = .{ .wg = &wg, .ctx = ctx, .start = start, .end = end };
            pool.schedule(&jobs[i].task);
            start = end;
        }

        pool.wait(&wg);
    }
};
