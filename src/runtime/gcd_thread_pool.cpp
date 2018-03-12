#include "HalideRuntime.h"

extern "C" {

typedef long dispatch_once_t;
extern void dispatch_once_f(dispatch_once_t *, void *context, void (*initializer)(void *));

typedef struct dispatch_queue_s *dispatch_queue_t;
typedef long dispatch_queue_priority_t;

extern dispatch_queue_t dispatch_get_global_queue(
    dispatch_queue_priority_t priority, unsigned long flags);

extern void dispatch_apply_f(size_t iterations, dispatch_queue_t queue,
                             void *context, void (*work)(void *, size_t));
extern void dispatch_async_f(dispatch_queue_t queue, void *context, void (*work)(void *));

typedef struct dispatch_semaphore_s *dispatch_semaphore_t;
typedef uint64_t dispatch_time_t;
#define DISPATCH_TIME_FOREVER (~0ull)

extern dispatch_semaphore_t dispatch_semaphore_create(long value);
extern long dispatch_semaphore_wait(dispatch_semaphore_t dsema, dispatch_time_t timeout);
extern long dispatch_semaphore_signal(dispatch_semaphore_t dsema);
extern void dispatch_release(void *object);

}

namespace Halide { namespace Runtime { namespace Internal {
struct spawned_thread {
    void (*f)(void *);
    void *closure;
    dispatch_semaphore_t join_semaphore;
};
WEAK void spawn_thread_helper(void *arg) {
    spawned_thread *t = (spawned_thread *)arg;
    t->f(t->closure);
    dispatch_semaphore_signal(t->join_semaphore);
}
}}} // namespace Halide::Runtime::Internal


WEAK halide_thread *halide_spawn_thread(void (*f)(void *), void *closure) {
    spawned_thread *thread = (spawned_thread *)malloc(sizeof(spawned_thread));
    thread->f = f;
    thread->closure = closure;
    thread->join_semaphore = dispatch_semaphore_create(0);
    dispatch_async_f(dispatch_get_global_queue(0, 0), thread, spawn_thread_helper);
    return (halide_thread *)thread;
}

WEAK void halide_join_thread(halide_thread *thread_arg) {
    spawned_thread *thread = (spawned_thread *)thread_arg;
    dispatch_semaphore_wait(thread->join_semaphore, DISPATCH_TIME_FOREVER);
    free(thread);
}

// Join thread and condition variables intentionally unimplemented for
// now on OS X. Use of them will result in linker errors. Currently
// only the common thread pool uses them.

namespace Halide { namespace Runtime { namespace Internal {

WEAK int custom_num_threads = 0;

struct gcd_mutex {
    dispatch_once_t once;
    dispatch_semaphore_t semaphore;
};

WEAK void init_mutex(void *mutex_arg) {
    gcd_mutex *mutex = (gcd_mutex *)mutex_arg;
    mutex->semaphore = dispatch_semaphore_create(1);
}

#define MAKE_HALIDE_GCD_JOB(name, dtype) \
struct name { \
    int (*f)(void *, dtype, uint8_t *); \
    void *user_context; \
    uint8_t *closure; \
    dtype min; \
    int exit_status; \
};

MAKE_HALIDE_GCD_JOB(halide_gcd_job, int)
MAKE_HALIDE_GCD_JOB(halide_gcd_job_64, int64_t)

// Take a call from grand-central-dispatch's parallel for loop, and
// make a call to Halide's do task
#define MAKE_HALIDE_DO_GCD_TASK(name, job_name, dtype, halide_do_task_name) \
WEAK void name(void *job, size_t idx) { \
    job_name *j = (job_name *)job; \
    j->exit_status = halide_do_task_name(j->user_context, j->f, j->min + (dtype)idx, \
                                    j->closure); \
}

MAKE_HALIDE_DO_GCD_TASK(halide_do_gcd_task, halide_gcd_job, int, halide_do_task)
MAKE_HALIDE_DO_GCD_TASK(halide_do_gcd_task_64, halide_gcd_job_64, int64_t, halide_do_task_64)

}}}  // namespace Halide::Runtime::Internal

extern "C" {

#define MAKE_HALIDE_DEFAULT_DO_TASK(name, dtype) \
WEAK int name(void *user_context, int (*f)(void *, dtype, uint8_t *), \
                                dtype idx, uint8_t *closure) { \
    return f(user_context, idx, closure); \
}
MAKE_HALIDE_DEFAULT_DO_TASK(halide_default_do_task, int)
MAKE_HALIDE_DEFAULT_DO_TASK(halide_default_do_task_64, int64_t)

#define MAKE_HALIDE_DEFAULT_DO_PAR_FOR(name, dtype, halide_task_t_name, halide_do_task_name, halide_gcd_job_name, halide_do_gcd_task_name) \
WEAK int name(void *user_context, halide_task_t_name f, \
                                   dtype min, dtype size, uint8_t *closure) { \
    if (custom_num_threads == 1 || size == 1) { \
        /* GCD doesn't really allow us to limit the threads, */ \
        /* so ensure that there's no parallelism by executing serially. */ \
        for (dtype x = min; x < min + size; x++) { \
            int result = halide_do_task_name(user_context, f, x, closure); \
            if (result) { \
                return result; \
            } \
        } \
        return 0; \
    } \
\
    halide_gcd_job_name job; \
    job.f = f; \
    job.user_context = user_context; \
    job.closure = closure; \
    job.min = min; \
    job.exit_status = 0; \
\
    dispatch_apply_f(size, dispatch_get_global_queue(0, 0), &job, &halide_do_gcd_task_name); \
    return job.exit_status; \
}

MAKE_HALIDE_DEFAULT_DO_PAR_FOR(halide_default_do_par_for, int, halide_task_t, halide_do_task, halide_gcd_job, halide_do_gcd_task)
MAKE_HALIDE_DEFAULT_DO_PAR_FOR(halide_default_do_par_for_64, int64_t, halide_task_64_t, halide_do_task_64,
                               halide_gcd_job_64, halide_do_gcd_task_64)

}  // extern "C"

namespace Halide { namespace Runtime { namespace Internal {

WEAK halide_do_task_t custom_do_task = halide_default_do_task;
WEAK halide_do_task_64_t custom_do_task_64 = halide_default_do_task_64;
WEAK halide_do_par_for_t custom_do_par_for = halide_default_do_par_for;
WEAK halide_do_par_for_64_t custom_do_par_for_64 = halide_default_do_par_for_64;

}}}  // namespace Halide::Runtime::Internal

extern "C" {

WEAK void halide_mutex_destroy(halide_mutex *mutex_arg) {
    gcd_mutex *mutex = (gcd_mutex *)mutex_arg;
    if (mutex->once != 0) {
        dispatch_release(mutex->semaphore);
        memset(mutex_arg, 0, sizeof(halide_mutex));
    }
}

WEAK void halide_mutex_lock(halide_mutex *mutex_arg) {
    gcd_mutex *mutex = (gcd_mutex *)mutex_arg;
    dispatch_once_f(&mutex->once, mutex, init_mutex);
    dispatch_semaphore_wait(mutex->semaphore, DISPATCH_TIME_FOREVER);
}

WEAK void halide_mutex_unlock(halide_mutex *mutex_arg) {
    gcd_mutex *mutex = (gcd_mutex *)mutex_arg;
    dispatch_semaphore_signal(mutex->semaphore);
}

WEAK void halide_shutdown_thread_pool() {
}

WEAK int halide_set_num_threads(int n) {
    if (n < 0) {
        halide_error(NULL, "halide_set_num_threads: must be >= 0.");
    }
    int old_custom_num_threads = custom_num_threads;
    custom_num_threads = n;
    return old_custom_num_threads;
}

WEAK halide_do_task_t halide_set_custom_do_task(halide_do_task_t f) {
    halide_do_task_t result = custom_do_task;
    custom_do_task = f;
    return result;
}

WEAK halide_do_task_64_t halide_set_custom_do_task_64(halide_do_task_64_t f) {
    halide_do_task_64_t result = custom_do_task_64;
    custom_do_task_64 = f;
    return result;
}

WEAK halide_do_par_for_t halide_set_custom_do_par_for(halide_do_par_for_t f) {
    halide_do_par_for_t result = custom_do_par_for;
    custom_do_par_for = f;
    return result;
}

WEAK halide_do_par_for_64_t halide_set_custom_do_par_for_64(halide_do_par_for_64_t f) {
    halide_do_par_for_64_t result = custom_do_par_for_64;
    custom_do_par_for_64 = f;
    return result;
}

WEAK int halide_do_task(void *user_context, halide_task_t f, int idx,
                        uint8_t *closure) {
    return (*custom_do_task)(user_context, f, idx, closure);
}

WEAK int halide_do_task_64(void *user_context, halide_task_64_t f, int64_t idx,
                        uint8_t *closure) {
    return (*custom_do_task_64)(user_context, f, idx, closure);
}

WEAK int halide_do_par_for(void *user_context, halide_task_t f,
                           int min, int size, uint8_t *closure) {
  return (*custom_do_par_for)(user_context, f, min, size, closure);
}

WEAK int halide_do_par_for_64(void *user_context, halide_task_64_t f,
                           int64_t min, int64_t size, uint8_t *closure) {
    return (*custom_do_par_for_64)(user_context, f, min, size, closure);
}

}
