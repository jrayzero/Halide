
namespace Halide { namespace Runtime { namespace Internal {

#define MAKE_WORK(name, dtype) \
struct name { \
    name *next_job; \
    int (*f)(void *, dtype, uint8_t *); \
    void *user_context; \
    int next, max; \
    uint8_t *closure; \
    int active_workers; \
    int exit_status; \
    bool running() { return next < max || active_workers > 0; } \
};

MAKE_WORK(work, int)
MAKE_WORK(work_64, int64_t)

// The work queue and thread pool is weak, so one big work queue is shared by all halide functions
#define MAX_THREADS 64
#define MAKE_WORK_QUEUE(name, work_name) \
struct name { \
    /* all fields are protected by this mutex. */ \
    halide_mutex mutex; \
    /* Singly linked list for job stack */ \
    work_name *jobs; \
    /* Worker threads are divided into an 'A' team and a 'B' team. The */ \
    /* B team sleeps on the wakeup_b_team condition variable. The A */ \
    /* team does work. Threads transition to the B team if they wake */ \
    /* up and find that a_team_size > target_a_team_size.  Threads */ \
    /* move into the A team whenever they wake up and find that */ \
    /* a_team_size < target_a_team_size. */ \
    int a_team_size, target_a_team_size; \
    /* Broadcast when a job completes. */ \
    halide_cond wakeup_owners; \
    /* Broadcast whenever items are added to the work queue. */ \
    halide_cond wakeup_a_team; \
    /* May also be broadcast when items are added to the work queue if */ \
    /* more threads are required than are currently in the A team. */ \
    halide_cond wakeup_b_team; \
    /* Keep track of threads so they can be joined at shutdown*/ \
    halide_thread *threads[MAX_THREADS]; \
    /* The number threads created */ \
    int threads_created; \
    /* The desired number threads doing work. */ \
    int desired_num_threads; \
    /* Global flags indicating the threadpool should shut down, and */ \
    /* whether the thread pool has been initialized. */ \
    bool shutdown, initialized; \
    bool running() { \
        return !shutdown; \
    } \
};
MAKE_WORK_QUEUE(work_queue_t, work)
MAKE_WORK_QUEUE(work_queue_64_t, work_64)
WEAK work_queue_t work_queue;
WEAK work_queue_64_t work_queue_64;

WEAK int clamp_num_threads(int desired_num_threads) {
    if (desired_num_threads > MAX_THREADS) {
        desired_num_threads = MAX_THREADS;
    } else if (desired_num_threads < 1) {
        desired_num_threads = 1;
    }
    return desired_num_threads;
}

WEAK int default_desired_num_threads() {
    int desired_num_threads = 0;
    char *threads_str = getenv("HL_NUM_THREADS");
    if (!threads_str) {
        // Legacy name for HL_NUM_THREADS
        threads_str = getenv("HL_NUMTHREADS");
    }
    if (threads_str) {
        desired_num_threads = atoi(threads_str);
    } else {
        desired_num_threads = halide_host_cpu_count();
    }
    return desired_num_threads;
}

#define MAKE_WORKER_THREAD_ALREADY_LOCKED(name, work_name, work_queue_name, halide_do_task_name) \
WEAK void name(work_name *owned_job) { \
    /* If I'm a job owner, then I was the thread that called */ \
    /* do_par_for, and I should only stay in this function until my */ \
    /* job is complete. If I'm a lowly worker thread, I should stay in */ \
    /* this function as long as the work queue is running. */ \
    while (owned_job != NULL ? owned_job->running() \
           : (work_queue_name).running()) { \
\
        if ((work_queue_name).jobs == NULL) { \
            if (owned_job) { \
                /* There are no jobs pending. Wait for the last worker */ \
                /* to signal that the job is finished. */ \
                halide_cond_wait(&(work_queue_name).wakeup_owners, &work_queue_name.mutex); \
            } else if (work_queue_name.a_team_size <= work_queue_name.target_a_team_size) { \
                /* There are no jobs pending. Wait until more jobs are enqueued. */ \
                halide_cond_wait(&(work_queue_name).wakeup_a_team, &work_queue_name.mutex); \
            } else { \
                /* There are no jobs pending, and there are too many */ \
                /* threads in the A team. Transition to the B team */ \
                /* until the wakeup_b_team condition is fired. */ \
                (work_queue_name).a_team_size--; \
                halide_cond_wait(&(work_queue_name).wakeup_b_team, &(work_queue_name).mutex); \
                (work_queue_name).a_team_size++; \
            } \
        } else { \
            /* Grab the next job. */ \
            work_name *job = (work_queue_name).jobs; \
\
            /* Claim a task from it. */ \
            work_name myjob = *job; \
            job->next++; \
\
            /* If there were no more tasks pending for this job, */ \
            /* remove it from the stack. */ \
            if (job->next == job->max) { \
                (work_queue_name).jobs = job->next_job; \
            } \
\
            /* Increment the active_worker count so that other threads */ \
            /* are aware that this job is still in progress even */ \
            /* though there are no outstanding tasks for it. */ \
            job->active_workers++; \
\
            /* Release the lock and do the task. */ \
            halide_mutex_unlock(&(work_queue_name).mutex); \
            int result = halide_do_task_name(myjob.user_context, myjob.f, myjob.next, \
                                        myjob.closure); \
            halide_mutex_lock(&work_queue_name.mutex); \
\
            /* If this task failed, set the exit status on the job. */ \
            if (result) { \
                job->exit_status = result; \
            } \
\
            /* We are no longer active on this job */ \
            job->active_workers--; \
\
            /* If the job is done and I'm not the owner of it, wake up */ \
            /* the owner. */ \
            if (!job->running() && job != owned_job) { \
                halide_cond_broadcast(&(work_queue_name).wakeup_owners); \
            } \
        } \
    } \
} \

MAKE_WORKER_THREAD_ALREADY_LOCKED(worker_thread_already_locked, work, work_queue, halide_do_task)
MAKE_WORKER_THREAD_ALREADY_LOCKED(worker_thread_already_locked_64, work_64, work_queue_64, halide_do_task_64)

#define MAKE_WORKER_THREAD(name, work_queue_name, worker_thread_already_locked_name) \
WEAK void name(void *) {  \
    halide_mutex_lock(&(work_queue_name).mutex); \
    worker_thread_already_locked_name(NULL); \
    halide_mutex_unlock(&(work_queue_name).mutex); \
}

MAKE_WORKER_THREAD(worker_thread, work_queue, worker_thread_already_locked)
MAKE_WORKER_THREAD(worker_thread_64, work_queue_64, worker_thread_already_locked_64)


}}}  // namespace Halide::Runtime::Internal

using namespace Halide::Runtime::Internal;

extern "C" {

#define MAKE_HALIDE_DEFAULT_DO_TASK(name, halide_task_t_type, dtype) \
WEAK int name(void *user_context, halide_task_t_type f, dtype idx, \
                                uint8_t *closure) { \
    return f(user_context, idx, closure); \
}

MAKE_HALIDE_DEFAULT_DO_TASK(halide_default_do_task, halide_task_t, int)
MAKE_HALIDE_DEFAULT_DO_TASK(halide_default_do_task_64, halide_task_64_t, int64_t)

#define MAKE_HALIDE_DEFAULT_DO_PAR_FOR(name, halide_task_t_type, dtype, work_name, work_queue_name, worker_thread_already_locked_name, worker_thread_name) \
WEAK int name(void *user_context, halide_task_t_type f, \
                                   dtype min, dtype size, uint8_t *closure) { \
    /* Our for loops are expected to gracefully handle sizes <= 0 */ \
    if (size <= 0) { \
        return 0; \
    } \
\
    /* Grab the lock. If it hasn't been initialized yet, then the */ \
    /* field will be zero-initialized because it's a static global. */ \
    halide_mutex_lock(&(work_queue_name).mutex); \
\
    if (!(work_queue_name).initialized) { \
        (work_queue_name).shutdown = false; \
        halide_cond_init(&(work_queue_name).wakeup_owners); \
        halide_cond_init(&(work_queue_name).wakeup_a_team); \
        halide_cond_init(&(work_queue_name).wakeup_b_team); \
        (work_queue_name).jobs = NULL; \
\
        /* Compute the desired number of threads to use. Other code */ \
        /* can also mess with this value, but only when the work queue */ \
        /* is locked. */ \
        if (!(work_queue_name).desired_num_threads) { \
            (work_queue_name).desired_num_threads = default_desired_num_threads(); \
        } \
        (work_queue_name).desired_num_threads = clamp_num_threads((work_queue_name).desired_num_threads); \
        (work_queue_name).threads_created = 0; \
\
        /* Everyone starts on the a team. */ \
        (work_queue_name).a_team_size = work_queue_name.desired_num_threads; \
\
        (work_queue_name).initialized = true; \
    } \
\
    while ((work_queue_name).threads_created < work_queue_name.desired_num_threads - 1) { \
        /* We might need to make some new threads, if work_queue_name.desired_num_threads has */ \
        /* increased. */ \
      (work_queue_name).threads[work_queue_name.threads_created++] =	\
            halide_spawn_thread(worker_thread_name, NULL); \
    } \
\
    /* Make the job. */ \
    work_name job; \
    job.f = f;               /* The job should call this function. It takes an index and a closure. */ \
    job.user_context = user_context; \
    job.next = min;          /* Start at this index. */ \
    job.max  = min + size;   /* Keep going until one less than this index. */ \
    job.closure = closure;   /* Use this closure. */ \
    job.exit_status = 0;     /* The job hasn't failed yet */ \
    job.active_workers = 0;  /* Nobody is working on this yet */ \
\
    if (!(work_queue_name).jobs && size < (work_queue_name).desired_num_threads) { \
        /* If there's no nested parallelism happening and there are */ \
        /* fewer tasks to do than threads, then set the target A team */ \
        /* size so that some threads will put themselves to sleep */ \
        /* until a larger job arrives. */ \
        (work_queue_name).target_a_team_size = size; \
    } else { \
        /* Otherwise the target A team size is */ \
        /* desired_num_threads. This may still be less than */ \
        /* threads_created if desired_num_threads has been reduced by */ \
        /* other code. */ \
        (work_queue_name).target_a_team_size = work_queue_name.desired_num_threads; \
    } \
\
    /* Push the job onto the stack. */ \
    job.next_job = (work_queue_name).jobs; \
    (work_queue_name).jobs = &job; \
\
    /* Wake up our A team. */ \
    halide_cond_broadcast(&(work_queue_name).wakeup_a_team); \
\
    /* If there are fewer threads than we would like on the a team, */ \
    /* wake up the b team too. */ \
    if ((work_queue_name).target_a_team_size > (work_queue_name).a_team_size) { \
        halide_cond_broadcast(&(work_queue_name).wakeup_b_team); \
    } \
\
    /* Do some work myself. */ \
    worker_thread_already_locked_name(&job); \
\
    halide_mutex_unlock(&work_queue_name.mutex); \
\
    /* Return zero if the job succeeded, otherwise return the exit */ \
    /* status of one of the failing jobs (whichever one failed last). */ \
    return job.exit_status; \
}

MAKE_HALIDE_DEFAULT_DO_PAR_FOR(halide_default_do_par_for, halide_task_t, int, work, work_queue, worker_thread_already_locked, worker_thread)
MAKE_HALIDE_DEFAULT_DO_PAR_FOR(halide_default_do_par_for_64, halide_task_64_t, int64_t , work_64, work_queue_64, worker_thread_already_locked_64, worker_thread_64)

#define MAKE_HALIDE_SET_NUM_THREADS(name, work_queue_name) \
WEAK int name(int n) { \
    if (n < 0) { \
        halide_error(NULL, "halide_set_num_threads: must be >= 0."); \
    } \
    /* Don't make this an atomic swap - we don't want to be changing */ \
    /* the desired number of threads while another thread is in the */ \
    /* middle of a sequence of non-atomic operations. */ \
    halide_mutex_lock(&(work_queue_name).mutex); \
    if (n == 0) { \
        n = default_desired_num_threads(); \
    } \
    int old = (work_queue_name).desired_num_threads; \
    (work_queue_name).desired_num_threads = clamp_num_threads(n); \
    halide_mutex_unlock(&(work_queue_name).mutex); \
    return old; \
}

MAKE_HALIDE_SET_NUM_THREADS(halide_set_num_threads, work_queue)
MAKE_HALIDE_SET_NUM_THREADS(halide_set_num_threads_64, work_queue_64)

#define MAKE_HALIDE_SHUTDOWN_THREAD_POOL(name, work_queue_name) \
WEAK void name() { \
    if (!(work_queue_name).initialized) return; \
\
    /* Wake everyone up and tell them the party's over and it's time */ \
    /* to go home */ \
    halide_mutex_lock(&(work_queue_name).mutex); \
    (work_queue_name).shutdown = true; \
    halide_cond_broadcast(&(work_queue_name).wakeup_owners); \
    halide_cond_broadcast(&(work_queue_name).wakeup_a_team); \
    halide_cond_broadcast(&(work_queue_name).wakeup_b_team); \
    halide_mutex_unlock(&(work_queue_name).mutex); \
\
    /* Wait until they leave */ \
    for (int i = 0; i < (work_queue_name).threads_created; i++) { \
        halide_join_thread((work_queue_name).threads[i]); \
    } \
\
    /* Tidy up */ \
    halide_mutex_destroy(&(work_queue_name).mutex); \
    halide_cond_destroy(&(work_queue_name).wakeup_owners); \
    halide_cond_destroy(&(work_queue_name).wakeup_a_team); \
    halide_cond_destroy(&(work_queue_name).wakeup_b_team); \
    (work_queue_name).initialized = false; \
}

MAKE_HALIDE_SHUTDOWN_THREAD_POOL(halide_shutdown_thread_pool, work_queue)
MAKE_HALIDE_SHUTDOWN_THREAD_POOL(halide_shutdown_thread_pool_64, work_queue_64)

}
