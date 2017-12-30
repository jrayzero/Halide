#include "HalideRuntime.h"
#include "thread_pool_common.h"

extern "C" {

namespace {
__attribute__((destructor))
WEAK void halide_thread_pool_cleanup() {
    halide_shutdown_thread_pool();
    //    halide_64bit_shutdown_thread_pool();
}
}

  /*namespace Halide { namespace Runtime { namespace Internal {
      WEAK halide_do_task_t custom_do_task = halide_default_do_task;
WEAK halide_do_par_for_t custom_do_par_for = halide_default_do_par_for;
}}}

WEAK halide_do_task_t halide_set_custom_do_task(halide_do_task_t f) {
    halide_do_task_t result = custom_do_task;
    custom_do_task = f;
    return result;
}

WEAK halide_do_par_for_t halide_set_custom_do_par_for(halide_do_par_for_t f) {
    halide_do_par_for_t result = custom_do_par_for;
    custom_do_par_for = f;
    return result;
}

WEAK int halide_do_task(void *user_context, halide_task_t f, int idx,
                        uint8_t *closure) {
    return (*custom_do_task)(user_context, f, idx, closure);
}

WEAK int halide_do_par_for(void *user_context, halide_task_t f,
                           int min, int size, uint8_t *closure) {
  return (*custom_do_par_for)(user_context, f, min, size, closure);
  }*/

WEAK int halide_64bit_do_task(void *user_context, halide_64bit_task_t f, size_t idx,
                        uint8_t *closure) {
    return halide_64bit_default_do_task(user_context, f, idx, closure);
}

WEAK int halide_64bit_do_par_for(void *user_context, halide_64bit_task_t f,
                          size_t min, size_t size, uint8_t *closure) {
  return halide_64bit_default_do_par_for(user_context, f, min, size, closure);
}

} // extern "C"
