#ifndef HALIDE_INTERNAL_LOWER_H
#define HALIDE_INTERNAL_LOWER_H

/** \file
 *
 * Defines the function that generates a statement that computes a
 * Halide function using its schedule.
 */

#include "IR.h"
#include "Target.h"

namespace Halide {
namespace Internal {

class IRMutator;

/** Given a halide function with a schedule, create a statement that
 * evaluates it. Automatically pulls in all the functions f depends
 * on. Some stages of lowering may be target-specific. */
EXPORT Stmt lower(std::vector<Function> outputs, const std::string &pipeline_name, const Target &t,
                  std::vector<Function> &func_order, std::map<std::string, Schedule> &schedules,
                  const std::vector<IRMutator *> &custom_passes = std::vector<IRMutator *>(),
                  bool compile_to_coli = false);

void lower_test();

}
}

#endif
