#ifndef HALIDE_CODEGEN_TIRAMISU_H
#define HALIDE_CODEGEN_TIRAMISU_H

/** \file
 *
 * Defines an IRPrinter that emits C++ Tiramisu code equivalent to a halide stmt
 */

#include "IRVisitor.h"
#include "Function.h"
#include "Simplify.h"
#include "Schedule.h"


namespace Halide {

namespace Internal {

/** This class emits C++ code equivalent to a halide Stmt. It's
 * mostly the same as an IRPrinter, but it's wrapped in a function
 * definition, and some things are handled differently to be valid
 * C++.
 */
class CodeGen_Tiramisu : public IRVisitor {
public:
    /** Initialize a Tiramisu code generator pointing at a particular output
     * stream (e.g. a file, or std::cout) */
    CodeGen_Tiramisu(std::ostream &dest,
                     const std::string &pipeline_name,
                     const std::vector<Function> &outputs,
                     const std::vector<std::vector<int32_t>> &output_buffer_extents,
                     const std::vector<Type> &output_buffer_types,
                     const std::vector<std::string> &inputs,
                     const std::vector<std::vector<int32_t>> &input_buffer_extents,
                     const std::vector<Type> &input_buffer_types,
                     const std::vector<std::string> &input_params,
                     const std::vector<Type> &input_param_types,
                     const std::vector<std::string> &order,
                     const std::map<std::string, Function> &env);
    ~CodeGen_Tiramisu();

    EXPORT std::string print(Expr e);
    EXPORT void print(Stmt s);

    EXPORT static void test();

private:
    struct Loop {
        std::string name;
        Expr min, extent;
        std::string func;
        int stage;
        std::string var;

        std::string to_string() const {
            std::ostringstream ss;
            Expr max = simplify(min + extent - 1);
            ss << min << " <= " << name << " <= " << max;
            return ss.str();
        }
    };

    std::string expr;
    std::ostream &stream;
    int indent;

    std::string pipeline; // Represent one Halide pipeline
    const std::vector<std::string> &order;
    const std::map<std::string, Function> &env;

    std::vector<std::pair<std::string, Expr>> scope; // Scope of the variables
    std::vector<Loop> loop_dims;
    std::set<std::string> output_buffers;
    std::set<std::string> input_buffers;
    std::set<std::string> temporary_buffers;
    std::set<std::string> computation_list;
    std::set<std::string> constant_list;
    std::set<std::string> extent_list;
    int loop_depth;
    std::vector<std::string> buffer_str;
    std::string current_computation;

    std::string do_indent() const;

    void push_loop_dim(const std::string &name, Expr min, Expr extent,
                       const std::string &func_name, int stage, const std::string &var);
    void pop_loop_dim();
    std::string get_loop_bound_vars() const;
    std::string get_loop_bounds() const;
    std::string define_constant(const std::string &name, Expr value);
    std::string define_wrapper_let(const std::string &computation_name, const std::string &name, Expr value);
    void generate_schedules();
    void generate_schedule(const Function &func, int stage, const StageSchedule &schedule, size_t i);

    Expr substitute_in_scope(Expr expr) const;

    std::string get_current_func_name() const;
    int get_current_stage() const;

    void generate_buffer(const Realize *op);

    std::vector<std::string> get_stage_dims(const std::string &name, int stage, bool ignore_rvar) const;
    std::vector<std::string> get_stage_rvars(const std::string &name, int stage) const;

    template <typename T>
    void visit_binary(const T *op, const std::string &op_str);

protected:
    using IRVisitor::visit;

    void visit(const IntImm *);
    void visit(const UIntImm *);
    void visit(const FloatImm *);
    void visit(const StringImm *);
    void visit(const Cast *);
    void visit(const Variable *);
    void visit(const Add *);
    void visit(const Sub *);
    void visit(const Mul *);
    void visit(const Div *);
    void visit(const Mod *);
    void visit(const Min *);
    void visit(const Max *);
    void visit(const EQ *);
    void visit(const NE *);
    void visit(const LT *);
    void visit(const LE *);
    void visit(const GT *);
    void visit(const GE *);
    void visit(const And *);
    void visit(const Or *);
    void visit(const Not *);
    void visit(const Select *);
    void visit(const Load *);
    void visit(const Ramp *);
    void visit(const Broadcast *);
    void visit(const Call *);
    void visit(const Let *);
    void visit(const LetStmt *);
    void visit(const AssertStmt *);
    void visit(const ProducerConsumer *);
    void visit(const For *);
    void visit(const Store *);
    void visit(const Provide *);
    void visit(const Allocate *);
    void visit(const Free *);
    void visit(const Realize *);
    void visit(const Block *);
    void visit(const IfThenElse *);
    void visit(const Evaluate *);
    void visit(const Shuffle *);
    void visit(const Prefetch *);
};

/**
 * Dump an HTML-formatted print of a Stmt to 'dest'.
 */
EXPORT void print_to_tiramisu(
    Stmt s, std::ostream &dest,
    const std::string &pipeline_name,
    const std::vector<Function> &outputs,
    const std::vector<std::vector<int32_t>> &output_buffer_extents,
    const std::vector<Type> &output_buffer_types,
    const std::vector<std::string> &inputs,
    const std::vector<std::vector<int32_t>> &input_buffer_extents,
    const std::vector<Type> &input_buffer_types,
    const std::vector<std::string> &input_params,
    const std::vector<Type> &input_param_types,
    const std::vector<std::string> &order,
    const std::map<std::string, Function> &env);

}
}

#endif
