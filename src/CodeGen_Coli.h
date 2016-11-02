#ifndef HALIDE_CODEGEN_COLI_H
#define HALIDE_CODEGEN_COLI_H

/** \file
 *
 * Defines an IRPrinter that emits C++ COLi code equivalent to a halide stmt
 */

#include "IRPrinter.h"
#include "Function.h"
#include "Scope.h"
#include "Schedule.h"
#include "Simplify.h"

namespace Halide {

namespace Internal {

/** This class emits C++ code equivalent to a halide Stmt. It's
 * mostly the same as an IRPrinter, but it's wrapped in a function
 * definition, and some things are handled differently to be valid
 * C++.
 */
class CodeGen_Coli : public IRPrinter {
public:
    /** Initialize a COLi code generator pointing at a particular output
     * stream (e.g. a file, or std::cout) */
    CodeGen_Coli(std::ostream &dest,
                 const std::string &pipeline_name,
                 const std::vector<Function> &outputs,
                 const std::vector<std::vector<int32_t>> &output_buffer_extents,
                 const std::vector<Type> &output_buffer_types,
                 const std::vector<std::string> &inputs,
                 const std::vector<std::vector<int32_t>> &input_buffer_extents,
                 const std::vector<Type> &input_buffer_types,
                 const std::vector<Function> &order,
                 const std::map<std::string, Schedule> &schedules);
    ~CodeGen_Coli();

    void print(Expr e);
    void print(Stmt s);

    EXPORT static void test();

private:
    struct Loop {
        std::string name;
        Expr min, extent;

        std::string to_string() const {
            std::ostringstream ss;
            Expr max = simplify(min + extent - 1);
            ss << min << " <= " << name << " <= " << max;
            return ss.str();
        }
    };

    std::string func; // Represent one Halide pipeline
    const std::vector<Function> &order;
    const std::map<std::string, Schedule> &schedules;

    std::set<std::string> output_buffers;
    std::set<std::string> input_buffers;
    Scope<Expr> scope; // Scope of the variables
    std::vector<Loop> loop_dims;
    std::set<std::string> temporary_buffers;
    std::set<std::string> computation_list;
    std::set<std::string> constant_list;

    void push_loop_dim(const std::string &name, Expr min, Expr extent);
    void pop_loop_dim();
    std::string get_loop_bound_vars() const;
    std::string get_loop_bounds() const;
    void define_constant(const std::string &name, Expr value);
    void generate_schedule();

    Expr substitute_in_lets(Expr expr) const;

protected:
    using IRPrinter::visit;

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
};

/**
 * Dump an HTML-formatted print of a Stmt to filename.
 */
EXPORT void print_to_coli(
    Stmt s, std::ostream &dest,
    const std::string &pipeline_name,
    const std::vector<Function> &outputs,
    const std::vector<std::vector<int32_t>> &output_buffer_extents,
    const std::vector<Type> &output_buffer_types,
    const std::vector<std::string> &inputs,
    const std::vector<std::vector<int32_t>> &input_buffer_extents,
    const std::vector<Type> &input_buffer_types,
    const std::vector<Function> &order,
    const std::map<std::string, Schedule> &schedules);

}
}

#endif
