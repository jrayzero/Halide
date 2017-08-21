#include "ScheduleFunctions.h"
#include "IROperator.h"
#include "Simplify.h"
#include "Substitute.h"
#include "ExprUsesVar.h"
#include "Var.h"
#include "Qualify.h"
#include "IRMutator.h"
#include "Target.h"
#include "Inline.h"
#include "CodeGen_GPU_Dev.h"
#include "IRPrinter.h"
#include "Func.h"
#include "ApplySplit.h"
#include "IREquality.h"

#include <iostream>

namespace Halide {
namespace Internal {

using std::string;
using std::map;
using std::vector;
using std::pair;
using std::set;
using std::ostringstream;

namespace {

string print_name(const string &name) {
    ostringstream oss;

    // Prefix an underscore to avoid reserved words (e.g. a variable named "while")
    if (isalpha(name[0])) {
        oss << '_';
    }

    for (size_t i = 0; i < name.size(); i++) {
        if (name[i] == '.') {
            oss << '_';
        } else if (name[i] == '$') {
            oss << "__";
        } else if (name[i] != '_' && !isalnum(name[i])) {
            oss << "___";
        }
        else oss << name[i];
    }
    return oss.str();
}

bool var_name_match(string candidate, string var) {
    internal_assert(var.find('.') == string::npos)
        << "var_name_match expects unqualified names for the second argument. "
        << "Name passed: " << var << "\n";
    return (candidate == var) || Internal::ends_with(candidate, "." + var);
}

int get_dim_index(const vector<Dim> &dims, const string &var) {
    const auto &iter = std::find_if(dims.begin(), dims.end(),
        [&var](const Dim &d) { return var_name_match(d.var, var); });
    internal_assert(iter != dims.end());
    return (iter - dims.begin());
}

string generate_split_schedule(const string &func_name, const string &prefix,
                               const Split &split, vector<string> &dims) {
    string schedule;
    // TODO(tiramisu): Rename and Fuse are not currently supported by Tiramisu,
    // so throw an error for now.
    if (split.is_split()) {
        int old_idx = get_dim_index(dims, split.old_var);
        int outer = get_dim_index(dims, split.outer);
        int inner = get_dim_index(dims, split.inner);

        const int64_t *factor = as_const_int(split.factor);
        internal_assert(factor);




        // TODO(tiramisu): Currently Tiramisu ignores the split strategy and
        // uses whatever the default is in
    } else if (split.is_fuse()) {
        user_error << "Fuse is not currently supported by Tiramisu\n";
    } else if (split.is_rename()) {
        user_error << "Rename is not currently supported by Tiramisu\n";
    }
    // Do nothing for purify

    return result;
}

// Build a loop nest about a provide node using a schedule
Stmt provide_loop_nest_schedule_helper(string func_name,
                                       string prefix,
                                       const vector<string> &dims,
                                       const vector<Expr> &site,
                                       const vector<Expr> &values,
                                       const vector<Expr> &predicates,
                                       const FuncSchedule &func_s,
                                       const StageSchedule &stage_s,
                                       bool is_update) {

    // TODO(tiramisu): Handle Halide func_s.bounds()

    // Make the (multi-dimensional multi-valued) store node.
    Stmt stmt = Provide::make(func_name, values, site); // TODO(tiramisu): computation expression (buffer mapping ???)

    vector<Split> splits = stage_s.splits();

    vector<string> pure_loop_dims;

    // Dim list in Halide is defined from innermost to outermost, while in
    // Tiramisu, it is defined from outermost to innermost.
    for (int i = (int)dims.size() - 1; i >= 0; --i) {
        pure_loop_dims.push_back(print_name(prefix + dims[i]));
    }

    // Put all the reduction domain predicates into the containers vector.
    for (Expr pred : predicates) {
        pred = qualify(prefix, pred);
        Container c = {Container::If, 0, "", likely(pred)};
        pred_container.push_back(c);
    }
    int n_predicates = pred_container.size();

    // Define the loop mins and extents for the reduction domain (if there is any)
    // in terms of the mins and maxs produced by bounds inference
    for (const ReductionVariable &rv : stage_s.rvars()) {
        string p = prefix + rv.var;
        Expr rmin = Variable::make(Int(32), p + ".min");
        Expr rmax = Variable::make(Int(32), p + ".max");
        stmt = LetStmt::make(p + ".loop_min", rmin, stmt);
        stmt = LetStmt::make(p + ".loop_max", rmax, stmt);
        stmt = LetStmt::make(p + ".loop_extent", rmax - rmin + 1, stmt);
    }

    return stmt;
}

// Build a loop nest about a provide node using a schedule
Stmt provide_loop_nest_schedule(string func_name,
                                string prefix,
                                const vector<string> &dims,
                                const FuncSchedule &f_sched,
                                const Definition &def,
                                bool is_update) {

    internal_assert(!is_update == def.is_init());

    // Default stored values
    vector<Expr> site(def.args().size());
    vector<Expr> values(def.values().size());
    for (size_t i = 0; i < values.size(); i++) {
        Expr v = def.values()[i];
        v = qualify(prefix, v);
        values[i] = v;
        debug(3) << "Value " << i << " = " << v << "\n";
    }

    // Default stored locations
    for (size_t i = 0; i < def.args().size(); i++) {
        Expr s = def.args()[i];
        s = qualify(prefix, s);
        site[i] = s;
        debug(3) << "Site " << i << " = " << s << "\n";
    }

    // Default schedule/values if there is no specialization
    Stmt stmt = provide_loop_nest_schedule_helper(
        func_name, prefix, dims, site, values, def.split_predicate(), f_sched,
        def.schedule(), is_update);

    // TODO(tiramisu): Handle specializations
    user_assert(def.specializations().empty())
        << "Specialization is not currently supported by Tiramisu.\n";

    return stmt;
}


vector<string> produce_schedule(Function f, const Target &t) {
    if (f.has_extern_definition()) {
        // TODO(tiramisu): Handle extern definition. Need to take care of
        // buffer allocation, annotation, etc.
        user_error << "Function \"" << f.name() << "\" is an extern function, "
                   << "which is not currently supported by Tiramisu.\n";
    } else {
        string prefix = f.name() + ".s0.";
        vector<string> dims = f.args(); // Pure dimensions
        return provide_loop_nest_schedule(f.name(), prefix, dims, f.schedule(),
                                          f.definition(), false);
    }
}

// Build the loop nests that update a function (assuming it's a reduction).
vector<Stmt> update_schedule(Function f) {
    vector<Stmt> updates;

    for (size_t i = 0; i < f.updates().size(); i++) {
        const Definition &def = f.update(i);

        string prefix = f.name() + ".s" + std::to_string(i+1) + ".";

        vector<string> dims = f.args();
        Stmt loop = provide_loop_nest_schedule(f.name(), prefix, dims, f.schedule(),
                                               def, true);
        updates.push_back(loop);
    }

    return updates;
}

pair<Stmt, Stmt> production_schedule(Function func, const Target &target) {
    Stmt produce = produce_schedule(func, target);
    vector<Stmt> updates = update_schedule(func);

    // Combine the update steps
    Stmt merged_updates = Block::make(updates);
    return { produce, merged_updates };
}

vector<string> generate_func_tiramisu_schedule(const Function &f, bool is_output, const Target &target) {

    const LoopLevel &compute_level = f.schedule().compute_level();
    const LoopLevel &store_level = f.schedule().store_level();

    // Generate tiramisu::computation including its schedule
    pair<Stmt, Stmt> realization = production_schedule(f, target);

    // Generate tiramisu::buffer
    if (!is_output) {
        Region bounds;
        string name = f.name();
        const vector<string> func_args = f.args();
        for (int i = 0; i < f.dimensions(); i++) {
            const string &arg = func_args[i];
            Expr min = Variable::make(Int(32), name + "." + arg + ".min_realized");
            Expr extent = Variable::make(Int(32), name + "." + arg + ".extent_realized");
            bounds.push_back(Range(min, extent));
        }

        s = Realize::make(name, f.output_types(), bounds, const_true(), s);
    }
}


bool inline_functions_in_env(const vector<Function> &outputs,
                             const vector<string> &order,
                             map<string, Function> &env) {
    bool is_inlined = false;
    // The very last few functions in 'order' are the last to be realized in the
    // pipeline (the final producers) so there is no point in checking it.
    for (int i = 0; i < (int)order.size() - (int)outputs.size(); ++i) {
        bool is_output = false;
        for (const Function &f : outputs) {
            if (order[i] == f.name()) {
                is_output = true;
                break;
            }
        }
        if (is_output) {
            // Should not inline output Func
            debug(5) << "Skip inlining " << order[i] << " since it is an output\n";
            continue;
        }
        Function f1 = env.at(order[i]);
        if (f.can_be_inlined() &&
            f.schedule().compute_level().is_inline()) {
            is_inlined = true;
            debug(4) << "Function \"" << order[i] << "\" is trivial to inline\n";
            for (int j = i + 1; j < (int)order.size() - (int)outputs.size(); ++j) {
                internal_assert(order[i] != order[j]);
                Function f2 = env.at(order[j]);
                debug(5) << "Inline trivial function \"" << f1.name()
                         << "\" inside \"" << f2.name() << "\"\n";
                inline_function(f2, f1);
            }
        }
    }
    return is_inlined;
}

} // anonymous namespace

map<string, vector<string>> schedule_functions_tiramisu(const vector<Function> &outputs,
                                                        const Target &target
                                                        vector<string> order,
                                                        map<string, Function> env) {

    // Inline inlined functions into the functions definition and
    // remove the inlined functions from the environment and
    // recompute the realization order.
    if (inline_functions_in_env(outputs, order, env)) {
        env.clear();
        for (Function f : outputs) {
            map<string, Function> more_funcs = find_transitive_calls(f);
            env.insert(more_funcs.begin(), more_funcs.end());
        }
        order = realization_order(outputs, env);
    }

    map<string, vector<string>> schedules;

    for (size_t i = order.size(); i > 0; i--) {
        Function f = env.find(order[i-1])->second;

        bool is_output = false;
        for (Function o : outputs) {
            is_output |= o.same_as(f);
        }

        internal_assert(!f.can_be_inlined() || !f.schedule().compute_level().is_inline());

        vector<string> f_schedule = generate_func_tiramisu_schedule(f, is_output, target);
        schedules.emplace(f.name(), f_schedule);
    }

    return schedules;
}

}
}
