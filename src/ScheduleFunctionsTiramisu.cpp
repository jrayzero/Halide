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

vector<ApplySplitResult> generate_split_schedule(const Split &split, bool is_update, string prefix) {
    // TODO(tiramisu): Rename and Fuse are not currently supported by Tiramisu,
    // so throw an error for now.
    if (split.is_split()) {
        Expr inner = Variable::make(Int(32), prefix + split.inner);
        Expr old_max = Variable::make(Int(32), prefix + split.old_var + ".loop_max");
        Expr old_min = Variable::make(Int(32), prefix + split.old_var + ".loop_min");
        Expr old_extent = Variable::make(Int(32), prefix + split.old_var + ".loop_extent");

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
Stmt build_provide_loop_nest_helper(string func_name,
                                    string prefix,
                                    const vector<string> &dims,
                                    const vector<Expr> &site,
                                    const vector<Expr> &values,
                                    const vector<Expr> &predicates,
                                    const FuncSchedule &func_s,
                                    const StageSchedule &stage_s,
                                    bool is_update,
                                    bool compile_to_tiramisu) {

    // We'll build it from inside out, starting from a store node,
    // then wrapping it in for loops.

    // Make the (multi-dimensional multi-valued) store node.
    Stmt stmt = Provide::make(func_name, values, site); // TODO(tiramisu): computation expression (buffer mapping ???)

    vector<Split> splits = stage_s.splits();

    // Define the function args in terms of the loop variables using the splits
    if (!compile_to_tiramisu) {
        for (const Split &split : splits) {
            vector<ApplySplitResult> splits_result = apply_split(split, is_update, prefix, dim_extent_alignment);
        }
    }

    // All containing lets and fors. Outermost first.
    vector<Container> nest;

    // Put the desired loop nest into the containers vector.
    for (int i = (int)stage_s.dims().size() - 1; i >= 0; i--) {
        // TODO(tiramisu): for-loop dimensions (outermost to innermost) (-> after splits etc are applied)
        const Dim &dim = stage_s.dims()[i];
        Container c = {Container::For, i, prefix + dim.var, Expr()};
        nest.push_back(c);
    }

    // Put all the reduction domain predicates into the containers vector.
    for (Expr pred : predicates) {
        pred = qualify(prefix, pred);
        Container c = {Container::If, 0, "", likely(pred)};
        pred_container.push_back(c);
    }
    int n_predicates = pred_container.size();

    // Define the loop mins and extents in terms of the mins and maxs produced by bounds inference
    for (const std::string &i : dims) { // TODO(psuriana): pure loop dimension (no splits, etc)
        string var = prefix + i;
        Expr max = Variable::make(Int(32), var + ".max");
        Expr min = Variable::make(Int(32), var + ".min"); // Inject instance name here? (compute instance names during lowering)
        stmt = LetStmt::make(var + ".loop_extent",
                             (max + 1) - min,
                             stmt);
        stmt = LetStmt::make(var + ".loop_min", min, stmt);
        stmt = LetStmt::make(var + ".loop_max", max, stmt);
    }

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
Stmt build_provide_loop_nest(string func_name,
                             string prefix,
                             const vector<string> &dims,
                             const FuncSchedule &f_sched,
                             const Definition &def,
                             bool is_update,
                             bool compile_to_tiramisu) {

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
    Stmt stmt = build_provide_loop_nest_helper(
        func_name, prefix, dims, site, values, def.split_predicate(), f_sched,
        def.schedule(), is_update, compile_to_tiramisu);

    // Make any specialized copies
    const vector<Specialization> &specializations = def.specializations();
    for (size_t i = specializations.size(); i > 0; i--) {
        const Specialization &s = specializations[i-1];
        Expr c = s.condition;
        const Definition &s_def = s.definition;
        Stmt then_case;
        if (s.failure_message.empty()) {
            then_case = build_provide_loop_nest(func_name, prefix, dims, f_sched,
                                                s_def, is_update, compile_to_tiramisu);
        } else {
            internal_assert(equal(c, const_true()));
            // specialize_fail() should only be possible on the final specialization
            internal_assert(i == specializations.size());
            Expr specialize_fail_error =
                Internal::Call::make(Int(32),
                                     "halide_error_specialize_fail",
                                     {StringImm::make(s.failure_message)},
                                     Internal::Call::Extern);
            then_case = AssertStmt::make(const_false(), specialize_fail_error);
        }
        stmt = IfThenElse::make(c, then_case, stmt);
    }

    return stmt;
}


Stmt build_produce(Function f, const Target &t) {
    if (f.has_extern_definition()) {
        // TODO(tiramisu): Handle extern definition. Need to take care of
        // buffer allocation, annotation, etc.
        user_error << "Function \"" << f.name() << "\" is an extern function, "
                   << "which is not currently supported by Tiramisu.\n";
    } else {
        string prefix = f.name() + ".s0.";
        vector<string> dims = f.args();
        return build_provide_loop_nest(f.name(), prefix, dims, f.schedule(),
                                       f.definition(), false, compile_to_tiramisu);
    }
}

// Build the loop nests that update a function (assuming it's a reduction).
vector<Stmt> build_update(Function f, bool compile_to_tiramisu) {

    vector<Stmt> updates;

    for (size_t i = 0; i < f.updates().size(); i++) {
        const Definition &def = f.update(i);

        string prefix = f.name() + ".s" + std::to_string(i+1) + ".";

        vector<string> dims = f.args();
        Stmt loop = build_provide_loop_nest(f.name(), prefix, dims, f.schedule(),
                                            def, true, compile_to_tiramisu);
        updates.push_back(loop);
    }

    return updates;
}

pair<Stmt, Stmt> build_production(Function func, const Target &target, bool compile_to_tiramisu) {
    Stmt produce = build_produce(func, target, compile_to_tiramisu);
    vector<Stmt> updates = build_update(func, compile_to_tiramisu);

    // Combine the update steps
    Stmt merged_updates = Block::make(updates);
    return { produce, merged_updates };
}

// Inject the allocation and realization of a function into an
// existing loop nest using its schedule
class GetFuncSchedule : public IRVisitor {
public:
    const Function &func;
    bool is_output;
    const Target &target;

    GetFuncSchedule(const Function &f, bool o, const Target &t) :
        func(f), is_output(o), target(t), ss(s) {}

private:

    ostringstream ss;
    string producing;

    Stmt build_pipeline(Stmt consumer) { // tiramisu::computation (each definition will need its own computation)
        pair<Stmt, Stmt> realization = build_production(func, target, compile_to_tiramisu);

        Stmt producer;
        if (realization.first.defined() && realization.second.defined()) {
            producer = Block::make(realization.first, realization.second);
        } else if (realization.first.defined()) {
            producer = realization.first;
        } else {
            internal_assert(realization.second.defined());
            producer = realization.second;
        }
        producer = ProducerConsumer::make_produce(func.name(), producer);

        // Outputs don't have consume nodes
        if (!is_output) {
            consumer = ProducerConsumer::make_consume(func.name(), consumer);
        }

        if (is_no_op(consumer)) {
            // For the very first output to be scheduled, the consumer
            // Stmt will be a no-op. No point in preserving it.
            return producer;
        } else {
            return Block::make(producer, consumer);
        }
    }

    Stmt build_realize(Stmt s) { // tiramisu::buffer
        if (!is_output) {
            Region bounds;
            string name = func.name();
            const vector<string> func_args = func.args();
            for (int i = 0; i < func.dimensions(); i++) {
                const string &arg = func_args[i];
                Expr min = Variable::make(Int(32), name + "." + arg + ".min_realized");
                Expr extent = Variable::make(Int(32), name + "." + arg + ".extent_realized");
                bounds.push_back(Range(min, extent));
            }

            s = Realize::make(name, func.output_types(), bounds, const_true(), s);
        }

        // This is also the point at which we inject explicit bounds
        // for this realization.
        inject_explicit_bounds(s, func);
    }

    using IRVisitor::visit;

    void visit(const For *for_loop) {
        debug(3) << "InjectRealization of " << func.name() << " entering for loop over " << for_loop->name << "\n";
        const LoopLevel &compute_level = func.schedule().compute_level();
        const LoopLevel &store_level = func.schedule().store_level();

        Stmt body = for_loop->body;

        // Dig through any let statements
        vector<pair<string, Expr>> lets;
        while (const LetStmt *l = body.as<LetStmt>()) {
            lets.push_back({ l->name, l->value });
            body = l->body;
        }

        // Can't schedule extern things inside a vector for loop
        if (func.has_extern_definition() &&
            func.schedule().compute_level().is_inline() &&
            for_loop->for_type == ForType::Vectorized &&
            function_is_used_in_stmt(func, for_loop)) {

            // If we're trying to inline an extern function, schedule it here and bail out
            debug(2) << "Injecting realization of " << func.name() << " around node " << Stmt(for_loop) << "\n";
            stmt = build_realize(build_pipeline(for_loop));
            found_store_level = found_compute_level = true;
            return;
        }

        body = mutate(body);

        if (compute_level.match(for_loop->name)) {
            debug(3) << "Found compute level\n";
            if (function_is_used_in_stmt(func, body) || is_output) {
                body = build_pipeline(body);
            }
            found_compute_level = true;
        }

        if (store_level.match(for_loop->name)) {
            debug(3) << "Found store level\n";
            internal_assert(found_compute_level)
                << "The compute loop level was not found within the store loop level!\n";

            if (function_is_used_in_stmt(func, body) || is_output) {
                body = build_realize(body);
            }

            found_store_level = true;
        }

        // Reinstate the let statements
        for (size_t i = lets.size(); i > 0; i--) {
            body = LetStmt::make(lets[i - 1].first, lets[i - 1].second, body);
        }

        if (body.same_as(for_loop->body)) {
            stmt = for_loop;
        } else {
            stmt = For::make(for_loop->name,
                             for_loop->min,
                             for_loop->extent,
                             for_loop->for_type,
                             for_loop->device_api,
                             body);
        }
    }

    // If we're an inline update or extern, we may need to inject a realization here
    virtual void visit(const Provide *op) {
        if (op->name != func.name() &&
            !func.is_pure() &&
            func.schedule().compute_level().is_inline() &&
            function_is_used_in_stmt(func, op)) {

            // Prefix all calls to func in op
            stmt = build_realize(build_pipeline(op));
            found_store_level = found_compute_level = true;
        } else {
            stmt = op;
        }
    }
};


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
                                                        vector<string> &order,
                                                        map<string, Function> &env) {

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

        GetFuncSchedule func_sched(f, is_output, target);
        func_sched.accept(&func_sched);
    }

    return schedules;
}

}
}
