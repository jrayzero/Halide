#include <iostream>
#include <limits>

#include "CodeGen_Tiramisu.h"
#include "IROperator.h"
#include "IRMutator.h"
#include "Substitute.h"
#include "Func.h"

namespace Halide {
namespace Internal {

using std::ostream;
using std::set;
using std::string;
using std::vector;
using std::ostringstream;
using std::map;
using std::pair;

namespace {

const string headers =
    "#include <isl/set.h>\n"
    "#include <isl/union_map.h>\n"
    "#include <isl/union_set.h>\n"
    "#include <isl/ast_build.h>\n"
    "#include <isl/schedule.h>\n"
    "#include <isl/schedule_node.h>\n\n"
    "#include <tiramisu/debug.h>\n"
    "#include <tiramisu/core.h>\n\n"
    "#include <string.h>\n"
    "#include <Halide.h>\n"
    "#include \"halide_image_io.h\"\n";

const int tab_size = 4;
}

template<typename T>
std::string to_string(const std::vector<T>& v) {
    std::ostringstream ss;
    ss << "[";
    for (size_t i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i != v.size() - 1) {
            ss << ", ";
        }
    }
    ss << "]";
    return ss.str();
}

string CodeGen_Tiramisu::print(Expr e) {
    internal_assert(e.defined()) << "CodeGen_Tiramisu can't convert undefined expr.\n";
    // For now, substitute in all lets to make life easier (does not substitute in lets in stmt though)
    e = substitute_in_all_lets(e);
    e.accept(this);
    return expr;
}

void CodeGen_Tiramisu::print(Stmt s) {
    internal_assert(s.defined()) << "CodeGen_Tiramisu can't convert undefined stmt.\n";
    // For now, substitute in all lets to make life easier (does not substitute in lets in stmt though)
    s = substitute_in_all_lets(s);
    s.accept(this);
}

namespace {
string print_name(const string &name) {
    ostringstream oss;

    for (size_t i = 0; i < name.size(); i++) {
        if (name[i] == '.') {
            oss << '_';
        } else if (name[i] == '$') {
            oss << "__";
        } else if (name[i] != '_' && !isalnum(name[i])) {
            oss << "___";
        } else {
            oss << name[i];
        }
    }
    return oss.str();
}

string halide_type_to_tiramisu_type_str(Type type) {
    if (type.is_uint()) {
        if (type.bits() == 8) {
            return "tiramisu::p_uint8";
        } else if (type.bits() == 16) {
            return "tiramisu::p_uint16";
        } else if (type.bits() == 32) {
            return "tiramisu::p_uint32";
        } else {
            internal_assert(type.bits() == 64);
            return "tiramisu::p_uint64";
        }
    } else if (type.is_int()) {
        if (type.bits() == 8) {
            return "tiramisu::p_int8";
        } else if (type.bits() == 16) {
            return "tiramisu::p_int16";
        } else if (type.bits() == 32) {
            return "tiramisu::p_int32";
        } else {
            internal_assert(type.bits() == 64);
            return "tiramisu::p_int64";
        }
    } else if (type.is_float()) {
        if (type.bits() == 32) {
            return "tiramisu::p_float32";
        } else if (type.bits() == 64) {
            return "tiramisu::p_float64";
        } else {
            user_error << "Floats other than 32 and 64 bits are not suppored in Tiramisu.\n";
        }
    } else if (type.is_bool()) {
        return "tiramisu::p_boolean";
    } else {
        user_error << "Halide type \"" << type << "\" cannot be translated to Tiramisu type.\n";
    }
    return "tiramisu::p_none";
}

string halide_type_to_c_type_str(Type type) {
    if (type.is_uint()) {
        if (type.bits() == 8) {
            return "uint8_t";
        } else if (type.bits() == 16) {
            return "uint16_t";
        } else if (type.bits() == 32) {
            return "uint32_t";
        } else {
            internal_assert(type.bits() == 64);
            return "uint64_t";
        }
    } else if (type.is_int()) {
        if (type.bits() == 8) {
            return "int8_t";
        } else if (type.bits() == 16) {
            return "int16_t";
        } else if (type.bits() == 32) {
            return "int32_t";
        } else {
            internal_assert(type.bits() == 64);
            return "int64_t";
        }
    } else if (type.is_float()) {
        if (type.bits() == 32) {
            return "float";
        } else if (type.bits() == 64) {
            return "double";
        } else {
            user_error << "Floats other than 32 and 64 bits are not suppored in C.\n";
        }
    } else if (type.is_bool()) {
        return "bool";
    } else {
        user_error << "Halide type \"" << type << "\" cannot be translated to C type.\n";
    }
    return "void";
}

class NormalizeVariableName : public IRMutator {
    using IRMutator::visit;

    void visit(const For *op) {
        string name = print_name(op->name);
        Expr min = mutate(op->min);
        Expr extent = mutate(op->extent);
        Stmt body = mutate(op->body);
        if ((name == op->name) &&
            min.same_as(op->min) &&
            extent.same_as(op->extent) &&
            body.same_as(op->body)) {
            stmt = op;
        } else {
            stmt = For::make(name, min, extent, op->for_type, op->device_api, body);
        }
    }

    void visit(const Let *op) {
        string name = print_name(op->name);
        Expr value = mutate(op->value);
        Expr body = mutate(op->body);
        if ((name == op->name) &&
            value.same_as(op->value) &&
            body.same_as(op->body)) {
            expr = op;
        } else {
            expr = Let::make(name, value, body);
        }
    }

    void visit(const LetStmt *op) {
        string name = print_name(op->name);
        Expr value = mutate(op->value);
        Stmt body = mutate(op->body);
        if ((name == op->name) &&
            value.same_as(op->value) &&
            body.same_as(op->body)) {
            stmt = op;
        } else {
            stmt = LetStmt::make(name, value, body);
        }
    }

    void visit(const Variable *op) {
        string name = print_name(op->name);
        if (name != op->name) {
            expr = Variable::make(op->type, name, op->image, op->param, op->reduction_domain);
        } else {
            expr = op;
        }
    }
};

} // anonymous namespace

CodeGen_Tiramisu::CodeGen_Tiramisu(ostream &dest, const string &pipeline_name,
                                   const vector<Function> &outputs,
                                   const vector<vector<int32_t>> &output_buffer_extents,
                                   const vector<Type> &output_buffer_types,
                                   const vector<string> &inputs,
                                   const vector<vector<int32_t>> &input_buffer_extents,
                                   const vector<Type> &input_buffer_types,
                                   const vector<string> &input_params,
                                   const vector<Type> &input_param_types,
                                   const vector<string> &order,
                                   const map<string, Function> &env)
        : stream(dest), indent(0), pipeline(pipeline_name), order(order), env(env), loop_depth(0),
          current_computation("") {

    internal_assert(outputs.size() == output_buffer_extents.size());
    internal_assert(output_buffer_extents.size() == output_buffer_types.size());
    internal_assert(inputs.size() == input_buffer_extents.size());
    internal_assert(input_buffer_extents.size() == input_buffer_types.size());

    stream << headers << "\n\n";
    stream << "using namespace tiramisu;\n\n";
    stream << "int main(int argc, char **argv)\n";
    stream << "{\n";

    indent += tab_size;

    stream << do_indent();
    stream << "// Set default tiramisu options.\n";
    stream << do_indent();
    stream << "global::set_default_tiramisu_options();\n\n";
    stream << do_indent();
    stream << "tiramisu::function " << pipeline << "(\"" << pipeline << "\")" << ";\n";

    // Define the input params
    if (!input_params.empty()) {
        stream << "\n";
        stream << do_indent();
        stream << "// Input params.\n";
        for (size_t k = 0; k < input_params.size(); ++k) {
            const string &param_name = input_params[k];
            const Type &type = input_param_types[k];

            string extent_str = print_name(param_name);
            extent_list.insert(extent_str);
            stream << do_indent();
            // Just assigned it to some random value. We need to replace it
            // manually later to the intended value.
            stream << halide_type_to_c_type_str(type) << " " << extent_str << " = 0;\n";
        }
    }

    // Allocate the output buffers
    stream << "\n";
    stream << do_indent();
    stream << "// Output buffers.\n";
    for (size_t k = 0; k < outputs.size(); ++k) {
        const Function &f = outputs[k];
        const vector<int32_t> &buffer_extents = output_buffer_extents[k];
        const Type &type = output_buffer_types[k];

        internal_assert(buffer_extents.size() == f.args().size());

        ostringstream sizes;
        sizes << "{";
        for (int i = buffer_extents.size()-1; i >= 0; --i) {
            string min_str = print_name(f.name() + ".min." + std::to_string(i));

            string extent_str = print_name(f.name() + ".extent." + std::to_string(i));
            extent_list.insert(extent_str);
            stream << do_indent();
            stream << "int " << extent_str << " = " << buffer_extents[i] << ";\n";

            sizes << "tiramisu::expr(" << extent_str << ")";
            // TODO(psuriana): Assume "min" always starts from 0 for now
            scope.push_back(std::make_pair(min_str, make_const(Int(32), 0)));
            //scope.push_back(std::make_pair(extent_str, make_const(Int(32), buffer_extents[i])));
            if (i != 0) {
                sizes << ", ";
            }
        }
        sizes << "}";

        string buffer_name = "buff_" + print_name(f.name());
        stream << do_indent();
        stream << "tiramisu::buffer " << buffer_name << "(\"" << buffer_name << "\", "
               << f.args().size() << ", " << sizes.str() << ", "
               << halide_type_to_tiramisu_type_str(type) << ", NULL, tiramisu::a_output, "
               << "&" << pipeline << ");\n";
        output_buffers.insert(buffer_name);
    }

    // Declare the input buffers
    stream << "\n";
    stream << do_indent();
    stream << "// Input buffers.\n";
    for (size_t k = 0; k < inputs.size(); ++k) {
        const string &input_name = inputs[k];
        const vector<int32_t> &buffer_extents = input_buffer_extents[k];
        const Type &type = input_buffer_types[k];

        vector<string> dummy_dims;

        ostringstream sizes;
        sizes << "{";
        for (int i = buffer_extents.size()-1; i >= 0; --i) {
            string dummy_dim = "i" + std::to_string(i);
            dummy_dims.push_back(dummy_dim);

            string extent_str = print_name(input_name + ".extent." + std::to_string(i));
            extent_list.insert(extent_str);
            stream << do_indent();
            stream << "int " << extent_str << " = " << buffer_extents[i] << ";\n";

            push_loop_dim(dummy_dim, make_const(Int(32), 0), Variable::make(Int(32), extent_str), input_name, 0, "");

            sizes << "tiramisu::expr(" << extent_str << ")";
            if (i != 0) {
                sizes << ", ";
            }
        }
        sizes << "}";

        string buffer_name = "buff_" + print_name(input_name);
        stream << do_indent();
        stream << "tiramisu::buffer " << buffer_name << "(\"" << buffer_name << "\", "
               << buffer_extents.size() << ", " << sizes.str() << ", "
               << halide_type_to_tiramisu_type_str(type) << ", NULL, tiramisu::a_input, "
               << "&" << pipeline << ");\n";
        input_buffers.insert(buffer_name);

        // Bind the input buffer to a computation
        string dims_str = to_string(dummy_dims);

        string symbolic_str = get_loop_bound_vars();
        string iter_space_str;
        if (!symbolic_str.empty()) {
            iter_space_str = symbolic_str + "->{" + input_name + dims_str + ": " + get_loop_bounds() + "}";
        } else {
            iter_space_str = "{" + input_name + dims_str + ": " + get_loop_bounds() + "}";
        }

        stream << do_indent();
        stream << "tiramisu::computation " << input_name << "(\"" << iter_space_str << "\", "
               << "expr(), false, " << halide_type_to_tiramisu_type_str(type)
               << ", &" << pipeline << ");\n";

        // 1-to-1 mapping to buffer
        string access_str = "{" + input_name + dims_str + "->" + "buff_" + input_name + dims_str + "}";
        stream << do_indent();
        stream << input_name << ".set_access(\"" << access_str << "\");\n";
        stream << "\n";

        computation_list.insert(input_name);

        for (int i = buffer_extents.size()-1; i >= 0; --i) {
            pop_loop_dim();
        }
    }
}

namespace {

int get_split_fuse_dim_index(const vector<Dim> &dims, const string &var) {
    // TODO(psuriana): we need to pass the "original" dim list (i.e. the one
    // before Halide applies any splits, reordering, etc.)
    return 0;
}

int get_compute_at_dim_index(const vector<Dim> &dims, const string &var) {
    const auto iter = std::find_if(dims.begin(), dims.end(),
        [&var](const Dim& d) { return (d.var == var); });
    internal_assert(iter != dims.end());
    return (iter - dims.begin());
}

} // anonymous namespace

void CodeGen_Tiramisu::generate_schedule(const Function &func, int stage,
                                         const StageSchedule &schedule, size_t order_index) {
    string name = print_name(func.name() + ".s" + std::to_string(stage));

    const vector<Dim> &dims = schedule.dims();
    size_t dim_size = dims.size() - 1; // Ignore __outermost

    // Add split/fuse schedule
    // TODO(psuriana): the mapping from dim name to index is a complete mess in Tiramisu. Not
    // sure how this will pan out will a more complex schedules (with reorder, etc)
    /*for (const Split &split : schedule.splits()) {
        if (split.is_split()) {
            int index = get_split_fuse_dim_index(dims, split.old_var);

            debug(0) << "Splitting " << func.name() << "." << split.old_var << "(" << index << ")\n";
            stream << do_indent();
            stream << func.name() << ".split(" << std::to_string(index) << ", "
                   << split.factor << ");\n";

        } else {
            user_error << "Tiramisu only handles split.\n";
        }
    }*/

    // Tag any parallel dimensions
    for (size_t i = 0; i < dim_size; ++i) {
        const Dim &d = dims[i];
        if (d.for_type == ForType::Parallel) {
            debug(5) << "...Parallelize " << name << "." << d.var << " at index " << i << "\n";
            stream << do_indent();
            stream << name << ".tag_parallel_dimension(" << std::to_string(dim_size - i - 1) << ");\n";
        } else if (d.for_type == ForType::Vectorized) {
            internal_error << "Does not currently support vectorization\n";
            /*debug(5) << "...Vectorize " << name << "." << d.var << "\n";
            stream << do_indent();
            stream << name << ".tag_vector_dimension(" << std::to_string(dim_size - i) << ");\n";*/
        } else {
            internal_assert(d.for_type == ForType::Serial)
                << "Can only emit serial/parallel/vectorized for loops to Tiramisu\n";
        }
    }

    // TODO(psuriana): add GPU schedules
}

void CodeGen_Tiramisu::generate_schedules() {
    internal_assert(!order.empty());
    stream << "\n";
    stream << do_indent();
    stream << "// Add schedules.\n";
    for (size_t i = 0; i < order.size(); ++i) {
        const string &func_name = order[i];
        internal_assert(env.count(func_name));
        const Function &func = env.find(func_name)->second;

        // TODO(psuriana): schedule the update definition as well
        generate_schedule(func, 0, func.definition().schedule(), i);
    }
}

CodeGen_Tiramisu::~CodeGen_Tiramisu() {
    generate_schedules();

    // Bind the output and input buffers
    ostringstream buffers_stream;
    buffers_stream << "{";
    int count = 0;
    // Input buffers
    for (auto it = input_buffers.begin(); it != input_buffers.end(); ++it) {
        count += 1;
        buffers_stream << "&" << (*it);
        if (it != (--input_buffers.end())) {
            buffers_stream << ", ";
        }
    }
    // Output buffers
    bool start = true;
    for (auto it = output_buffers.begin(); it != output_buffers.end(); ++it) {
        if ((count > 0) && start) {
            start = false;
            buffers_stream << ", ";
        }
        buffers_stream << "&" << (*it);
        if (it != (--output_buffers.end())) {
            buffers_stream << ", ";
        }
    }
    buffers_stream << "}";

    stream << "\n";
    stream << do_indent();
    stream << pipeline << ".set_arguments(" << buffers_stream.str() << ");\n";
    stream << do_indent();
    stream << pipeline << ".gen_time_processor_domain();\n";
    stream << do_indent();
    stream << pipeline << ".gen_isl_ast();\n";
    stream << do_indent();
    stream << pipeline << ".gen_halide_stmt();\n";
    stream << do_indent();
    stream << pipeline << ".dump_halide_stmt();\n";
    stream << do_indent();
    stream << pipeline << ".gen_halide_obj(\"build/generated_fct_" << pipeline << ".o\");\n";

    stream << "\n";
    stream << do_indent();
    stream << "return 0;\n";

    indent -= tab_size;

    stream << do_indent();
    stream << "}\n\n";
}

string CodeGen_Tiramisu::do_indent() const {
    ostringstream ss;
    for (int i = 0; i < indent; i++) {
        ss << ' ';
    }
    return ss.str();
}

void CodeGen_Tiramisu::push_loop_dim(const string &name, Expr min, Expr extent,
                                     const string &func, int stage, const string &var) {
    loop_dims.push_back({name, min, extent, func, stage, var});
}

void CodeGen_Tiramisu::pop_loop_dim() {
    loop_dims.pop_back();
}

string CodeGen_Tiramisu::get_current_func_name() const {
    internal_assert(!loop_dims.empty());
    return loop_dims[loop_dims.size()-1].func;
}

int CodeGen_Tiramisu::get_current_stage() const {
    internal_assert(!loop_dims.empty());
    return loop_dims[loop_dims.size()-1].stage;
}

vector<string> CodeGen_Tiramisu::get_stage_dims(const string &name, int stage, bool ignore_rvar) const {
    internal_assert(env.count(name));
    const Function &func = env.find(name)->second;
    const StageSchedule &schedule = (stage == 0) ? func.definition().schedule() : func.update(stage-1).schedule();
    const vector<Dim> &dims = schedule.dims();

    vector<string> dim_str;
    // Ignore __outermost
    for (int i = dims.size()-2; i >= 0; --i) {
        if (ignore_rvar && dims[i].is_rvar()) {
            continue;
        }
        dim_str.push_back(print_name(name + ".s" + std::to_string(stage) + "." + dims[i].var));
    }
    return dim_str;
}

vector<string> CodeGen_Tiramisu::get_stage_rvars(const string &name, int stage) const {
    internal_assert(env.count(name));
    const Function &func = env.find(name)->second;
    const StageSchedule &schedule = (stage == 0) ? func.definition().schedule() : func.update(stage-1).schedule();
    const vector<Dim> &dims = schedule.dims();

    vector<string> dim_str;
    // Ignore __outermost
    for (int i = dims.size()-2; i >= 0; --i) {
        if (dims[i].is_rvar()) {
            dim_str.push_back(print_name(name + ".s" + std::to_string(stage) + "." + dims[i].var));
        }
    }
    return dim_str;
}

// Return str representation of all symbolic loop bounds
string CodeGen_Tiramisu::get_loop_bound_vars() const {
    vector<Expr> relevant_exprs;
    for (size_t i = 0; i < loop_dims.size(); ++i) {
        if (!is_const(loop_dims[i].min)) {
            relevant_exprs.push_back(loop_dims[i].min);
        }
        if (!is_const(loop_dims[i].extent)) {
            relevant_exprs.push_back(loop_dims[i].extent);
        }
    }
    if (relevant_exprs.empty()) {
        return "";
    }
    return to_string(relevant_exprs);
}

string CodeGen_Tiramisu::get_loop_bounds() const {
    ostringstream ss;
    ss << "(";
    for (size_t i = 0; i < loop_dims.size(); ++i) {
        ss << loop_dims[i].to_string();
        if (i != loop_dims.size() - 1) {
            ss << ") and (";
        }
    }
    ss << ")";
    return ss.str();
}

void CodeGen_Tiramisu::visit(const StringImm *op) {
    user_error << "Conversion of StringImm to Tiramisu is not supported.\n";
}

void CodeGen_Tiramisu::visit(const AssertStmt *op) {
    user_error << "Conversion of AssertStmt to Tiramisu is not supported.\n";
}

void CodeGen_Tiramisu::visit(const Ramp *op) {
    user_error << "Conversion of Ramp to Tiramisu is not supported.\n";
}

void CodeGen_Tiramisu::visit(const Broadcast *op) {
    user_error << "Conversion of Broadcast to Tiramisu is not supported.\n";
}

void CodeGen_Tiramisu::visit(const IfThenElse *op) {
    user_error << "Conversion of IfThenElse to Tiramisu is not supported.\n";
}

void CodeGen_Tiramisu::visit(const Free *op) {
    user_error << "Conversion of Free to Tiramisu is not supported.\n";
}

void CodeGen_Tiramisu::visit(const Shuffle *op) {
    user_error << "Conversion of Shuffle to Tiramisu is not currently supported.\n";
}

void CodeGen_Tiramisu::visit(const Prefetch *op) {
    user_error << "Conversion of Prefetch to Tiramisu is not currently supported.\n";
}

void CodeGen_Tiramisu::visit(const Load *op) {
    user_error << "Should pass the unflatten version of Load (in the form of Call node) to Tiramisu\n.\n";
}

void CodeGen_Tiramisu::visit(const Store *op) {
    user_error << "Should pass the unflatten version of Store (in the form of Provide node) to Tiramisu\n.\n";
}

void CodeGen_Tiramisu::visit(const Allocate *op) {
    user_error << "Should pass the unflatten version of Allocate (in the form of Realize node) to Tiramisu\n.\n";
}

void CodeGen_Tiramisu::visit(const IntImm *op) {
    ostringstream ss;
    ss << "tiramisu::expr(";
    if (op->type.bits() == 8) {
        ss << "(int8_t)";
    } else if (op->type.bits() == 16) {
        ss << "(int16_t)";
    } else if (op->type.bits() == 32) {
        ss << "(int32_t)";
    } else {
        internal_assert(op->type.bits() == 64);
        ss << "(int64_t)";
    }
    ss << op->value << ")";
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const UIntImm *op) {
    ostringstream ss;
    ss << "tiramisu::expr(";
    if (op->type.bits() == 8) {
        ss << "(uint8_t)";
    } else if (op->type.bits() == 16) {
        ss << "(uint16_t)";
    } else if (op->type.bits() == 32) {
        ss << "(uint32_t)";
    } else {
        internal_assert(op->type.bits() == 64);
        ss << "(uint32_t)";
    }
    ss << op->value << ")";
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const FloatImm *op) {
    ostringstream ss;
    if (op->type.bits() == 32) {
        ss << "tiramisu::expr((float)" << op->value << ")";
    } else if (op->type.bits() == 64) {
        ss << "tiramisu::expr(" << op->value << ")";
    } else {
        // Only support 32- and 64-bit integer
        user_error << "Conversion of " << op->type << " to Tiramisu is not currently supported.\n";
    }
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Cast *op) {
    ostringstream ss;
    ss << "tiramisu::expr(tiramisu::o_cast, ";
    ss << halide_type_to_tiramisu_type_str(op->type);
    ss << ", ";
    ss << print(op->value);
    ss << ")";
    expr = ss.str();
}

//TODO(psuriana) : fix this
void CodeGen_Tiramisu::visit(const Variable *op) {
    /*user_assert(!op->param.defined() && !op->image.defined())
        << "Can only handle conversion of simple variable to Tiramisu for now.\n";*/

    ostringstream ss;
    const auto &iter = constant_list.find(op->name);
    if (iter != constant_list.end()) {
        // It is a reference to variable defined in Let/LetStmt
        //TODO(psuriana): when do we actually generate constant???
        ss << (*iter) << "(0)";
    } else {
        const auto &it = extent_list.find(op->name);
        if (it != extent_list.end()) {
            ss << "tiramisu::expr(" << op->name << ")";
        } else {
            // It is presumably a reference to loop variable
            ss << "tiramisu::idx(\"" << op->name << "\")";
        }
    }
    expr = ss.str();
}

template <typename T>
void CodeGen_Tiramisu::visit_binary(const T *op, const string &op_str) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " " << op_str << " ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Add *op) { visit_binary(op, "+"); }
void CodeGen_Tiramisu::visit(const Sub *op) { visit_binary(op, "-"); }
void CodeGen_Tiramisu::visit(const Mul *op) { visit_binary(op, "*"); }
void CodeGen_Tiramisu::visit(const Div *op) { visit_binary(op, "/"); }
void CodeGen_Tiramisu::visit(const Mod *op) { visit_binary(op, "%"); }
void CodeGen_Tiramisu::visit(const EQ *op)  { visit_binary(op, "=="); }
void CodeGen_Tiramisu::visit(const NE *op)  { visit_binary(op, "!="); }
void CodeGen_Tiramisu::visit(const LT *op)  { visit_binary(op, "<"); }
void CodeGen_Tiramisu::visit(const LE *op)  { visit_binary(op, "<="); }
void CodeGen_Tiramisu::visit(const GT *op)  { visit_binary(op, ">"); }
void CodeGen_Tiramisu::visit(const GE *op)  { visit_binary(op, ">="); }
void CodeGen_Tiramisu::visit(const And *op) { visit_binary(op, "&&"); }
void CodeGen_Tiramisu::visit(const Or *op)  { visit_binary(op, "||"); }

void CodeGen_Tiramisu::visit(const Min *op) {
    ostringstream ss;
    ss << "tiramisu::expr(tiramisu::o_min, ";
    ss << print(op->a);
    ss << ", ";
    ss << print(op->b);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Max *op) {
    ostringstream ss;
    ss << "tiramisu::expr(tiramisu::o_max, ";
    ss << print(op->a);
    ss << ", ";
    ss << print(op->b);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Not *op) {
    ostringstream ss;
    ss << '!';
    ss << print(op->a);
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Select *op) {
    ostringstream ss;
    ss << "tiramisu::expr(tiramisu::o_select, ";
    ss << print(op->condition);
    ss << ", ";
    ss << print(op->true_value);
    ss << ", ";
    ss << print(op->false_value);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Let *op) {
    user_error << "Should not have encountered Let expr since we've called substitute_in_all_lets.\n";
}

void CodeGen_Tiramisu::visit(const LetStmt *op) {
    scope.push_back(std::make_pair(op->name, op->value));
    print(op->body);
    scope.pop_back();
}

void CodeGen_Tiramisu::visit(const ProducerConsumer *op) {
    internal_assert(!op->is_producer || (computation_list.find(op->name) == computation_list.end()))
        << "Found another computation with the same name \"" << op->name << "\".\n";

    vector<Loop> old_loop_dims = loop_dims;
    int old_loop_depth = loop_depth;
    loop_depth = 0;
    print(op->body);
    loop_dims = old_loop_dims;
    loop_depth = old_loop_depth;

    if (op->is_producer) {
        stream << "\n";
        stream << do_indent();
        stream << "// Define compute level for \"" << op->name << "\".\n";

        internal_assert(env.count(op->name));
        const Function &func = env.find(op->name)->second;

        const LoopLevel &compute_at = func.schedule().compute_level();
        // TODO(psuriana): Do you need to do .after() for inlined?
        if (compute_at.is_root() || compute_at.is_inline()) {
            if (order[0] == op->name) {
                // Initial definition
                stream << do_indent();
                stream << print_name(op->name + ".s0") << ".first(computation::root_dimension);\n";
                // Update definitions
                for (size_t i = 0; i < func.updates().size(); ++i) {
                    stream << do_indent();
                    stream << print_name(op->name + ".s" + std::to_string(i+1)) << ".after("
                           << print_name(op->name + ".s" + std::to_string(i)) << ", computation::root_dimension);\n";
                }
            } else {
                const auto iter = std::find_if(order.begin(), order.end(),
                    [&op](const string& f) { return (f == op->name); });
                internal_assert(iter != order.end());
                int order_index = iter - order.begin();

                internal_assert(env.count(order[order_index-1]));
                const Function &prev = env.find(order[order_index-1])->second;

                // Initial definition
                stream << do_indent();
                stream << print_name(op->name + ".s0") << ".after(" << print_name(prev.name() + ".s")
                       << std::to_string(prev.updates().size()) << ", computation::root_dimension);\n";
                // Update definitions
                for (size_t i = 0; i < func.updates().size(); ++i) {
                    stream << do_indent();
                    stream << print_name(op->name + ".s" + std::to_string(i+1)) << ".after("
                           << print_name(op->name + ".s" + std::to_string(i)) << ", computation::root_dimension);\n";
                }
            }
        } else {
            string enclosing_func_name = get_current_func_name();
            internal_assert(compute_at.func() == enclosing_func_name);
            int enclosing_func_stage = get_current_stage();
            string parent = print_name(enclosing_func_name + ".s" + std::to_string(enclosing_func_stage));

            debug(5) << "Compute Func " << op->name << " at " << compute_at.var().name() << "(" << loop_depth << ")\n";

            // Initial definition
            stream << do_indent();
            stream << print_name(op->name + ".s0") << ".after(" << parent << ", " << std::to_string(loop_depth) << ");\n";
            // Update definitions
            for (size_t i = 0; i < func.updates().size(); ++i) {
                stream << do_indent();
                stream << print_name(op->name + ".s" + std::to_string(i+1)) << ".after(" << parent
                       << ", " << std::to_string(loop_depth) << ");\n";
            }
        }
    }
}

string CodeGen_Tiramisu::define_constant(const string &name, Expr val) {
    internal_assert(constant_list.find(name) == constant_list.end())
        << "Redefinition of lets \"" << name << "\" is not supported right now.\n";

    ostringstream ss;

    val = simplify(val);

    // Global constant
    ss << "tiramisu::constant " << name << "(\"" << name << "\", ";
    ss << print(val);
    ss << ", " << halide_type_to_tiramisu_type_str(val.type())
       << ", true, NULL, 0, &" << pipeline << ");\n";

    constant_list.insert(name);

    return ss.str();
}

string CodeGen_Tiramisu::define_wrapper_let(const string &computation_name,
                                        const string &name, Expr val) {
    internal_assert(constant_list.find(name) == constant_list.end())
        << "Redefinition of lets \"" << name << "\" is not supported right now.\n";
    internal_assert(!current_computation.empty())
        << "The computation name should not be empty.\n";

    ostringstream ss;

    val = simplify(val);

    ss << "tiramisu::constant " << name << "(\"" << name << "\", ";
    ss << print(val);
    ss << ", " << halide_type_to_tiramisu_type_str(val.type())
       << ", false, &" << computation_name << ", 1, &" << pipeline << ");\n";

    constant_list.insert(name);

    return ss.str();
}

void CodeGen_Tiramisu::visit(const For *op) {
    loop_depth += 1;

    vector<string> v = split_string(op->name, "_");
    internal_assert((v.size() > 2) && (!v[0].empty()) && (!v[v.size()-1].empty()));
    string func_name = v[0];
    string var = v[v.size()-1];

    int producer_stage = -1;
    for (size_t i = 1; i < v.size() - 1; ++i) {
        if (v[i].substr(0, 1) == "s") {
            string str = v[i].substr(1, v[i].size() - 1);
            bool has_only_digits = (str.find_first_not_of( "0123456789" ) == string::npos);
            if (has_only_digits) {
                producer_stage = atoi(str.c_str());
            }
        }
    }
    internal_assert(producer_stage >= 0);

    push_loop_dim(op->name, op->min, op->extent, func_name, producer_stage, var);

    const Variable *min = op->min.as<Variable>();
    internal_assert(min != NULL) << "Min value of a loop should have been a variable.\n";
    const Variable *extent = op->extent.as<Variable>();
    internal_assert(extent != NULL) << "Extent of a loop should have been a variable.\n";

    Expr min_val;
    for (int i = scope.size() - 1; i >= 0; --i) {
        if (scope[i].first == min->name) {
            min_val = scope[i].second;
            break;
        }
    }
    internal_assert(min_val.defined()) << "min->name: " << min->name << "\n";

    Expr extent_val;
    for (int i = scope.size() - 1; i >= 0; --i) {
        if (scope[i].first == extent->name) {
            extent_val = scope[i].second;
            break;
        }
    }
    internal_assert(extent_val.defined());

    // Substitute it in all references to some other variables in the min/extent val
    min_val = substitute_in_scope(min_val);
    extent_val = substitute_in_scope(extent_val);

    stream << "\n";
    stream << do_indent();
    stream << "// Define loop bounds for dimension \"" << op->name << "\".\n";
    stream << do_indent();
    stream << define_constant(min->name, min_val);
    stream << do_indent();
    stream << define_constant(extent->name, extent_val);

    print(op->body);
    pop_loop_dim();

    loop_depth -= 1;
}

void CodeGen_Tiramisu::visit(const Evaluate *op) {
    // TODO(psuriana): do nothing for now
}

void CodeGen_Tiramisu::visit(const Provide *op) {
    int stage = get_current_stage();
    string name = print_name(op->name + "_s" + std::to_string(stage));
    string buffer_name = "buff_" + print_name(op->name);

    string old_computation = current_computation;
    current_computation = name;

    internal_assert(computation_list.find(name) == computation_list.end())
        << "Duplicate computation \"" << op->name << "\" is not currently supported.\n";
    internal_assert(temporary_buffers.count(buffer_name) || output_buffers.count(buffer_name))
        << "The buffer should have been allocated previously.\n";

    for (size_t i = 0; i < op->args.size(); ++i) {
        user_assert(op->args[i].as<Variable>() != NULL)
            << "Expect args of provide to be loop dims for now (doesn't currently handle update).\n";
    }
    user_assert(op->values.size() == 1) << "Expect 1D store (no tuple) in the Provide node for now.\n";

    ostringstream ss;

    ss << do_indent();
    ss << "tiramisu::computation " << name << "(\n";
    indent += tab_size;

    string symbolic_str = get_loop_bound_vars();
    string dim_str = to_string(get_stage_dims(op->name, stage, false));
    ss << do_indent() << "\"";
    if (!symbolic_str.empty()) {
        ss << symbolic_str + "->{" << name + dim_str << ": \"\n";
    } else {
        ss << "{" << name << dim_str + ": \"\n";
    }

    ss << do_indent();
    ss << "\"" << get_loop_bounds() << "}\", \n";
    ss << do_indent();
    if (stage == 0) {
        ss << print(op->values[0]);
    } else {
        ss << "tiramisu::expr()";
    }

    ss << ", true, " << halide_type_to_tiramisu_type_str(op->values[0].type())
           << ", &" << pipeline << ");\n";
    indent -= tab_size;

    computation_list.insert(name);

    if (stage > 0) {
        ss << do_indent();
        ss << name << ".set_expression(" << print(op->values[0]) << ");\n";
    }

    stream << ss.str();

    if (!buffer_str.empty()) {
        for (const auto &str : buffer_str) {
            stream << do_indent();
            stream << str;
        }
        buffer_str.clear();
    }

    // 1-to-1 mapping to buffer
    string access_str = "{" + name + dim_str + "->" + buffer_name +
                        to_string(get_stage_dims(op->name, stage, true)) + "}";
    stream << do_indent();
    stream << name << ".set_access(\"" << access_str << "\");\n";

    current_computation = old_computation;
}

Expr CodeGen_Tiramisu::substitute_in_scope(Expr expr) const {
    for (int i = scope.size() - 1; i >= 0; --i) {
        expr = substitute(scope[i].first, scope[i].second, expr);
    }
    return simplify(expr);
}

void CodeGen_Tiramisu::generate_buffer(const Realize *op) {
    string name = print_name(op->name);

    // TODO(psuriana): Assume we don't have duplicate compute for now
    user_assert(temporary_buffers.find("buff_" + name) == temporary_buffers.end())
        << "Duplicate allocation (i.e. duplicate compute) is not currently supported.\n";

    // Assert that the types of all buffer dimensions are the same for now.
    for (size_t i = 1; i < op->types.size(); ++i) {
        user_assert(op->types[i-1] == op->types[i])
            << "Realize node should have the same types for all dimensions for now.\n";
    }

    vector<Range> bounds;

    // Substitute it in all references to some other variables in the bounds
    // Need to reverse it since Tiramisu is from outermost to innermost
    for (int i = op->bounds.size() - 1; i >= 0; --i) {
        Range r = op->bounds[i];
        r.min = substitute_in_scope(r.min);
        r.extent = substitute_in_scope(r.extent);
        bounds.push_back(r);
    }

    // TODO(psuriana): Assert that the bounds on the dimensions start from 0 for now.
    for (size_t i = 0; i < bounds.size(); ++i) {
        user_assert(is_zero(bounds[i].min))
            << "Bound of realize node \"" << name << "\" should start from 0 for now.\n"
            << "Got " << bounds[i].min << " instead.\n"
            << "Original bound min: " << op->bounds[i].min << "\n";
    }

    // Create a temporary buffer

    string buffer_name = "buff_" + name;
    stream << do_indent();
    stream << "tiramisu::buffer " << buffer_name << "(\"" << buffer_name << "\", "
           << bounds.size() << ", ";

    stream << "{";
    for (size_t i = 0; i < bounds.size(); ++i) {
        stream << print(bounds[i].extent);
        if (i != bounds.size() - 1) {
            stream << ", ";
        }
    }
    stream << "}, ";

    stream << halide_type_to_tiramisu_type_str(op->types[0]) << ", NULL, tiramisu::a_temporary, "
           << "&" << pipeline << ");\n";

    temporary_buffers.insert(buffer_name);
}

void CodeGen_Tiramisu::visit(const Realize *op) {
    // TODO(psuriana): We will ignore the condition on the Realize node for now.

    stream << "\n";
    stream << do_indent();
    stream << "// Define temporary buffers for \"" << print_name(op->name) << "\".\n";

    generate_buffer(op);

    print(op->body);
}

void CodeGen_Tiramisu::visit(const Call *op) {
    ostringstream ss;
    if (op->is_intrinsic(Call::shift_right)) {
        internal_assert(op->args.size() == 2);
        ss << '(';
        ss << print(op->args[0]);
        ss << " >> ";
        ss << print(op->args[1]);
        ss << ')';
    } else if (op->is_intrinsic(Call::shift_left)) {
        internal_assert(op->args.size() == 2);
        ss << '(';
        ss << print(op->args[0]);
        ss << " << ";
        ss << print(op->args[1]);
        ss << ')';
    } else if ((op->name == "floor_f16") || (op->name == "floor_f32") || (op->name == "floor_f64")) {
        internal_assert(op->args.size() == 1);
        ss << "tiramisu::expr(o_floor, ";
        ss << print(op->args[0]);
        ss << ')';
    } else {
        user_assert((op->call_type == Call::CallType::Halide) || (op->call_type == Call::CallType::Image))
            << "Only handle call to halide func or image for now.\n"
            << Expr(op) << "\n"
            << "is pure? " << op->is_pure() << "\n";

        //TODO(psuriana): need to normalize the index
        vector<Expr> normalized_args;
        for (int i = op->args.size()-1; i >= 0; --i) {
            Expr arg = op->args[i];
            // Normalize it by introducing let Expr if it is not a Sub or Add or Variable
            if (!(arg.as<Variable>() || is_const(arg) || arg.as<Add>() || arg.as<Sub>())) {
                string var_name = unique_name('t');
                Expr var = Variable::make(Int(32), var_name);
                buffer_str.push_back(define_wrapper_let(current_computation, var_name, substitute_in_scope(arg)));
                arg = var;
            }
            normalized_args.push_back(arg);
        }

        string call_name;
        if (op->call_type == Call::CallType::Halide) {
            string current_producer = get_current_func_name();
            int current_stage = get_current_stage();

            internal_assert(env.count(op->name));
            Function func = env.find(op->name)->second;

            if (op->name != current_producer) {
                call_name = print_name(op->name + ".s0");
            } else {
                call_name = print_name(op->name + ".s" + std::to_string(current_stage));
            }

            internal_assert(computation_list.find(call_name) != computation_list.end())
                << "Call to computation \"" << call_name << "\" that does not exist.\n";

            if ((current_producer == op->name) && (current_stage > 0)) {
                // It's a reduction, we need to add the "new" reduction dimension
                vector<string> rvars = get_stage_rvars(current_producer, current_stage);
                for (const auto &var_name : rvars) {
                    Expr var = Variable::make(Int(32), var_name);
                    normalized_args.push_back(var - 1);
                }
            }

        } else {
            call_name = print_name(op->name);
        }

        ss << call_name << "(";
        for (size_t i = 0; i < normalized_args.size(); i++) {
            ss << print(normalized_args[i]);
            if (i < normalized_args.size() - 1) {
                ss << ", ";
            }
        }
        ss << ")";
    }
    expr = ss.str();
}

void CodeGen_Tiramisu::visit(const Block *op) {
    print(op->first);
    if (op->rest.defined()) {
        print(op->rest);
    }
}

void CodeGen_Tiramisu::test() {

    std::cout << "CodeGen_Tiramisu test passed\n";
}

void print_to_tiramisu(Stmt s, ostream &dest, const string &pipeline_name,
                       const vector<Function> &outputs,
                       const vector<vector<int32_t>> &output_buffer_extents,
                       const vector<Type> &output_buffer_types,
                       const vector<string> &inputs,
                       const vector<vector<int32_t>> &input_buffer_extents,
                       const vector<Type> &input_buffer_types,
                       const vector<string> &input_params,
                       const vector<Type> &input_param_types,
                       const vector<string> &order,
                       const map<string, Function> &env) {

    NormalizeVariableName normalize;
    s = normalize.mutate(s);
    debug(3) << "After normalization:\n" << s << "\n\n";

    // TODO(psuriana): Need to re-normalize the allocation bound to start from 0

    CodeGen_Tiramisu cg(dest, pipeline_name, outputs, output_buffer_extents,
                        output_buffer_types, inputs, input_buffer_extents,
                        input_buffer_types, input_params, input_param_types,
                        order, env);
    cg.print(s);
}

}
}
