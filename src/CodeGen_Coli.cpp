#include <iostream>
#include <limits>

#include "CodeGen_Coli.h"
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
    "#include <coli/debug.h>\n"
    "#include <coli/core.h>\n\n"
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

string CodeGen_Coli::print(Expr e) {
    internal_assert(e.defined()) << "CodeGen_Coli can't convert undefined expr.\n";
    // For now, substitute in all lets to make life easier (does not substitute in lets in stmt though)
    e = substitute_in_all_lets(e);
    e.accept(this);
    return expr;
}

void CodeGen_Coli::print(Stmt s) {
    internal_assert(s.defined()) << "CodeGen_Coli can't convert undefined stmt.\n";
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
        }
        else oss << name[i];
    }
    return oss.str();
}

string halide_type_to_coli_type_str(Type type) {
    if (type.is_uint()) {
        if (type.bits() == 8) {
            return "coli::p_uint8";
        } else if (type.bits() == 16) {
            return "coli::p_uint16";
        } else if (type.bits() == 32) {
            return "coli::p_uint32";
        } else {
            return "coli::p_uint64";
        }
    } else if (type.is_int()) {
        if (type.bits() == 8) {
            return "coli::p_int8";
        } else if (type.bits() == 16) {
            return "coli::p_int16";
        } else if (type.bits() == 32) {
            return "coli::p_int32";
        } else {
            return "coli::p_int64";
        }
    } else if (type.is_float()) {
        if (type.bits() == 32) {
            return "coli::p_float32";
        } else if (type.bits() == 64) {
            return "coli::p_float64";
        } else {
            user_error << "Floats other than 32 and 64 bits are not suppored in Coli.\n";
        }
    } else if (type.is_bool()) {
        return "coli::p_boolean";
    } else {
        user_error << "Halide type cannot be translated to Coli type.\n";
    }
    return "coli::p_none";
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
}

CodeGen_Coli::CodeGen_Coli(ostream &dest, const string &pipeline_name,
                           const vector<Function> &outputs,
                           const vector<vector<int32_t>> &output_buffer_extents,
                           const vector<Type> &output_buffer_types,
                           const vector<string> &inputs,
                           const vector<vector<int32_t>> &input_buffer_extents,
                           const vector<Type> &input_buffer_types,
                           const vector<string> &order,
                           const map<string, Function> &env)
        : stream(dest), indent(0), func(pipeline_name), order(order), env(env), loop_depth(0),
          current_computation("") {

    internal_assert(outputs.size() == output_buffer_extents.size());
    internal_assert(output_buffer_extents.size() == output_buffer_types.size());
    internal_assert(inputs.size() == input_buffer_extents.size());
    internal_assert(input_buffer_extents.size() == input_buffer_types.size());

    stream << headers << "\n\n";
    stream << "using namespace coli;\n\n";
    stream << "int main(int argc, char **argv)\n";
    stream << "{\n";

    indent += tab_size;

    stream << do_indent();
    stream << "// Set default coli options.\n";
    stream << do_indent();
    stream << "global::set_default_coli_options();\n\n";
    stream << do_indent();
    stream << "coli::function " << func << "(\"" << func << "\")" << ";\n";

    // Allocate the output buffers
    for (size_t k = 0; k < outputs.size(); ++k) {
        const Function &f = outputs[k];
        const vector<int32_t> &buffer_extents = output_buffer_extents[k];
        const Type type = output_buffer_types[k];

        internal_assert(buffer_extents.size() == f.args().size());

        ostringstream sizes;
        sizes << "{";
        for (size_t i = 0; i < buffer_extents.size(); ++i) {
            sizes << "coli::expr(" << buffer_extents[i] << ")";
            string min_str = print_name(f.name() + ".min." + std::to_string(i));
            string extent_str = print_name(f.name() + ".extent." + std::to_string(i));
            scope.push_back(std::make_pair(min_str, make_const(Int(32), 0)));
            scope.push_back(std::make_pair(extent_str, make_const(Int(32), buffer_extents[i])));
            if (i != buffer_extents.size() - 1) {
                sizes << ", ";
            }
        }
        sizes << "}";

        string buffer_name = "buff_" + print_name(f.name());
        stream << do_indent();
        stream << "coli::buffer " << buffer_name << "(\"" << buffer_name << "\", "
               << f.args().size() << ", " << sizes.str() << ", "
               << halide_type_to_coli_type_str(type) << ", NULL, coli::a_output, "
               << "&" << func << ");\n";
        output_buffers.insert(buffer_name);
    }

    // Bind to the input buffers
    for (size_t k = 0; k < inputs.size(); ++k) {
        const string &input_name = inputs[k];
        const vector<int32_t> &buffer_extents = input_buffer_extents[k];
        const Type type = input_buffer_types[k];

        vector<string> dummy_dims(buffer_extents.size());

        ostringstream sizes;
        sizes << "{";
        for (size_t i = 0; i < buffer_extents.size(); ++i) {
            dummy_dims[i] = "i" + std::to_string(i);
            push_loop_dim(dummy_dims[i], make_const(Int(32), 0), buffer_extents[i], input_name, 0, "");

            sizes << "coli::expr(" << buffer_extents[i] << ")";
            if (i != buffer_extents.size() - 1) {
                sizes << ", ";
            }
        }
        sizes << "}";

        string buffer_name = "buff_" + print_name(input_name);
        stream << do_indent();
        stream << "coli::buffer " << buffer_name << "(\"" << buffer_name << "\", "
               << buffer_extents.size() << ", " << sizes.str() << ", "
               << halide_type_to_coli_type_str(type) << ", NULL, coli::a_input, "
               << "&" << func << ");\n";
        input_buffers.insert(buffer_name);

        // Bind the input buffer to a computation
        string dims_str = to_string(dummy_dims);

        string symbolic_str = get_loop_bound_vars();
        string iter_space_str;
        if (!symbolic_str.empty()) {
            iter_space_str = get_loop_bound_vars() + "->{" + input_name + dims_str + ": " + get_loop_bounds() + "}";
        } else {
            iter_space_str = "{" + input_name + dims_str + ": " + get_loop_bounds() + "}";
        }

        stream << do_indent();
        stream << "coli::computation " << input_name << "(\"" << iter_space_str << "\", "
               << "expr(), false, " << halide_type_to_coli_type_str(type)
               << ", &" << func << ");\n";

        // 1-to-1 mapping to buffer
        string access_str = "{" + input_name + dims_str + "->" + "buff_" + input_name + dims_str + "}";
        stream << do_indent();
        stream << input_name << ".set_access(\"" << access_str << "\");\n";
        stream << "\n";

        computation_list.insert(input_name);

        for (size_t i = 0; i < buffer_extents.size(); ++i) {
            pop_loop_dim();
        }
    }
}

namespace {

int get_split_fuse_dim_index(const vector<Dim> &dims, const string &var) {
    // TODO(psuriana): we need to pass the *original* dim list (i.e. the one
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

void CodeGen_Coli::generate_schedule(const Function &func, int stage,
                                     const Schedule &schedule, size_t order_index) {
    string name = print_name(func.name() + ".s" + std::to_string(stage));

    const vector<Dim> &dims = schedule.dims();
    size_t dim_size = dims.size() - 1; // Ignore __outermost

    // Add split/fuse schedule
    // TODO(psuriana): the mapping from dim name to index is a complete mess in COLi. Not
    // sure how this will pan out will a more complex schedules (with reorder, etc)
    /*for (const Split &split : schedule.splits()) {
        if (split.is_split()) {
            int index = get_split_fuse_dim_index(dims, split.old_var);

            debug(0) << "Splitting " << func.name() << "." << split.old_var << "(" << index << ")\n";
            stream << do_indent();
            stream << func.name() << ".split(" << std::to_string(index) << ", "
                   << split.factor << ");\n";

        } else {
            user_error << "COLi only handles split.\n";
        }
    }*/

    // Tag any parallel dimensions
    for (size_t i = 0; i < dim_size; ++i) {
        const Dim &d = dims[i];
        if (d.for_type == ForType::Parallel) {
            debug(5) << "...Parallelize " << name << "." << d.var << "\n";
            stream << do_indent();
            stream << name << ".tag_parallel_dimension(" << std::to_string(dim_size - i) << ");\n";
        } else if (d.for_type == ForType::Vectorized) {
            internal_error << "Does not currently support vectorization\n";
            /*debug(5) << "...Vectorize " << name << "." << d.var << "\n";
            stream << do_indent();
            stream << name << ".tag_vector_dimension(" << std::to_string(dim_size - i) << ");\n";*/
        } else {
            internal_assert(d.for_type == ForType::Serial)
                << "Can only emit serial/parallel/vectorized for loops to COLi\n";
        }
    }

    // TODO(psuriana): add GPU schedules
}

void CodeGen_Coli::generate_schedules() {
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

CodeGen_Coli::~CodeGen_Coli() {
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
    stream << func << ".set_arguments(" << buffers_stream.str() << ");\n";
    stream << do_indent();
    stream << func << ".gen_time_processor_domain();\n";
    stream << do_indent();
    stream << func << ".gen_isl_ast();\n";
    stream << do_indent();
    stream << func << ".gen_halide_stmt();\n";
    stream << do_indent();
    stream << func << ".dump_halide_stmt();\n";
    stream << do_indent();
    stream << func << ".gen_halide_obj(\"build/generated_" << func << "_test.o\");\n";

    stream << "\n";
    stream << do_indent();
    stream << "return 0;\n";

    indent -= tab_size;

    stream << do_indent();
    stream << "}\n\n";
}

string CodeGen_Coli::do_indent() const {
    ostringstream ss;
    for (int i = 0; i < indent; i++) {
        ss << ' ';
    }
    return ss.str();
}

void CodeGen_Coli::push_loop_dim(const string &name, Expr min, Expr extent,
                                 const string &func, int stage, const string &var) {
    loop_dims.push_back({name, min, extent, func, stage, var});
}

void CodeGen_Coli::pop_loop_dim() {
    loop_dims.pop_back();
}

string CodeGen_Coli::get_current_func_name() const {
    internal_assert(!loop_dims.empty());
    return loop_dims[loop_dims.size()-1].func;
}

int CodeGen_Coli::get_current_stage() const {
    internal_assert(!loop_dims.empty());
    return loop_dims[loop_dims.size()-1].stage;
}

string CodeGen_Coli::get_loop_bound_vars() const {
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

string CodeGen_Coli::get_loop_bounds() const {
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

void CodeGen_Coli::visit(const StringImm *op) {
    user_error << "Conversion of StringImm to COLi is not supported.\n";
}

void CodeGen_Coli::visit(const AssertStmt *op) {
    user_error << "Conversion of AssertStmt to COLi is not supported.\n";
}

void CodeGen_Coli::visit(const Ramp *op) {
    user_error << "Conversion of Ramp to COLi is not supported.\n";
}

void CodeGen_Coli::visit(const Broadcast *op) {
    user_error << "Conversion of Broadcast to COLi is not supported.\n";
}

void CodeGen_Coli::visit(const IfThenElse *op) {
    user_error << "Conversion of IfThenElse to COLi is not supported.\n";
}

void CodeGen_Coli::visit(const Free *op) {
    user_error << "Conversion of Free to COLi is not supported.\n";
}

void CodeGen_Coli::visit(const Store *op) {
    user_error << "Should pass the unflatten version of Store to COLi\n.\n";
}

void CodeGen_Coli::visit(const Allocate *op) {
    user_error << "Should pass the unflatten version of Allocate to COLi\n.\n";
}

void CodeGen_Coli::visit(const IntImm *op) {
    ostringstream ss;
    ss << "coli::expr(";
    if (op->type.bits() == 8) {
        ss << "(int8_t)";
    } else if (op->type.bits() == 16) {
        ss << "(int16_t)";
    } else if (op->type.bits() == 32) {
        ss << "(int32_t)";
    }
    ss << op->value << ")";
    expr = ss.str();
}

void CodeGen_Coli::visit(const UIntImm *op) {
    ostringstream ss;
    ss << "coli::expr(";
    if (op->type.bits() == 8) {
        ss << "(uint8_t)";
    } else if (op->type.bits() == 16) {
        ss << "(uint16_t)";
    } else if (op->type.bits() == 32) {
        ss << "(uint32_t)";
    }
    ss << op->value << ")";
    expr = ss.str();
}

void CodeGen_Coli::visit(const FloatImm *op) {
    ostringstream ss;
    if (op->type.bits() == 32) {
        ss << "coli::expr((float)" << op->value << ")";
    } else if (op->type.bits() == 64) {
        ss << "coli::expr(" << op->value << ")";
    } else {
        // Only support 32- and 64-bit integer
        user_error << "Conversion of float " << op->type.bits() << "_t to COLi is not currently supported.\n";
    }
    expr = ss.str();
}

void CodeGen_Coli::visit(const Cast *op) {
    ostringstream ss;
    ss << "coli::expr(coli::o_cast, ";
    ss << halide_type_to_coli_type_str(op->type);
    ss << ", ";
    ss << print(op->value);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Coli::visit(const Variable *op) {
    user_assert(!op->param.defined() && !op->image.defined())
        << "Can only handle conversion of simple variable to COLi for now.\n";

    ostringstream ss;
    const auto &iter = constant_list.find(op->name);
    if (iter != constant_list.end()) {
        // It is a reference to variable defined in Let/LetStmt
        //TODO(psuriana): when do we actually generate constant???
        ss << (*iter) << "(0)";
    } else {
        // It is presumably a reference to loop variable
        ss << "coli::idx(\"" << op->name << "\")";
    }
    expr = ss.str();
}

void CodeGen_Coli::visit(const Add *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " + ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Sub *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " - ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Mul *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << "*";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Div *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << "/";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Mod *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " % ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Min *op) {
    ostringstream ss;
    ss << "coli::expr(coli::o_min, ";
    ss << print(op->a);
    ss << ", ";
    ss << print(op->b);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Coli::visit(const Max *op) {
    ostringstream ss;
    ss << "coli::expr(coli::o_max, ";
    ss << print(op->a);
    ss << ", ";
    ss << print(op->b);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Coli::visit(const EQ *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " == ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const NE *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " != ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const LT *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " < ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const LE *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " <= ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const GT *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " > ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const GE *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " >= ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const And *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " && ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Or *op) {
    ostringstream ss;
    ss << '(';
    ss << print(op->a);
    ss << " || ";
    ss << print(op->b);
    ss << ')';
    expr = ss.str();
}

void CodeGen_Coli::visit(const Not *op) {
    ostringstream ss;
    ss << '!';
    ss << print(op->a);
    expr = ss.str();
}

void CodeGen_Coli::visit(const Select *op) {
    ostringstream ss;
    ss << "coli::expr(coli::o_cond, ";
    ss << print(op->condition);
    ss << ", ";
    ss << print(op->true_value);
    ss << ", ";
    ss << print(op->false_value);
    ss << ")";
    expr = ss.str();
}

void CodeGen_Coli::visit(const Let *op) {
    user_error << "Should not have encountered Let expr since we've called substitute_in_all_lets.\n";
}

void CodeGen_Coli::visit(const LetStmt *op) {
    scope.push_back(std::make_pair(op->name, op->value));
    print(op->body);
    scope.pop_back();
}

void CodeGen_Coli::visit(const ProducerConsumer *op) {
    internal_assert(!op->is_producer || (computation_list.find(op->name) == computation_list.end()))
        << "Found another computation with the same name.\n";

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
        Function func = env.find(op->name)->second;

        const LoopLevel &compute_at = func.schedule().compute_level();
        if (compute_at.is_root() || compute_at.is_inline()) {
            if (order[0] == op->name) {
                // Initial definition
                stream << do_indent();
                stream << print_name(op->name + ".s0") << ".first(computation::root_dimension);\n";
                // Update definitions
                for (size_t i = 0; i < func.updates().size(); ++i) {
                    stream << do_indent();
                    stream << print_name(op->name + ".s" + std::to_string(i+1)) << ".first(computation::root_dimension);\n";
                }
            } else {
                const auto iter = std::find_if(order.begin(), order.end(),
                    [&op](const string& f) { return (f == op->name); });
                internal_assert(iter != order.end());
                int order_index = iter - order.begin();

                // Initial definition
                stream << do_indent();
                stream << print_name(op->name + ".s0") << ".after(" << print_name(order[order_index-1] + "_s")
                       << std::to_string(func.updates().size()) << ", computation::root_dimension);\n";
                // Update definitions
                for (size_t i = 0; i < func.updates().size(); ++i) {
                    stream << do_indent();
                    stream << print_name(op->name + ".s" + std::to_string(i+1)) << ".after("
                           << print_name(order[order_index-1] + "_s") << ", computation::root_dimension);\n";
                }
            }
        } else {
            string enclosing_func_name = get_current_func_name();
            internal_assert(compute_at.func().name() == enclosing_func_name);
            int enclosing_func_stage = get_current_stage();
            string parent = print_name(enclosing_func_name + "_s" + std::to_string(enclosing_func_stage));

            // TODO(psuriana): should compute at stage as well
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

string CodeGen_Coli::define_constant(const string &name, Expr val) {
    internal_assert(constant_list.find(name) == constant_list.end())
        << "Redefinition of lets \"" << name << "\" is not supported right now.\n";

    ostringstream ss;

    val = simplify(val);

    ss << "coli::constant " << name << "(\"" << name << "\", ";
    ss << print(val);
    ss << ", " << halide_type_to_coli_type_str(val.type())
       << ", true, NULL, 0, &" << func << ");\n";

    constant_list.insert(name);

    return ss.str();
}

string CodeGen_Coli::define_wrapper_let(const string &computation_name,
                                        const string &name, Expr val) {
    internal_assert(constant_list.find(name) == constant_list.end())
        << "Redefinition of lets \"" << name << "\" is not supported right now.\n";
    internal_assert(!current_computation.empty())
        << "The computation name should not be empty.\n";

    ostringstream ss;

    val = simplify(val);

    ss << "coli::constant " << name << "(\"" << name << "\", ";
    ss << print(val);
    ss << ", " << halide_type_to_coli_type_str(val.type())
       << ", false, &" << computation_name << ", 1, &" << func << ");\n";

    constant_list.insert(name);

    return ss.str();
}

void CodeGen_Coli::visit(const For *op) {
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

void CodeGen_Coli::visit(const Evaluate *op) {
    //TODO(psuriana): do nothing for now
}

void CodeGen_Coli::visit(const Load *op) {
    user_error << "Conversion of Load to COLi is not currently supported.\n";
}

void CodeGen_Coli::visit(const Provide *op) {
    string name = print_name(op->name + "_s" + std::to_string(get_current_stage()));
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
    ss << "coli::computation " << name << "(\"";
    indent += 5*tab_size;

    string dims_str = to_string(op->args);
    string symbolic_str = get_loop_bound_vars();
    if (!symbolic_str.empty()) {
        ss << get_loop_bound_vars() + "->{" << name + dims_str << ": \"\n";
    } else {
        ss << "{" << name << dims_str + ": \"\n";
    }

    ss << do_indent();
    ss << "\"" << get_loop_bounds() << "}\", \n";
    ss << do_indent();
    ss << print(op->values[0]);
    ss << ", true, " << halide_type_to_coli_type_str(op->values[0].type())
           << ", &" << func << ");\n";
    indent -= 5*tab_size;

    stream << ss.str();
    computation_list.insert(name);

    if (!buffer_str.empty()) {
        for (const auto &str : buffer_str) {
            stream << do_indent();
            stream << str;
        }
        buffer_str.clear();
    }

    // 1-to-1 mapping to buffer
    string access_str = "{" + name + dims_str + "->" + buffer_name + dims_str + "}";
    stream << do_indent();
    stream << name << ".set_access(\"" << access_str << "\");\n";

    current_computation = old_computation;
}

Expr CodeGen_Coli::substitute_in_scope(Expr expr) const {
    for (int i = scope.size() - 1; i >= 0; --i) {
        expr = substitute(scope[i].first, scope[i].second, expr);
    }
    return simplify(expr);
}

void CodeGen_Coli::generate_buffer(const Realize *op) {
    string name = print_name(op->name);

    user_assert(temporary_buffers.find("buff_" + name) == temporary_buffers.end())
        << "Duplicate allocation (i.e. duplicate compute) is not currently supported.\n";

    // Assert that the types of all buffer dimensions are the same for now.
    for (size_t i = 1; i < op->types.size(); ++i) {
        user_assert(op->types[i-1] == op->types[i])
            << "Realize node should have the same types for all dimensions for now.\n";
    }

    vector<Range> bounds(op->bounds);

    // Substitute it in all references to some other variables in the bounds
    for (Range &r : bounds) {
        r.min = substitute_in_scope(r.min);
        r.extent = substitute_in_scope(r.extent);
    }

    // Assert that the bounds on the dimensions start from 0 for now.
    for (size_t i = 0; i < bounds.size(); ++i) {
        user_assert(is_zero(bounds[i].min))
            << "Bound of realize node \"" << name << "\" should start from 0 for now.\n"
            << "Got " << bounds[i].min << " instead.\n"
            << "Original bound min: " << op->bounds[i].min << "\n";
    }

    // Create a temporary buffer

    string buffer_name = "buff_" + name;
    stream << do_indent();
    stream << "coli::buffer " << buffer_name << "(\"" << buffer_name << "\", "
           << bounds.size() << ", ";

    stream << "{";
    for (size_t i = 0; i < bounds.size(); ++i) {
        stream << print(bounds[i].extent);
        if (i != bounds.size() - 1) {
            stream << ", ";
        }
    }
    stream << "}, ";

    stream << halide_type_to_coli_type_str(op->types[0]) << ", NULL, coli::a_temporary, "
           << "&" << func << ");\n";

    temporary_buffers.insert(buffer_name);
}

void CodeGen_Coli::visit(const Realize *op) {
    // We will ignore the condition on the Realize node for now.

    stream << "\n";
    stream << do_indent();
    stream << "// Define temporary buffers for \"" << op->name << "\".\n";

    internal_assert(env.count(op->name));
    Function func = env.find(op->name)->second;

    generate_buffer(op);

    print(op->body);
}

void CodeGen_Coli::visit(const Call *op) {
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
        ss << "coli::expr(o_floor, ";
        ss << print(op->args[0]);
        ss << ')';
    } else {
        user_assert((op->call_type == Call::CallType::Halide) || (op->call_type == Call::CallType::Image))
            << "Only handle call to halide func or image for now.\n"
            << Expr(op) << "\n"
            << "is pure? " << op->is_pure() << "\n";

        //TODO(psuriana): need to normalize the index
        vector<Expr> normalized_args(op->args);
        for (auto &arg : normalized_args) {
            // Normalize it by introduction let Expr if it is not a Sub or Add or Variable
            if (!(arg.as<Variable>() || is_const(arg) || arg.as<Add>() || arg.as<Sub>())) {
                string var_name = unique_name('t');
                Expr var = Variable::make(Int(32), var_name);
                buffer_str.push_back(define_wrapper_let(current_computation, var_name, substitute_in_scope(arg)));
                arg = var;
            }
        }

        string call_name;
        if (op->call_type == Call::CallType::Halide) {
            internal_assert(env.count(op->name));
            Function func = env.find(op->name)->second;
            call_name = print_name(op->name + ".s" + std::to_string(func.updates().size()));
            internal_assert(computation_list.find(call_name) != computation_list.end())
                << "Call to computation \"" << call_name << "\" that does not exist.\n";
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

void CodeGen_Coli::visit(const Block *op) {
    print(op->first);
    if (op->rest.defined()) print(op->rest);
}

void CodeGen_Coli::test() {

    std::cout << "CodeGen_Coli test passed\n";
}

void print_to_coli(Stmt s, ostream &dest, const string &pipeline_name,
                   const vector<Function> &outputs,
                   const vector<vector<int32_t>> &output_buffer_extents,
                   const vector<Type> &output_buffer_types,
                   const vector<string> &inputs,
                   const vector<vector<int32_t>> &input_buffer_extents,
                   const vector<Type> &input_buffer_types,
                   const vector<string> &order,
                   const map<string, Function> &env) {

    NormalizeVariableName normalize;
    s = normalize.mutate(s);
    debug(0) << "After normalization:\n" << s << "\n\n";

    // TODO(psuriana): Need to re-normalize the allocation bound to start from 0

    CodeGen_Coli cg(dest, pipeline_name, outputs, output_buffer_extents,
                    output_buffer_types, inputs, input_buffer_extents,
                    input_buffer_types, order, env);
    cg.print(s);
}

}
}
