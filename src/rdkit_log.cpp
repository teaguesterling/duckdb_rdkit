#include "rdkit_log.hpp"
#include "duckdb/function/scalar_function.hpp"
#include <RDGeneral/RDLog.h>

namespace duckdb_rdkit {

// rdkit_log_disable() - Disable all RDKit logging
static void rdkit_log_disable(DataChunk &args, ExpressionState &state,
                              Vector &result) {
  auto count = args.size();

  // Disable all RDKit loggers
  boost::logging::disable_logs("rdApp.*");

  // Return true to indicate success
  auto result_data = FlatVector::GetData<bool>(result);
  for (idx_t i = 0; i < count; i++) {
    result_data[i] = true;
  }
}

// rdkit_log_enable() - Enable all RDKit logging
static void rdkit_log_enable(DataChunk &args, ExpressionState &state,
                             Vector &result) {
  auto count = args.size();

  // Enable all RDKit loggers
  boost::logging::enable_logs("rdApp.*");

  // Return true to indicate success
  auto result_data = FlatVector::GetData<bool>(result);
  for (idx_t i = 0; i < count; i++) {
    result_data[i] = true;
  }
}

// rdkit_log_status() - Get current logging status
static void rdkit_log_status(DataChunk &args, ExpressionState &state,
                             Vector &result) {
  auto count = args.size();

  std::string status = boost::logging::log_status();

  for (idx_t i = 0; i < count; i++) {
    result.SetValue(i, Value(status));
  }
}

void RegisterLogFunctions(ExtensionLoader &loader) {
  // rdkit_log_disable() - returns BOOLEAN
  ScalarFunctionSet disable_set("rdkit_log_disable");
  disable_set.AddFunction(
      ScalarFunction({}, LogicalType::BOOLEAN, rdkit_log_disable));
  loader.RegisterFunction(disable_set);

  // rdkit_log_enable() - returns BOOLEAN
  ScalarFunctionSet enable_set("rdkit_log_enable");
  enable_set.AddFunction(
      ScalarFunction({}, LogicalType::BOOLEAN, rdkit_log_enable));
  loader.RegisterFunction(enable_set);

  // rdkit_log_status() - returns VARCHAR with log status
  ScalarFunctionSet status_set("rdkit_log_status");
  status_set.AddFunction(
      ScalarFunction({}, LogicalType::VARCHAR, rdkit_log_status));
  loader.RegisterFunction(status_set);
}

} // namespace duckdb_rdkit
