#include "mol_descriptors.hpp"
#include "common.hpp"
#include "duckdb/common/assert.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/function_set.hpp"
#include "mol_formats.hpp"
#include "qed.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"

namespace duckdb_rdkit {

// ============================================================================
// Helper template for descriptor functions
// ============================================================================

// Template for float descriptors from pure Mol (RDKit pickle)
template <typename DescriptorFn>
void mol_descriptor_float(DataChunk &args, ExpressionState &state, Vector &result,
                          DescriptorFn descriptor_fn) {
  D_ASSERT(args.data.size() == 1);
  auto &mol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      mol_vec, result, count, [&](string_t pickle) {
        auto mol = rdkit_binary_mol_to_mol(pickle.GetString());
        return descriptor_fn(*mol);
      });
}

// Template for float descriptors from UmbraMol
template <typename DescriptorFn>
void umbramol_descriptor_float(DataChunk &args, ExpressionState &state, Vector &result,
                               DescriptorFn descriptor_fn) {
  D_ASSERT(args.data.size() == 1);
  auto &umbramol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, float>(
      umbramol_vec, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto mol = rdkit_binary_mol_to_mol(umbra_mol.GetBinaryMol());
        return descriptor_fn(*mol);
      });
}

// Template for int32 descriptors from pure Mol
template <typename DescriptorFn>
void mol_descriptor_int(DataChunk &args, ExpressionState &state, Vector &result,
                        DescriptorFn descriptor_fn) {
  D_ASSERT(args.data.size() == 1);
  auto &mol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, int32_t>(
      mol_vec, result, count, [&](string_t pickle) {
        auto mol = rdkit_binary_mol_to_mol(pickle.GetString());
        return descriptor_fn(*mol);
      });
}

// Template for int32 descriptors from UmbraMol
template <typename DescriptorFn>
void umbramol_descriptor_int(DataChunk &args, ExpressionState &state, Vector &result,
                             DescriptorFn descriptor_fn) {
  D_ASSERT(args.data.size() == 1);
  auto &umbramol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, int32_t>(
      umbramol_vec, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto mol = rdkit_binary_mol_to_mol(umbra_mol.GetBinaryMol());
        return descriptor_fn(*mol);
      });
}

// ============================================================================
// Descriptor function implementations - Mol variants
// ============================================================================

void mol_logp_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    double logp, _;
    RDKit::Descriptors::calcCrippenDescriptors(mol, logp, _);
    return static_cast<float>(logp);
  });
}

void mol_qed_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    auto qed = QED();
    return static_cast<float>(qed.CalcQED(mol));
  });
}

void mol_amw_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<float>(RDKit::Descriptors::calcAMW(mol));
  });
}

void mol_exactmw_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<float>(RDKit::Descriptors::calcExactMW(mol));
  });
}

void mol_tpsa_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<float>(RDKit::Descriptors::calcTPSA(mol));
  });
}

void mol_hbd_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_int(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<int32_t>(RDKit::Descriptors::calcNumHBD(mol));
  });
}

void mol_hba_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_int(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<int32_t>(RDKit::Descriptors::calcNumHBA(mol));
  });
}

void mol_num_rotatable_bonds_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  mol_descriptor_int(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<int32_t>(RDKit::Descriptors::calcNumRotatableBonds(mol));
  });
}

// ============================================================================
// Descriptor function implementations - UmbraMol variants
// ============================================================================

void mol_logp_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    double logp, _;
    RDKit::Descriptors::calcCrippenDescriptors(mol, logp, _);
    return static_cast<float>(logp);
  });
}

void mol_qed_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    auto qed = QED();
    return static_cast<float>(qed.CalcQED(mol));
  });
}

void mol_amw_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<float>(RDKit::Descriptors::calcAMW(mol));
  });
}

void mol_exactmw_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<float>(RDKit::Descriptors::calcExactMW(mol));
  });
}

void mol_tpsa_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_float(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<float>(RDKit::Descriptors::calcTPSA(mol));
  });
}

void mol_hbd_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_int(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<int32_t>(RDKit::Descriptors::calcNumHBD(mol));
  });
}

void mol_hba_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_int(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<int32_t>(RDKit::Descriptors::calcNumHBA(mol));
  });
}

void mol_num_rotatable_bonds_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  umbramol_descriptor_int(args, state, result, [](const RDKit::ROMol &mol) {
    return static_cast<int32_t>(RDKit::Descriptors::calcNumRotatableBonds(mol));
  });
}

void RegisterDescriptorFunctions(ExtensionLoader &loader) {
  // mol_amw - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_amw("mol_amw");
  set_mol_amw.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_amw_mol));
  set_mol_amw.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::FLOAT, mol_amw_umbramol));
  loader.RegisterFunction(set_mol_amw);

  // mol_exactmw - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_exactmw("mol_exactmw");
  set_mol_exactmw.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_exactmw_mol));
  set_mol_exactmw.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::FLOAT, mol_exactmw_umbramol));
  loader.RegisterFunction(set_mol_exactmw);

  // mol_tpsa - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_tpsa("mol_tpsa");
  set_mol_tpsa.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_tpsa_mol));
  set_mol_tpsa.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::FLOAT, mol_tpsa_umbramol));
  loader.RegisterFunction(set_mol_tpsa);

  // mol_qed - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_qed("mol_qed");
  set_mol_qed.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_qed_mol));
  set_mol_qed.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::FLOAT, mol_qed_umbramol));
  loader.RegisterFunction(set_mol_qed);

  // mol_logp - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_logp("mol_logp");
  set_mol_logp.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::FLOAT, mol_logp_mol));
  set_mol_logp.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::FLOAT, mol_logp_umbramol));
  loader.RegisterFunction(set_mol_logp);

  // mol_hbd - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_hbd("mol_hbd");
  set_mol_hbd.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::INTEGER, mol_hbd_mol));
  set_mol_hbd.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::INTEGER, mol_hbd_umbramol));
  loader.RegisterFunction(set_mol_hbd);

  // mol_hba - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_hba("mol_hba");
  set_mol_hba.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::INTEGER, mol_hba_mol));
  set_mol_hba.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::INTEGER, mol_hba_umbramol));
  loader.RegisterFunction(set_mol_hba);

  // mol_num_rotatable_bonds - Mol and UmbraMol variants
  ScalarFunctionSet set_mol_num_rotatable_bonds("mol_num_rotatable_bonds");
  set_mol_num_rotatable_bonds.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol()}, LogicalType::INTEGER, mol_num_rotatable_bonds_mol));
  set_mol_num_rotatable_bonds.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol()}, LogicalType::INTEGER, mol_num_rotatable_bonds_umbramol));
  loader.RegisterFunction(set_mol_num_rotatable_bonds);
}
} // namespace duckdb_rdkit
