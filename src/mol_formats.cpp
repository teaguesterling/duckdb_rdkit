#include "mol_formats.hpp"
#include "common.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/function_set.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

namespace duckdb_rdkit {
// Expects a SMILES string and returns a RDKit pickled molecule
std::unique_ptr<RDKit::ROMol> rdkit_mol_from_smiles(std::string s) {
  std::string smiles = s;
  std::unique_ptr<RDKit::ROMol> mol;
  try {
    mol.reset(RDKit::SmilesToMol(smiles));
  } catch (std::exception &e) {
    std::string msg = StringUtil::Format("%s", typeid(e).name());
    throw InvalidInputException(msg);
  }

  if (mol) {
    return mol;
  } else {
    std::string msg = StringUtil::Format("Could not convert %s to mol", smiles);
    throw InvalidInputException(msg);
  }
}

// Serialize a molecule to binary using RDKit's MolPickler
std::string rdkit_mol_to_binary_mol(const RDKit::ROMol mol) {
  std::string buf;
  try {
    RDKit::MolPickler::pickleMol(mol, buf);
  } catch (...) {
    std::string msg = "Could not serialize mol to binary";
    throw InvalidInputException(msg);
  }
  return buf;
}

// Deserialize a binary mol to RDKit mol
std::unique_ptr<RDKit::ROMol> rdkit_binary_mol_to_mol(std::string bmol) {
  std::unique_ptr<RDKit::ROMol> mol(new RDKit::ROMol());
  RDKit::MolPickler::molFromPickle(bmol, *mol);

  return mol;
}

std::string rdkit_mol_to_smiles(RDKit::ROMol mol) {
  std::string smiles = RDKit::MolToSmiles(mol);
  return smiles;
}

// ============================================================================
// mol_to_smiles - convert to SMILES string
// ============================================================================

// mol_to_smiles for pure Mol (RDKit pickle)
void mol_to_smiles_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &mol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      mol_vec, result, count, [&](string_t pickle) {
        auto mol = rdkit_binary_mol_to_mol(pickle.GetString());
        auto smiles = rdkit_mol_to_smiles(*mol);
        return StringVector::AddString(result, smiles);
      });
}

// mol_to_smiles for UmbraMol
void mol_to_smiles_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &umbramol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      umbramol_vec, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto mol = rdkit_binary_mol_to_mol(umbra_mol.GetBinaryMol());
        auto smiles = rdkit_mol_to_smiles(*mol);
        return StringVector::AddString(result, smiles);
      });
}

// ============================================================================
// mol_from_smiles - parse SMILES to Mol (pure RDKit pickle)
// ============================================================================

void mol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          // Return pure RDKit pickle for Mol type
          auto pickle = rdkit_mol_to_binary_mol(*mol);
          return StringVector::AddStringOrBlob(result, pickle);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

// ============================================================================
// umbramol_from_smiles - parse SMILES to UmbraMol ([8B DalkeFP][RDKit Pickle])
// ============================================================================

void umbramol_from_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      smiles, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          // Return UmbraMol format with DalkeFP prefix
          auto umbramol = get_umbra_mol_string(*mol);
          return StringVector::AddStringOrBlob(result, umbramol);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

// ============================================================================
// Conversion functions between Mol and UmbraMol
// ============================================================================

// mol_to_umbramol: Convert pure Mol to UmbraMol format
void mol_to_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &mol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      mol_vec, result, count,
      [&](string_t pickle, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_binary_mol_to_mol(pickle.GetString());
          auto umbramol = get_umbra_mol_string(*mol);
          return StringVector::AddStringOrBlob(result, umbramol);
        } catch (...) {
          mask.SetInvalid(idx);
          return string_t();
        }
      });
}

// umbramol_to_mol: Extract pure RDKit pickle from UmbraMol
void umbramol_to_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &umbramol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      umbramol_vec, result, count,
      [&](string_t b_umbra_mol, ValidityMask &mask, idx_t idx) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        return StringVector::AddStringOrBlob(result, umbra_mol.GetBinaryMol());
      });
}

// ============================================================================
// mol_to_smarts - convert to SMARTS string
// ============================================================================

void mol_to_smarts_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &mol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      mol_vec, result, count, [&](string_t pickle) {
        auto mol = rdkit_binary_mol_to_mol(pickle.GetString());
        auto smarts = RDKit::MolToSmarts(*mol);
        return StringVector::AddString(result, smarts);
      });
}

void mol_to_smarts_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &umbramol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, string_t>(
      umbramol_vec, result, count, [&](string_t b_umbra_mol) {
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto mol = rdkit_binary_mol_to_mol(umbra_mol.GetBinaryMol());
        auto smarts = RDKit::MolToSmarts(*mol);
        return StringVector::AddString(result, smarts);
      });
}

// ============================================================================
// is_valid_smiles - check if string is a valid SMILES
// ============================================================================

void is_valid_smiles(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smiles_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, bool>(
      smiles_vec, result, count, [&](string_t smiles) {
        try {
          std::unique_ptr<RDKit::ROMol> mol(RDKit::SmilesToMol(smiles.GetString()));
          return mol != nullptr;
        } catch (...) {
          return false;
        }
      });
}

// ============================================================================
// is_valid_smarts - check if string is a valid SMARTS
// ============================================================================

void is_valid_smarts(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &smarts_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::Execute<string_t, bool>(
      smarts_vec, result, count, [&](string_t smarts) {
        try {
          std::unique_ptr<RDKit::ROMol> mol(RDKit::SmartsToMol(smarts.GetString()));
          return mol != nullptr;
        } catch (...) {
          return false;
        }
      });
}

void RegisterFormatFunctions(ExtensionLoader &loader) {
  // mol_from_smiles: SMILES -> Mol (pure RDKit pickle)
  ScalarFunctionSet mol_from_smiles_set("mol_from_smiles");
  mol_from_smiles_set.AddFunction(
      ScalarFunction({LogicalType::VARCHAR}, Mol(), mol_from_smiles));
  loader.RegisterFunction(mol_from_smiles_set);

  // umbramol_from_smiles: SMILES -> UmbraMol
  ScalarFunctionSet umbramol_from_smiles_set("umbramol_from_smiles");
  umbramol_from_smiles_set.AddFunction(
      ScalarFunction({LogicalType::VARCHAR}, UmbraMol(), umbramol_from_smiles));
  loader.RegisterFunction(umbramol_from_smiles_set);

  // mol_to_smiles: Mol/UmbraMol -> SMILES
  ScalarFunctionSet mol_to_smiles_set("mol_to_smiles");
  mol_to_smiles_set.AddFunction(
      ScalarFunction({Mol()}, LogicalType::VARCHAR, mol_to_smiles_mol));
  mol_to_smiles_set.AddFunction(
      ScalarFunction({UmbraMol()}, LogicalType::VARCHAR, mol_to_smiles_umbramol));
  loader.RegisterFunction(mol_to_smiles_set);

  // mol_to_umbramol: Mol -> UmbraMol
  ScalarFunctionSet mol_to_umbramol_set("mol_to_umbramol");
  mol_to_umbramol_set.AddFunction(
      ScalarFunction({Mol()}, UmbraMol(), mol_to_umbramol));
  loader.RegisterFunction(mol_to_umbramol_set);

  // umbramol_to_mol: UmbraMol -> Mol
  ScalarFunctionSet umbramol_to_mol_set("umbramol_to_mol");
  umbramol_to_mol_set.AddFunction(
      ScalarFunction({UmbraMol()}, Mol(), umbramol_to_mol));
  loader.RegisterFunction(umbramol_to_mol_set);

  // mol_to_smarts: Mol/UmbraMol -> SMARTS
  ScalarFunctionSet mol_to_smarts_set("mol_to_smarts");
  mol_to_smarts_set.AddFunction(
      ScalarFunction({Mol()}, LogicalType::VARCHAR, mol_to_smarts_mol));
  mol_to_smarts_set.AddFunction(
      ScalarFunction({UmbraMol()}, LogicalType::VARCHAR, mol_to_smarts_umbramol));
  loader.RegisterFunction(mol_to_smarts_set);

  // is_valid_smiles: VARCHAR -> BOOLEAN
  ScalarFunctionSet is_valid_smiles_set("is_valid_smiles");
  is_valid_smiles_set.AddFunction(
      ScalarFunction({LogicalType::VARCHAR}, LogicalType::BOOLEAN, is_valid_smiles));
  loader.RegisterFunction(is_valid_smiles_set);

  // is_valid_smarts: VARCHAR -> BOOLEAN
  ScalarFunctionSet is_valid_smarts_set("is_valid_smarts");
  is_valid_smarts_set.AddFunction(
      ScalarFunction({LogicalType::VARCHAR}, LogicalType::BOOLEAN, is_valid_smarts));
  loader.RegisterFunction(is_valid_smarts_set);
}

} // namespace duckdb_rdkit
