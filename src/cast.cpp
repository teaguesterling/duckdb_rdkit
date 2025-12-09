#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/exception/conversion_exception.hpp"
#include "duckdb/common/operator/cast_operators.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/cast/default_casts.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmartsWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <duckdb/parser/parsed_data/create_scalar_function_info.hpp>

namespace duckdb_rdkit {

// VARCHAR -> Mol: Parse SMILES and create pure RDKit pickle
// This is consistent with the RDKit Postgres cartridge behavior
bool VarcharToMolCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  bool all_converted = true;
  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      source, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          // Create pure RDKit pickle (no DalkeFP prefix)
          auto pickle = rdkit_mol_to_binary_mol(*mol);
          return StringVector::AddStringOrBlob(result, pickle);
        } catch (...) {
          std::string error_msg = StringUtil::Format(
              "Could not convert string '%s' to Mol", smiles.GetString());
          if (parameters.strict) {
            throw ConversionException(error_msg);
          }
          HandleCastError::AssignError(error_msg, parameters);
          all_converted = false;
          mask.SetInvalid(idx);
          return string_t();
        }
      });
  return all_converted;
}

// VARCHAR -> UmbraMol: Parse SMILES and create UmbraMol format [8B DalkeFP][RDKit Pickle]
bool VarcharToUmbraMolCast(Vector &source, Vector &result, idx_t count,
                           CastParameters &parameters) {
  bool all_converted = true;
  UnaryExecutor::ExecuteWithNulls<string_t, string_t>(
      source, result, count,
      [&](string_t smiles, ValidityMask &mask, idx_t idx) {
        try {
          auto mol = rdkit_mol_from_smiles(smiles.GetString());
          auto umbra_mol = get_umbra_mol_string(*mol);
          return StringVector::AddStringOrBlob(result, umbra_mol);
        } catch (...) {
          std::string error_msg = StringUtil::Format(
              "Could not convert string '%s' to UmbraMol", smiles.GetString());
          if (parameters.strict) {
            throw ConversionException(error_msg);
          }
          HandleCastError::AssignError(error_msg, parameters);
          all_converted = false;
          mask.SetInvalid(idx);
          return string_t();
        }
      });
  return all_converted;
}

// Mol -> VARCHAR: Pure RDKit pickle to SMILES
bool MolToVarcharCast(Vector &source, Vector &result, idx_t count,
                      CastParameters &parameters) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t pickle) {
        // Mol is a pure RDKit pickle - convert directly to SMILES
        auto rdkit_mol = rdkit_binary_mol_to_mol(pickle.GetString());
        auto smiles = rdkit_mol_to_smiles(*rdkit_mol);
        return StringVector::AddString(result, smiles);
      });
  return true;
}

// UmbraMol -> VARCHAR: Extract pickle from UmbraMol format and convert to SMILES
bool UmbraMolToVarcharCast(Vector &source, Vector &result, idx_t count,
                           CastParameters &parameters) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t b_umbra_mol) {
        // UmbraMol format: [8B DalkeFP][RDKit Pickle]
        auto umbra_mol = umbra_mol_t(b_umbra_mol);
        auto bmol = umbra_mol.GetBinaryMol();
        auto rdkit_mol = rdkit_binary_mol_to_mol(bmol);
        auto smiles = rdkit_mol_to_smiles(*rdkit_mol);
        return StringVector::AddString(result, smiles);
      });
  return true;
}

// Mol -> BLOB: Pass through raw binary data (Mol is already a BLOB alias)
bool MolToBlobCast(Vector &source, Vector &result, idx_t count,
                   CastParameters &parameters) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t pickle) {
        return StringVector::AddStringOrBlob(result, pickle.GetString());
      });
  return true;
}

// UmbraMol -> BLOB: Pass through raw binary data (UmbraMol is already a BLOB alias)
bool UmbraMolToBlobCast(Vector &source, Vector &result, idx_t count,
                        CastParameters &parameters) {
  UnaryExecutor::Execute<string_t, string_t>(
      source, result, count, [&](string_t umbra_mol) {
        return StringVector::AddStringOrBlob(result, umbra_mol.GetString());
      });
  return true;
}

void RegisterCasts(ExtensionLoader &loader) {
  // Mol casts (pure RDKit pickle)
  loader.RegisterCastFunction(LogicalType::VARCHAR, duckdb_rdkit::Mol(),
                              BoundCastInfo(VarcharToMolCast), 1);
  loader.RegisterCastFunction(duckdb_rdkit::Mol(), LogicalType::VARCHAR,
                              BoundCastInfo(MolToVarcharCast), 1);
  loader.RegisterCastFunction(duckdb_rdkit::Mol(), LogicalType::BLOB,
                              BoundCastInfo(MolToBlobCast), 1);

  // UmbraMol casts ([8B DalkeFP][RDKit Pickle])
  loader.RegisterCastFunction(LogicalType::VARCHAR, duckdb_rdkit::UmbraMol(),
                              BoundCastInfo(VarcharToUmbraMolCast), 1);
  loader.RegisterCastFunction(duckdb_rdkit::UmbraMol(), LogicalType::VARCHAR,
                              BoundCastInfo(UmbraMolToVarcharCast), 1);
  loader.RegisterCastFunction(duckdb_rdkit::UmbraMol(), LogicalType::BLOB,
                              BoundCastInfo(UmbraMolToBlobCast), 1);
}

} // namespace duckdb_rdkit
