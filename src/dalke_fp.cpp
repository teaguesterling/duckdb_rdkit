#include "common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function_set.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"

namespace duckdb_rdkit {

// dalke_fp(Mol) -> DalkeFP
// Compute the 64-bit Dalke fingerprint from a Mol (pure RDKit pickle)
void dalke_fp_from_mol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &mol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, uint64_t>(
      mol_vec, result, count,
      [&](string_t mol_blob, ValidityMask &mask, idx_t idx) {
        try {
          // Parse the RDKit pickle
          auto mol = rdkit_binary_mol_to_mol(mol_blob.GetString());
          if (!mol) {
            mask.SetInvalid(idx);
            return uint64_t(0);
          }
          return make_dalke_fp(*mol);
        } catch (...) {
          mask.SetInvalid(idx);
          return uint64_t(0);
        }
      });
}

// dalke_fp(UmbraMol) -> DalkeFP
// Extract the embedded DalkeFP from an UmbraMol (stored in first 8 bytes)
void dalke_fp_from_umbramol(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 1);
  auto &umbramol_vec = args.data[0];
  auto count = args.size();

  UnaryExecutor::ExecuteWithNulls<string_t, uint64_t>(
      umbramol_vec, result, count,
      [&](string_t umbramol_blob, ValidityMask &mask, idx_t idx) {
        if (umbramol_blob.GetSize() < 8) {
          mask.SetInvalid(idx);
          return uint64_t(0);
        }
        // Extract the first 8 bytes as uint64_t
        uint64_t fp = 0;
        std::memcpy(&fp, umbramol_blob.GetData(), sizeof(uint64_t));
        return fp;
      });
}

// dalke_fp_contains(target DalkeFP, query DalkeFP) -> BOOLEAN
// Returns true if target fingerprint might contain query as substructure
void dalke_fp_contains_func(DataChunk &args, ExpressionState &state, Vector &result) {
  D_ASSERT(args.data.size() == 2);
  auto &target_vec = args.data[0];
  auto &query_vec = args.data[1];
  auto count = args.size();

  BinaryExecutor::Execute<uint64_t, uint64_t, bool>(
      target_vec, query_vec, result, count,
      [&](uint64_t target_fp, uint64_t query_fp) {
        return dalke_fp_contains(target_fp, query_fp);
      });
}

void RegisterDalkeFPFunctions(ExtensionLoader &loader) {
  // dalke_fp(Mol) -> DalkeFP
  ScalarFunctionSet dalke_fp_set("dalke_fp");
  dalke_fp_set.AddFunction(
      ScalarFunction({Mol()}, DalkeFP(), dalke_fp_from_mol));
  dalke_fp_set.AddFunction(
      ScalarFunction({UmbraMol()}, DalkeFP(), dalke_fp_from_umbramol));
  loader.RegisterFunction(dalke_fp_set);

  // dalke_fp_contains(target, query) -> BOOLEAN
  ScalarFunctionSet dalke_fp_contains_set("dalke_fp_contains");
  dalke_fp_contains_set.AddFunction(
      ScalarFunction({DalkeFP(), DalkeFP()}, LogicalType::BOOLEAN, dalke_fp_contains_func));
  loader.RegisterFunction(dalke_fp_contains_set);
}

} // namespace duckdb_rdkit
