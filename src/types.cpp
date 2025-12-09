// This file defines the new types in a way for duckdb to recognize

#include "types.hpp"
#include "common.hpp"
#include "duckdb/common/constants.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb_rdkit {

// Mol: Pure RDKit MolPickler binary (NEW default)
// This is a standard RDKit pickle format, directly interoperable
LogicalType Mol() {
  auto blob_type = LogicalType(LogicalTypeId::BLOB);
  blob_type.SetAlias("Mol");
  return blob_type;
}

// UmbraMol: [8B DalkeFP prefix][RDKit Pickle]
// Single-column format with embedded fingerprint for optimized substructure search
LogicalType UmbraMol() {
  auto blob_type = LogicalType(LogicalTypeId::BLOB);
  blob_type.SetAlias("UmbraMol");
  return blob_type;
}

// DalkeFP: 64-bit substructure screening fingerprint
// Bits 0-54: Dalke fragment patterns
// Bits 55-58: Heavy atom count bucket
// Bits 59-60: Ring count (0, 1, 2, 3+)
// Bit 61: Has stereocenters
// Bit 62: Has charges
// Bit 63: Reserved
LogicalType DalkeFP() {
  auto uint64_type = LogicalType(LogicalTypeId::UBIGINT);
  uint64_type.SetAlias("DalkeFP");
  return uint64_type;
}

void RegisterTypes(ExtensionLoader &loader) {
  // Register all molecule-related types
  loader.RegisterType("Mol", Mol());
  loader.RegisterType("UmbraMol", UmbraMol());
  loader.RegisterType("DalkeFP", DalkeFP());
}

} // namespace duckdb_rdkit
