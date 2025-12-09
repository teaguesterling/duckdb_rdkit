#pragma once
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include <cstddef>
#include <cstdint>
#include <sys/types.h>

namespace duckdb_rdkit {

// Mol: Pure RDKit MolPickler binary (NEW default)
LogicalType Mol();

// UmbraMol: [8B DalkeFP][RDKit Pickle] - single-column optimized format
LogicalType UmbraMol();

// DalkeFP: 64-bit substructure screening fingerprint
LogicalType DalkeFP();

// MolStruct: STRUCT(mol BLOB, dalke_fp UBIGINT)
// Stores molecule (RDKit pickle) and fingerprint separately for flexible columnar access
LogicalType MolStruct();

void RegisterTypes(ExtensionLoader &loader);

} // namespace duckdb_rdkit
