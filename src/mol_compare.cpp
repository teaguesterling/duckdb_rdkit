#include "common.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/execution/expression_executor_state.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "mol_formats.hpp"
#include "types.hpp"
#include "umbra_mol.hpp"
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdio>
#include <memory>

namespace duckdb_rdkit {

// credit: code is from chemicalite
// https://github.com/rvianello/chemicalite
// See mol_search.test for an example of
// a molecule which can return false negative, if the SMILES is different from
// the query
bool mol_cmp(std::string m1_bmol, std::string m2_bmol) {
  std::unique_ptr<RDKit::ROMol> m1(new RDKit::ROMol());
  std::unique_ptr<RDKit::ROMol> m2(new RDKit::ROMol());

  RDKit::MolPickler::molFromPickle(m1_bmol, *m1);
  RDKit::MolPickler::molFromPickle(m2_bmol, *m2);

  // credit: code is from chemicalite
  // https://github.com/rvianello/chemicalite
  // See mol_search.test for an example of
  // a molecule which can return false negative, if the SMILES is different
  // from the query if m1 is substruct of m2 and m2 is substruct of m1,
  // likely to be the same molecule
  RDKit::MatchVectType matchVect;
  bool recursion_possible = false;
  bool do_chiral_match = false; /* FIXME: make configurable getDoChiralSSS(); */
  bool ss1 = RDKit::SubstructMatch(*m1, *m2, matchVect, recursion_possible,
                                   do_chiral_match);
  bool ss2 = RDKit::SubstructMatch(*m2, *m1, matchVect, recursion_possible,
                                   do_chiral_match);
  if (ss1 && !ss2) {
    return false;
  } else if (!ss1 && ss2) {
    return false;
  }

  // the above can still fail in some chirality cases
  std::string smi1 = RDKit::MolToSmiles(*m1, do_chiral_match);
  std::string smi2 = RDKit::MolToSmiles(*m2, do_chiral_match);
  return smi1 == smi2;
}

// is_exact_match for pure Mol type (no prefix optimization)
static void is_exact_match_mol(DataChunk &args, ExpressionState &state,
                               Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_pickle, string_t &right_pickle) {
        // Pure Mol: direct RDKit comparison
        return mol_cmp(left_pickle.GetString(), right_pickle.GetString());
      });
}

// is_exact_match for UmbraMol type (with prefix optimization)
static void is_exact_match_umbramol(DataChunk &args, ExpressionState &state,
                                    Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        auto left = umbra_mol_t(left_umbra_blob);
        auto right = umbra_mol_t(right_umbra_blob);

        // The prefix of a umbra_mol contains a bit vector for substructure
        // screens. We also use this to check exact match. If the molecules
        // being compared do not have the same substructures marked by the
        // dalke_fp, they cannot be an exact match
        if (memcmp(left.GetPrefix(), right.GetPrefix(),
                   umbra_mol_t::PREFIX_BYTES) != 0) {
          return false;
        };

        // otherwise, do the more extensive check with rdkit
        return mol_cmp(left.GetBinaryMol(), right.GetBinaryMol());
      });
}

// Direct RDKit substructure match (used by both Mol and UmbraMol)
bool _is_substruct_rdkit(const std::string &target_pickle, const std::string &query_pickle) {
  std::unique_ptr<RDKit::ROMol> target_mol(new RDKit::ROMol());
  std::unique_ptr<RDKit::ROMol> query_mol(new RDKit::ROMol());

  RDKit::MolPickler::molFromPickle(target_pickle, *target_mol);
  RDKit::MolPickler::molFromPickle(query_pickle, *query_mol);

  RDKit::MatchVectType matchVect;
  bool recursion_possible = true;
  bool do_chiral_match = false; /* FIXME: make configurable getDoChiralSSS(); */
  return RDKit::SubstructMatch(*target_mol, *query_mol, matchVect,
                               recursion_possible, do_chiral_match);
}

// UmbraMol substructure match with DalkeFP optimization
bool _is_substruct_umbramol(umbra_mol_t target, umbra_mol_t query) {
  // Use the extended DalkeFP for early bailout
  auto q_dalke_fp = query.GetDalkeFP();
  auto t_dalke_fp = target.GetDalkeFP();

  // Use dalke_fp_contains for comprehensive screening
  if (!dalke_fp_contains(t_dalke_fp, q_dalke_fp)) {
    return false;
  }

  // query might be substructure of the target -- run a full substructure match
  return _is_substruct_rdkit(target.GetBinaryMol(), query.GetBinaryMol());
}

// is_substruct for pure Mol type (no fingerprint optimization)
static void is_substruct_mol(DataChunk &args, ExpressionState &state,
                             Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_pickle, string_t &right_pickle) {
        return _is_substruct_rdkit(left_pickle.GetString(), right_pickle.GetString());
      });
}

// is_substruct for UmbraMol type (with DalkeFP optimization)
static void is_substruct_umbramol(DataChunk &args, ExpressionState &state,
                                  Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, bool>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        auto left_umbra_mol = umbra_mol_t(left_umbra_blob);
        auto right_umbra_mol = umbra_mol_t(right_umbra_blob);
        return _is_substruct_umbramol(left_umbra_mol, right_umbra_mol);
      });
}

// ============================================================================
// substruct_count - count number of substructure matches
// ============================================================================

// Direct RDKit substructure count
int _substruct_count_rdkit(const std::string &target_pickle, const std::string &query_pickle) {
  std::unique_ptr<RDKit::ROMol> target_mol(new RDKit::ROMol());
  std::unique_ptr<RDKit::ROMol> query_mol(new RDKit::ROMol());

  RDKit::MolPickler::molFromPickle(target_pickle, *target_mol);
  RDKit::MolPickler::molFromPickle(query_pickle, *query_mol);

  std::vector<RDKit::MatchVectType> matches;
  bool uniquify = true;
  bool recursion_possible = true;
  bool do_chiral_match = false;

  RDKit::SubstructMatch(*target_mol, *query_mol, matches, uniquify, recursion_possible, do_chiral_match);
  return static_cast<int>(matches.size());
}

// substruct_count for pure Mol type
static void substruct_count_mol(DataChunk &args, ExpressionState &state,
                                Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, int32_t>(
      left, right, result, args.size(),
      [&](string_t &left_pickle, string_t &right_pickle) {
        return _substruct_count_rdkit(left_pickle.GetString(), right_pickle.GetString());
      });
}

// substruct_count for UmbraMol type (with DalkeFP optimization)
static void substruct_count_umbramol(DataChunk &args, ExpressionState &state,
                                     Vector &result) {
  D_ASSERT(args.ColumnCount() == 2);
  auto &left = args.data[0];
  auto &right = args.data[1];

  BinaryExecutor::Execute<string_t, string_t, int32_t>(
      left, right, result, args.size(),
      [&](string_t &left_umbra_blob, string_t &right_umbra_blob) {
        auto left_umbra_mol = umbra_mol_t(left_umbra_blob);
        auto right_umbra_mol = umbra_mol_t(right_umbra_blob);

        // Use DalkeFP for early bailout
        auto q_dalke_fp = right_umbra_mol.GetDalkeFP();
        auto t_dalke_fp = left_umbra_mol.GetDalkeFP();

        if (!dalke_fp_contains(t_dalke_fp, q_dalke_fp)) {
          return 0;
        }

        return _substruct_count_rdkit(left_umbra_mol.GetBinaryMol(), right_umbra_mol.GetBinaryMol());
      });
}

void RegisterCompareFunctions(ExtensionLoader &loader) {
  // is_exact_match: register for both Mol and UmbraMol
  ScalarFunctionSet set("is_exact_match");
  set.AddFunction(ScalarFunction({duckdb_rdkit::Mol(), duckdb_rdkit::Mol()},
                                 LogicalType::BOOLEAN, is_exact_match_mol));
  set.AddFunction(ScalarFunction({duckdb_rdkit::UmbraMol(), duckdb_rdkit::UmbraMol()},
                                 LogicalType::BOOLEAN, is_exact_match_umbramol));
  loader.RegisterFunction(set);

  // is_substruct: register for both Mol and UmbraMol
  ScalarFunctionSet set_is_substruct("is_substruct");
  set_is_substruct.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol(), duckdb_rdkit::Mol()},
                     LogicalType::BOOLEAN, is_substruct_mol));
  set_is_substruct.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol(), duckdb_rdkit::UmbraMol()},
                     LogicalType::BOOLEAN, is_substruct_umbramol));
  loader.RegisterFunction(set_is_substruct);

  // substruct_count: register for both Mol and UmbraMol
  ScalarFunctionSet set_substruct_count("substruct_count");
  set_substruct_count.AddFunction(
      ScalarFunction({duckdb_rdkit::Mol(), duckdb_rdkit::Mol()},
                     LogicalType::INTEGER, substruct_count_mol));
  set_substruct_count.AddFunction(
      ScalarFunction({duckdb_rdkit::UmbraMol(), duckdb_rdkit::UmbraMol()},
                     LogicalType::INTEGER, substruct_count_umbramol));
  loader.RegisterFunction(set_substruct_count);
}

} // namespace duckdb_rdkit
