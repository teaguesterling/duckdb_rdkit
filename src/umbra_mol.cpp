#include "umbra_mol.hpp"
#include "common.hpp"
#include "duckdb/common/exception.hpp"
#include "mol_formats.hpp"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <cstdint>
#include <string>

namespace duckdb_rdkit {

// These are fragments and the number of times they appear in a molecule.
// Use a vector to keep the order. The order of the vector will inform the
// bit order.
// "O" 2 times, is bit 0,
// "O" 3 times is bit 1, etc...
// This was calculated by Andrew Dalke:
// http://www.dalkescientific.com/writings/diary/archive/2012/06/11/optimizing_substructure_keys.html
// And Greg Landrum tested out the 55 bits I use here.
// https://www.mail-archive.com/rdkit-discuss@lists.sourceforge.net/msg02078.html
//
// This is used for a substructure filter and placed in the prefix of a
// Umbra-mol so that short-circuiting can take place and save on computation
// cost when deserializing and a full substructure search is unnecessary.
std::vector<std::vector<string>> dalke_counts = {{"O", "2", "3", "1", "4", "5"},
                                                 {"Ccc", "2", "4"},
                                                 {"CCN", "1"},
                                                 {"cnc", "1"},
                                                 {"cN", "1"},
                                                 {"C=O", "1"},
                                                 {"CCC", "1"},
                                                 {"S", "1"},
                                                 {"c1ccccc1", "1", "2"},
                                                 {"N", "2", "3", "1"},
                                                 {"C=C", "1"},
                                                 {"nn", "1"},
                                                 {"CO", "2"},
                                                 {"Ccn", "1", "2"},
                                                 {"CCCCC", "1"},
                                                 {"cc(c)c", "1"},
                                                 {"CNC", "2"},
                                                 {"s", "1"},
                                                 {"CC(C)C", "1"},
                                                 {"o", "1"},
                                                 {"cncnc", "1"},
                                                 {"C=N", "1"},
                                                 {"CC=O", "2", "3"},
                                                 {"Cl", "1"},
                                                 {"ccncc", "2"},
                                                 {"CCCCCC", "6"},
                                                 {"F", "1"},
                                                 {"CCOC", "3"},
                                                 {"c(cn)n", "1"},
                                                 {"C", "9", "6", "1"},
                                                 {"CC=C(C)C", "1"},
                                                 {"c1ccncc1", "1"},
                                                 {"CC(C)N", "1"},
                                                 {"CC", "1"},
                                                 {"CCC(C)O", "4"},
                                                 {"ccc(cc)n", "2"},
                                                 {"C1CCCC1", "1"},
                                                 {"CNCN", "1"},
                                                 {"cncn", "3"},
                                                 {"CSC", "1"},
                                                 {"CCNCCCN", "1"},
                                                 {"CccC", "1"},
                                                 {"ccccc(c)c", "3"}};

// Convert heavy atom count to a 4-bit bucket (0-15)
// Ranges designed for typical drug-like molecules
static uint8_t heavy_atom_bucket(unsigned int count) {
  if (count <= 5) return 0;
  if (count <= 10) return 1;
  if (count <= 15) return 2;
  if (count <= 20) return 3;
  if (count <= 25) return 4;
  if (count <= 30) return 5;
  if (count <= 35) return 6;
  if (count <= 40) return 7;
  if (count <= 50) return 8;
  if (count <= 60) return 9;
  if (count <= 75) return 10;
  if (count <= 90) return 11;
  if (count <= 110) return 12;
  if (count <= 140) return 13;
  if (count <= 180) return 14;
  return 15;  // 181+
}

// Convert ring count to 2-bit value (0, 1, 2, 3+)
static uint8_t ring_count_bucket(unsigned int count) {
  if (count >= 3) return 3;
  return static_cast<uint8_t>(count);
}

uint64_t make_dalke_fp(const RDKit::ROMol &mol) {
  uint64_t fp = 0;

  // === Bits 0-54: Original Dalke fragment patterns ===
  RDKit::SubstructMatchParameters params;
  params.uniquify = true;
  params.useQueryQueryMatches = false;
  params.recursionPossible = true;
  params.useChirality = false;
  params.maxMatches = 10;
  params.numThreads = 1;

  uint8_t curBit = 0;
  for (const auto &frag : dalke_counts) {
    std::unique_ptr<RDKit::ROMol> dalke_fp_mol;

    try {
      dalke_fp_mol.reset(RDKit::SmilesToMol(frag[0], 0, false));
    } catch (std::exception &e) {
      std::string msg = StringUtil::Format("%s", typeid(e).name());
      throw InvalidInputException(msg);
    }

    auto matchVect = RDKit::SubstructMatch(mol, *dalke_fp_mol, params);

    for (auto i = 1; i < frag.size(); i++) {
      if (matchVect.size() >= std::stoi(frag[i])) {
        fp |= (1ULL << curBit);
      }
      curBit++;
    }
  }
  D_ASSERT(curBit == 55);

  // === Bits 55-58: Heavy atom count bucket (4 bits) ===
  unsigned int numHeavyAtoms = mol.getNumHeavyAtoms();
  uint8_t ha_bucket = heavy_atom_bucket(numHeavyAtoms);
  fp |= (static_cast<uint64_t>(ha_bucket) << 55);

  // === Bits 59-60: Ring count bucket (2 bits) ===
  unsigned int numRings = mol.getRingInfo()->numRings();
  uint8_t ring_bucket = ring_count_bucket(numRings);
  fp |= (static_cast<uint64_t>(ring_bucket) << 59);

  // === Bit 61: Has stereocenters ===
  bool hasStereo = false;
  for (const auto atom : mol.atoms()) {
    if (atom->getChiralTag() != RDKit::Atom::CHI_UNSPECIFIED) {
      hasStereo = true;
      break;
    }
  }
  if (hasStereo) {
    fp |= (1ULL << 61);
  }

  // === Bit 62: Has formal charges ===
  bool hasCharges = false;
  for (const auto atom : mol.atoms()) {
    if (atom->getFormalCharge() != 0) {
      hasCharges = true;
      break;
    }
  }
  if (hasCharges) {
    fp |= (1ULL << 62);
  }

  // === Bit 63: Reserved (leave as 0) ===

  return fp;
}

// "Umbra-mol" has more than just the binary molecule
// There is a prefix in front of the binary molecule, inspired by
// Umbra-style strings
std::string get_umbra_mol_string(const RDKit::ROMol &mol) {
  auto binary_mol = rdkit_mol_to_binary_mol(mol);
  size_t total_size = umbra_mol_t::DALKE_FP_PREFIX_BYTES + binary_mol.size();

  uint64_t dalke_fp = make_dalke_fp(mol);

  // remember to keep endianess in mind if you print things out.
  // little endian on my machine
  std::string buffer;
  buffer.reserve(total_size);
  buffer.append(reinterpret_cast<const char *>(&dalke_fp),
                umbra_mol_t::DALKE_FP_PREFIX_BYTES);
  buffer.append(binary_mol);

  return buffer;
}
} // namespace duckdb_rdkit
