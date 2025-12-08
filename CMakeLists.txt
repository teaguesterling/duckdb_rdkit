cmake_minimum_required(VERSION 3.18)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# RDKit requires C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Set extension name here
set(TARGET_NAME duckdb_rdkit)

project(${TARGET_NAME})

include(ExternalProject)

# Option to use system RDKit or build from submodule
option(USE_SYSTEM_RDKIT "Use system-installed RDKit instead of submodule" OFF)
option(RDK_BUILD_CPP_TESTS "Build RDKit C++ tests" OFF)

# Find Boost (required by RDKit)
find_package(Boost REQUIRED COMPONENTS system serialization iostreams)

if(USE_SYSTEM_RDKIT)
    find_package(RDKit REQUIRED)
    set(DUCKDB_RDKIT_LIBRARIES
        RDKit::SmilesParse
        RDKit::GraphMol
        RDKit::Descriptors
    )
else()
    # Build RDKit from submodule using ExternalProject
    # This completely isolates RDKit's CMake from DuckDB's export system
    message(STATUS "Building RDKit from submodule using ExternalProject")

    set(RDKIT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/rdkit)
    set(RDKIT_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/rdkit-build)
    set(RDKIT_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/rdkit-install)

    ExternalProject_Add(rdkit_external
        SOURCE_DIR ${RDKIT_SOURCE_DIR}
        BINARY_DIR ${RDKIT_BINARY_DIR}
        INSTALL_DIR ${RDKIT_INSTALL_DIR}
        CMAKE_ARGS
            -DCMAKE_INSTALL_PREFIX=${RDKIT_INSTALL_DIR}
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
            -DCMAKE_CXX_STANDARD=20
            -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
            -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
            -DCMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH}
            -DBoost_DIR=${Boost_DIR}
            -DBoost_INCLUDE_DIR=${Boost_INCLUDE_DIR}
            -DBOOST_ROOT=${Boost_INCLUDE_DIR}/..
            -DRDK_BUILD_PYTHON_WRAPPERS=OFF
            -DRDK_BUILD_SWIG_WRAPPERS=OFF
            -DRDK_BUILD_PGSQL=OFF
            -DRDK_BUILD_CONTRIB=OFF
            -DRDK_BUILD_INCHI_SUPPORT=OFF
            -DRDK_BUILD_AVALON_SUPPORT=OFF
            -DRDK_BUILD_FREESASA_SUPPORT=OFF
            -DRDK_BUILD_YAEHMOP_SUPPORT=OFF
            -DRDK_BUILD_XYZ2MOL_SUPPORT=OFF
            -DRDK_BUILD_STRUCTCHECKER_SUPPORT=OFF
            -DRDK_BUILD_CFFI_LIB=OFF
            -DRDK_BUILD_MINIMAL_LIB=OFF
            -DRDK_BUILD_FUZZ_TARGETS=OFF
            -DRDK_BUILD_CPP_TESTS=${RDK_BUILD_CPP_TESTS}
            -DRDK_TEST_MMFF_COMPLIANCE=OFF
            -DRDK_BUILD_LONG_RUNNING_TESTS=OFF
            -DRDK_INSTALL_INTREE=OFF
            -DRDK_INSTALL_STATIC_LIBS=ON
            -DRDK_BUILD_THREADSAFE_SSS=ON
            -DRDK_BUILD_DESCRIPTORS3D=ON
            -DRDK_BUILD_COORDGEN_SUPPORT=ON
            -DRDK_BUILD_MAEPARSER_SUPPORT=ON
            -DRDK_BUILD_MOLINTERCHANGE_SUPPORT=ON
            -DRDK_USE_BOOST_SERIALIZATION=ON
            -DRDK_USE_BOOST_IOSTREAMS=ON
        BUILD_BYPRODUCTS
            ${RDKIT_INSTALL_DIR}/lib/libRDKitRDGeneral_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitRDGeometryLib_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitRDStreams_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitDataStructs_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitRingDecomposerLib_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitGraphMol_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitSmilesParse_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitFileParsers_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitSubgraphs_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitSubstructMatch_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitChemTransforms_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitFingerprints_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitDescriptors_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitPartialCharges_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitMolTransforms_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitForceField_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitForceFieldHelpers_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitEigenSolvers_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitAlignment_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitMolAlign_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitDistGeometry_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitDistGeomHelpers_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitDepictor_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitcoordgen_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitmaeparser_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitGenericGroups_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitCIPLabeler_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitMolEnumerator_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitCatalogs_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitFilterCatalog_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitFragCatalog_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitMolStandardize_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitTautomerQuery_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitMolInterchange_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitChemReactions_static.a
            ${RDKIT_INSTALL_DIR}/lib/libRDKitO3AAlign_static.a
    )

    # Include RDKit headers - use BEFORE to ensure they come before system includes
    include_directories(BEFORE ${RDKIT_SOURCE_DIR}/Code)
    include_directories(BEFORE ${RDKIT_BINARY_DIR}/Code)
    include_directories(BEFORE ${RDKIT_INSTALL_DIR}/include/rdkit)
    # Also include Boost headers before system to avoid version conflicts
    include_directories(BEFORE ${Boost_INCLUDE_DIR})

    # Create imported targets for RDKit libraries
    set(RDKIT_LIB_DIR ${RDKIT_INSTALL_DIR}/lib)

    # Helper macro to create imported static library
    # RDKit libraries have the "RDKit" prefix in their names (e.g., libRDKitSmilesParse_static.a)
    macro(add_rdkit_imported_lib name)
        add_library(RDKit_${name} STATIC IMPORTED)
        set_target_properties(RDKit_${name} PROPERTIES
            IMPORTED_LOCATION ${RDKIT_LIB_DIR}/libRDKit${name}_static.a)
        add_dependencies(RDKit_${name} rdkit_external)
    endmacro()

    # Core RDKit libraries
    add_rdkit_imported_lib(RDGeneral)
    add_rdkit_imported_lib(RDGeometryLib)
    add_rdkit_imported_lib(RDStreams)
    add_rdkit_imported_lib(DataStructs)
    add_rdkit_imported_lib(RingDecomposerLib)
    add_rdkit_imported_lib(GraphMol)
    add_rdkit_imported_lib(SmilesParse)
    add_rdkit_imported_lib(FileParsers)
    add_rdkit_imported_lib(Subgraphs)
    add_rdkit_imported_lib(SubstructMatch)
    add_rdkit_imported_lib(ChemTransforms)
    add_rdkit_imported_lib(Fingerprints)
    add_rdkit_imported_lib(Descriptors)
    add_rdkit_imported_lib(PartialCharges)
    add_rdkit_imported_lib(MolTransforms)
    add_rdkit_imported_lib(ForceField)
    add_rdkit_imported_lib(ForceFieldHelpers)
    add_rdkit_imported_lib(EigenSolvers)
    add_rdkit_imported_lib(Alignment)
    add_rdkit_imported_lib(MolAlign)
    add_rdkit_imported_lib(DistGeometry)
    add_rdkit_imported_lib(DistGeomHelpers)
    add_rdkit_imported_lib(Depictor)
    add_rdkit_imported_lib(coordgen)
    add_rdkit_imported_lib(maeparser)
    add_rdkit_imported_lib(GenericGroups)
    add_rdkit_imported_lib(CIPLabeler)
    add_rdkit_imported_lib(MolEnumerator)
    add_rdkit_imported_lib(Catalogs)
    add_rdkit_imported_lib(FilterCatalog)
    add_rdkit_imported_lib(FragCatalog)
    add_rdkit_imported_lib(MolStandardize)
    add_rdkit_imported_lib(TautomerQuery)
    add_rdkit_imported_lib(MolInterchange)
    add_rdkit_imported_lib(ChemReactions)
    add_rdkit_imported_lib(O3AAlign)

    # Use linker group to handle circular dependencies between RDKit static libraries
    # This is more robust than manually ordering libraries
    set(RDKIT_STATIC_LIBS
        RDKit_ChemReactions
        RDKit_MolInterchange
        RDKit_TautomerQuery
        RDKit_MolStandardize
        RDKit_FilterCatalog
        RDKit_FragCatalog
        RDKit_Catalogs
        RDKit_MolEnumerator
        RDKit_CIPLabeler
        RDKit_DistGeomHelpers
        RDKit_DistGeometry
        RDKit_O3AAlign
        RDKit_MolAlign
        RDKit_Alignment
        RDKit_EigenSolvers
        RDKit_ForceFieldHelpers
        RDKit_ForceField
        RDKit_Descriptors
        RDKit_PartialCharges
        RDKit_Fingerprints
        RDKit_ChemTransforms
        RDKit_SubstructMatch
        RDKit_Subgraphs
        RDKit_FileParsers
        RDKit_GenericGroups
        RDKit_Depictor
        RDKit_coordgen
        RDKit_maeparser
        RDKit_SmilesParse
        RDKit_MolTransforms
        RDKit_GraphMol
        RDKit_RingDecomposerLib
        RDKit_DataStructs
        RDKit_RDStreams
        RDKit_RDGeometryLib
        RDKit_RDGeneral
    )

    # Wrap RDKit libraries in a linker group to resolve circular dependencies
    if(UNIX AND NOT APPLE)
        set(DUCKDB_RDKIT_LIBRARIES
            -Wl,--start-group
            ${RDKIT_STATIC_LIBS}
            -Wl,--end-group
            Boost::serialization
            Boost::iostreams
            Boost::system
        )
    else()
        # On macOS/Windows, repeat libs to resolve circular dependencies
        set(DUCKDB_RDKIT_LIBRARIES
            ${RDKIT_STATIC_LIBS}
            ${RDKIT_STATIC_LIBS}
            Boost::serialization
            Boost::iostreams
            Boost::system
        )
    endif()
endif()

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

include_directories(src/include)

set(EXTENSION_SOURCES
    src/sdf_scanner/sdf_functions.cpp
    src/sdf_scanner/sdf_scan.cpp
    src/cast.cpp
    src/mol_compare.cpp
    src/mol_formats.cpp
    src/types.cpp
    src/duckdb_rdkit_extension.cpp
    src/umbra_mol.cpp
    src/mol_descriptors.cpp
    src/qed.cpp
    src/rdkit_log.cpp
)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

# Link RDKit libraries
target_link_libraries(${EXTENSION_NAME} ${DUCKDB_RDKIT_LIBRARIES})
target_link_libraries(${LOADABLE_EXTENSION_NAME} ${DUCKDB_RDKIT_LIBRARIES})

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
