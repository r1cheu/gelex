# --------------------------------------------------
# ------------ Core Library ------------------------
# --------------------------------------------------

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS "*.cpp")
list(REMOVE_ITEM CORE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")

add_library(gelex_core ${CORE_SOURCES})
add_library(gelex::core ALIAS gelex_core)

set_target_properties(gelex_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(
  gelex_core
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> $<INSTALL_INTERFACE:include>
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_include_directories(
  gelex_core SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/ext
                           ${PROJECT_SOURCE_DIR}/ext/eigen)

target_compile_definitions(
  gelex_core PUBLIC BARKEEP_ENABLE_FMT_FORMAT
                    $<$<CONFIG:Release>:EIGEN_NO_DEBUG>)

if(GELEX_USE_MKL)
  target_compile_definitions(gelex_core PUBLIC EIGEN_USE_MKL_ALL USE_MKL)
else()
  target_compile_definitions(gelex_core PUBLIC EIGEN_USE_BLAS EIGEN_USE_LAPACK)
endif()

target_compile_options(
  gelex_core PRIVATE -Wall -Wextra -Wpedantic
                     $<$<BOOL:${GELEX_NATIVE_OPTIMIZATION}>:-march=native -O3>)

target_link_libraries(
  gelex_core PUBLIC OpenMP::OpenMP_CXX BLAS::BLAS spdlog::spdlog
                    $<$<BOOL:${MINGW}>:ws2_32>)

if(GELEX_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(gelex_core PUBLIC -fsanitize=undefined
                                           -fsanitize=address)
  target_link_options(gelex_core PUBLIC -fsanitize=undefined -fsanitize=address)
endif()

if(GELEX_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(gelex_core PUBLIC -fprofile-instr-generate
                                           -fcoverage-mapping)
  target_link_options(gelex_core PUBLIC -fprofile-instr-generate
                      -fcoverage-mapping)
endif()

# --------------------------------------------------
# ------------ Main Executable ---------------------
# --------------------------------------------------
add_executable(gelex_cli main.cpp)
set_target_properties(gelex_cli PROPERTIES OUTPUT_NAME "gelex")

target_link_libraries(gelex_cli PRIVATE gelex::core)
