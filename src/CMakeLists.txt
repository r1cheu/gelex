# --------------------------------------------------
# ------------ Core Library ------------------------
# --------------------------------------------------

file(GLOB_RECURSE CORE_SOURCES CONFIGURE_DEPENDS "*.cpp")

add_library(gelex ${CORE_SOURCES})
add_library(gelex::gelex ALIAS gelex)

set_target_properties(gelex PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(
  gelex
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> $<INSTALL_INTERFACE:include>
  PRIVATE "${PROJECT_SOURCE_DIR}/src")
target_include_directories(gelex SYSTEM PRIVATE ${PROJECT_SOURCE_DIR}/ext/mio)

target_compile_definitions(gelex PUBLIC $<$<CONFIG:Release>:EIGEN_NO_DEBUG>)

if(GELEX_USE_MKL)
  target_compile_definitions(gelex PUBLIC EIGEN_USE_MKL_ALL USE_MKL)
else()
  target_compile_definitions(gelex PUBLIC EIGEN_USE_BLAS EIGEN_USE_LAPACKE)
endif()

target_compile_options(
  gelex PRIVATE -Wall -Wextra -Wpedantic
                $<$<BOOL:${GELEX_NATIVE_OPTIMIZATION}>:-march=native -O3>)

target_link_libraries(
  gelex PUBLIC OpenMP::OpenMP_CXX BLAS::BLAS LAPACKE::LAPACKE fmt::fmt
               Eigen3::Eigen $<$<BOOL:${MINGW}>:ws2_32>)

if(GELEX_ENABLE_SANITIZERS AND CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  target_compile_options(gelex PUBLIC -fsanitize=undefined -fsanitize=address)
  target_link_options(gelex PUBLIC -fsanitize=undefined -fsanitize=address)
endif()

if(GELEX_ENABLE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  target_compile_options(gelex PUBLIC -fprofile-instr-generate
                                      -fcoverage-mapping)
  target_link_options(gelex PUBLIC -fprofile-instr-generate -fcoverage-mapping)
endif()
