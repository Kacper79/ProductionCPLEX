# CMakeList.txt : CMake project for ProductionCPLEX, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("ProductionCPLEX")
set(CMAKE_CXX_STANDARD 20)

add_executable (ProductionCPLEX "src/ProductionColGenCPLEX.cpp") #   "src/ProductionBaseCPLEX.cpp"

set(CPLEX_ROOT_PATH "E:/CPLEX")

target_include_directories(ProductionCPLEX PUBLIC "${CMAKE_SOURCE_DIR}/include")
target_include_directories(ProductionCPLEX PUBLIC "${CPLEX_ROOT_PATH}/concert/include")
target_include_directories(ProductionCPLEX PUBLIC "${CPLEX_ROOT_PATH}/cplex/include")

target_link_directories(ProductionCPLEX PUBLIC "${CPLEX_ROOT_PATH}/concert/lib/x64_windows_msvc14/stat_mdd")
target_link_directories(ProductionCPLEX PUBLIC "${CPLEX_ROOT_PATH}/cplex/lib/x64_windows_msvc14/stat_mdd")

target_link_libraries(ProductionCPLEX PUBLIC concert cplex2212 ilocplex)


