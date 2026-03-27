cmake_minimum_required(VERSION 4.0.3)
project (HCFD LANGUAGES CXX)
set(CMAKE_SHARED_LIBRARY_PREFIX "")
set(CMAKE_STATIC_LIBRARY_PREFIX "")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(includes)

find_package(OpenMP)

add_compile_options(
    -march=native 
    -O3 
    -fopenmp 
    -fno-math-errno
    # -flto 
    # -fopt-info-inline-missed
    -Winline
    # -fopt-info-vec-all
)
#    -fopt-info-vec-optimized -finline-limit=1000 -fopt-info-vec-all -fno-math-errno

set (     CORE      "src/Test.cpp" "src/Mesh.cpp" "src/Writer.cpp")

file(GLOB BASE_BCS  "src/Base/BoundaryConditions/*.cpp")
file(GLOB BASE_COMM "src/Base/*.cpp")

file(GLOB SOA_BCS   "src/SoA/BoundaryConditions/*.cpp")
file(GLOB SOA_COMM  "src/SoA/*.cpp")

file(GLOB VEC_BCS   "src/VEC/BoundaryConditions/*.cpp")
file(GLOB VEC_COMM  "src/VEC/*.cpp")

file(GLOB CAOS_BCS   "src/CAoS/BoundaryConditions/*.cpp")
file(GLOB CAOS_COMM  "src/CAoS/*.cpp")

file(GLOB FAOSCAOS_BCS   "src/FAoS_CAoS/BoundaryConditions/*.cpp")
file(GLOB FAOSCAOS_COMM  "src/FAoS_CAoS/*.cpp")

file(GLOB AOSOA_BCS   "src/AoSoA/BoundaryConditions/*.cpp")
file(GLOB AOSOA_COMM  "src/AoSoA/*.cpp")

add_executable(
    Test ${CORE}  
    ${BASE_BCS} ${BASE_COMM}  
    ${SOA_BCS} ${SOA_COMM}  
    ${VEC_COMM} ${VEC_BCS}  
    ${CAOS_COMM} ${CAOS_BCS}  
    ${FAOSCAOS_COMM} ${FAOSCAOS_BCS}  
    ${AOSOA_COMM} ${AOSOA_BCS}
)
target_link_libraries(Test m)

