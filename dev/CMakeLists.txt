cmake_minimum_required(VERSION 3.4)

project("curve")

set(ViennaHRLE_DIR $ENV{VIENNAHRLE_DIR})
set(ViennaLS_DIR "${ViennaLS_BINARY_DIR}")
set(ViennaLS_COMMON_TARGET ${PROJECT_NAME})

find_package(ViennaHRLE REQUIRED)
find_package(ViennaLS REQUIRED)

add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNALS_LIBRARIES})