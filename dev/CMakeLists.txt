cmake_minimum_required(VERSION 3.4)

project("curve")

set(ViennaHRLE_DIR $ENV{VIENNAHRLE_DIR})
set(ViennaLS_DIR "${ViennaLS_BINARY_DIR}")
set(ViennaLS_COMMON_TARGET ${PROJECT_NAME})

find_package(ViennaHRLE REQUIRED)
find_package(ViennaLS REQUIRED)

#add_executable(${PROJECT_NAME} ${PROJECT_NAME}.cpp)
#target_include_directories(${PROJECT_NAME} PUBLIC ${VIENNALS_INCLUDE_DIRS})
#target_link_libraries(${PROJECT_NAME} PRIVATE ${VIENNALS_LIBRARIES})

#My test calculations
# add_executable(dMtodE dMtodE.cpp)
# target_include_directories(dMtodE PUBLIC ${VIENNALS_INCLUDE_DIRS})
# target_link_libraries(dMtodE PRIVATE ${VIENNALS_LIBRARIES})

# add_executable(normalizeFMM normalizeFMM.cpp)
# target_include_directories(normalizeFMM PUBLIC ${VIENNALS_INCLUDE_DIRS})
# target_link_libraries(normalizeFMM PRIVATE ${VIENNALS_LIBRARIES})

add_executable(compareCurvatures compareCurvatures.cpp)
target_include_directories(compareCurvatures PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(compareCurvatures PRIVATE ${VIENNALS_LIBRARIES})

add_executable(generalTest generalTest.cpp)
target_include_directories(generalTest PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(generalTest PRIVATE ${VIENNALS_LIBRARIES})

add_executable(flaggingTest flaggingTest.cpp)
target_include_directories(flaggingTest PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(flaggingTest PRIVATE ${VIENNALS_LIBRARIES})

add_executable(iteratorTests iteratorTests.cpp)
target_include_directories(iteratorTests PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(iteratorTests PRIVATE ${VIENNALS_LIBRARIES})

add_executable(StackedNanowire StackedNanowire.cpp)
target_include_directories(StackedNanowire PUBLIC ${VIENNALS_INCLUDE_DIRS})
target_link_libraries(StackedNanowire PRIVATE ${VIENNALS_LIBRARIES})





