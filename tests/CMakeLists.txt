set(SRC_FILES test_binomial_expectation test_binomial_expectation.cpp ${CMAKE_SOURCE_DIR}/src/binomial_expectation.cpp)
add_executable(test_binomial_expectation ${SRC_FILES})
target_include_directories(test_binomial_expectation PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/extern/catch2)
add_test(NAME test_binomial_expectation COMMAND test_binomial_expectation)