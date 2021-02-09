project (mr_tests)

find_package(GTest CONFIG REQUIRED)

add_executable(mr_tests 
	MeshReaderTests.cpp
)

target_link_libraries(mr_tests PRIVATE algo GTest::gtest GTest::gtest_main)
target_compile_definitions(mr_tests PRIVATE "DATA_PATH=\"${CMAKE_CURRENT_SOURCE_DIR}/../data/\"")

include(GoogleTest)
gtest_discover_tests(mr_tests)