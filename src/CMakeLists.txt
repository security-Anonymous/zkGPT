aux_source_directory(. gpt_src)
list(FILTER gpt_src EXCLUDE REGEX "main*")
option(ENABLE_MY_ASSERT "Turn on MY_ASSERT checks" ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
add_library(gpt_lib ${gpt_src})


add_executable(demo_llm_run main_demo_llm.cpp)
target_link_libraries(demo_llm_run gpt_lib mcl mclbn256)

# add_executable(test_range_prover_run main_test_range_proof.cpp)
# target_link_libraries(test_range_prover_run gpt_lib mcl mclbn256)