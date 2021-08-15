PROJECT(SPK)
add_subdirectory(goal_implementations)
add_executable (spk spkmeans.c spkmeans.h goal_implementations/spk.h goal_implementations/spk.c)
