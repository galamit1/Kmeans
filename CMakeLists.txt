PROJECT(SPK)
add_subdirectory(goal_implementations)
add_executable (spk spkmeans.c spkmeans.h
        utils.c utils.h
        goal_implementations/wam.c goal_implementations/wam.h)
