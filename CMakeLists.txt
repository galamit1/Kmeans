PROJECT(SPK)

add_subdirectory(goal_implementations)
add_executable (spk spkmeans.c spkmeans.h
        utils.c utils.h
        goal_implementations/wam.c goal_implementations/wam.h
        goal_implementations/ddg.c goal_implementations/ddg.h
        goal_implementations/lnorm.c goal_implementations/lnorm.h
        )
