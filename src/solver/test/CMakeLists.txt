package_add_test(burer_solver_gtest burer_solver_gtest.cpp
        ../../solver/burer_solver.cpp
)

package_add_test(enumeration_solver_gtest enumeration_solver_gtest.cpp
        ../../solver/enumeration_solver.cpp
)

package_add_test(unit_weight_solver_gtest unit_weight_solver_gtest.cpp
        ../../solver/unit_weight_solver.cpp
        ../../auxiliary/graphs.cpp
)