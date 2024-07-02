package_add_test(weight_stable_reducer_gtest weight_stable_reducer_gtest.cpp
        ../../prepro/weight_stable_reducer.cpp
)

package_add_test(weighted_reducer_gtest weighted_reducer_gtest.cpp
        ../../prepro/weighted_reducer.cpp
)

package_add_test(data_reduction_gtest data_reduction_gtest.cpp
        ../../auxiliary/graphs.cpp
        ../../prepro/weighted_reducer.cpp
        ../../prepro/weight_stable_reducer.cpp
        ../../prepro/spqr.cpp
        ../../solver/enumeration_solver.cpp
        ../../prepro/similar_subgraph_reduction.cpp
        ../../solver/unit_weight_solver.cpp)

package_add_test(similar_subgraph_reduction_gtest similar_subgraph_reduction_gtest.cpp
        ../../prepro/similar_subgraph_reduction.cpp
        ../../auxiliary/graphs.cpp
)

package_add_test(spqr_gtest spqr_gtest.cpp
        ../../prepro/spqr.cpp
        ../../solver/enumeration_solver.cpp)
