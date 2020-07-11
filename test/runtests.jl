test_examples()=include("./test_examples.jl")


#Polynomial optimization
test_random_dense_quadratic_on_sphere()=include("./poly_opt/test_random_dense_quadratic_on_sphere.jl")
test_random_dense_equality_constrained_QCQP_on_sphere_first_order()=include("./poly_opt/test_random_dense_equality_constrained_QCQP_on_sphere_first_order.jl")
test_random_dense_equality_constrained_QCQP_on_sphere_second_order()=include("./poly_opt/test_random_dense_equality_constrained_QCQP_on_sphere_second_order.jl")
test_random_dense_QCQP_unique_inequality_ball_constraint()=include("./poly_opt/test_random_dense_QCQP_unique_inequality_(ball)_constraint.jl")
test_random_dense_QCQP_on_ball()=include("./poly_opt/test_random_dense_QCQP_on_ball.jl")
test_random_dense_quartics_on_sphere()=include("./poly_opt/test_random_dense_quartics_on_sphere.jl")
Evaluation_comparisons()=include("./poly_opt/Evaluation_comparisons.jl")
Norm_Subgrad()=include("./poly_opt/Norm_Subgrad.jl")

#Polynomial systems
test_chemkin()=include("./poly_sys/chemkin.jl")
test_d1()=include("./poly_sys/d1.jl")
test_des22_24()=include("./poly_sys/des22_24.jl")
test_i1()=include("./poly_sys/i1.jl")
test_katsura10()=include("./poly_sys/katsura10.jl")
test_katsura9()=include("./poly_sys/katsura9.jl")
test_katsura8()=include("./poly_sys/katsura8.jl")
test_katsura7()=include("./poly_sys/katsura7.jl")
test_katsura6()=include("./poly_sys/katsura6.jl")
test_kin1()=include("./poly_sys/kin1.jl")
test_ku10()=include("./poly_sys/ku10.jl")
test_pole27sys()=include("./poly_sys/pole27sys.jl")
test_pole28sys()=include("./poly_sys/pole28sys.jl")
test_stewgou()=include("./poly_sys/stewgou.jl")
