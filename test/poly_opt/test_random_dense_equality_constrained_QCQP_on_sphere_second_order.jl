using DynamicPolynomials
using SpectralPOP


function test_random_dense_equality_constrained_QCQP_on_sphere_second_order(n::Int64)


    println("***Problem setting***")

    l=ceil(Int32, n/4)

    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    # random quadratic objective function f
    v=reverse(monomials(x,0:2))
    c=2*rand(Float64,length(v)).-1
    f=c'*v

    R=1.0
    # unit sphere constraint
    h=[R-sum(x.^2)] #type of coefficients of each polynomial must be float

    # random quadratic equality constraints
    randx=2*rand(n).-1# create a feasible solution
    randx=randx./sqrt(sum(randx.^2))


    for j in 1:l
        a=2*rand(Float64,length(v[2:end])).-1
        push!(h,a'*v[2:end])
        h[end]-=h[end](x => randx) #make constraints feasible
    end
    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")

    k=Int64(2)

    println("Relaxed order: k=",k)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=30
        g=Vector{Polynomial{true,Float64}}([])

        opt_val = SpectralPOP.SumofSquares_POP(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5) #Limited memory bundle method

    if n<=30
        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()

        opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-3)
    end
    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


N=[5;10;15;20;25;30;35;40;45;50;55;60;65]

for n in N
    test_random_dense_quadratic_on_sphere_second_order(n)
end
