using DynamicPolynomials
using SpectralPOP


function test_random_dense_quartics_on_sphere(n::Int64)


    println("***Problem setting***")


    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables

    # random quadratic objective function f
    v=reverse(monomials(x,0:4))
    c=2*rand(Float64,length(v)).-1
    f=c'*v


    # unit sphere constraint
    R=1.0
    h=[R-sum(x.^2)] #type of coefficients of each polynomial must be float

    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")

    k=Int64(2)

    println("Relaxed order: k=",k)


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()



    g=Vector{Polynomial{true,Float64}}([])

    opt_val = SpectralPOP.SumofSquares_POP(x,f,g,h,k) # SumOfSquares.jl + Mosek

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5) #Limited memory bundle method

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-3)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


N=[5;10;15;20;25;30]

for n in N
    test_random_dense_quartics_on_sphere(n)
end
