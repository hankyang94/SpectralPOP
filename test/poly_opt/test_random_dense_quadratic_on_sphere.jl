using DynamicPolynomials
using SpectralPOP


function test_random_dense_quadratic_on_sphere(n::Int64)


    println("***Problem setting***")


    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n] #variables


    v=reverse(monomials(x,0:2))
    c=2*rand(Float64,length(v)).-1
    f= c'*v #objective function



    R=1.0
    h=[R-sum(x.^2)] #sphere constraints

    l=length(h)

    println("Number of equality constraints: l=",l)
    println("====================")

    k=1 # relaxed order

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


N=[50;75;100;125;150;175;200;250;300;350;400;500;700;900;1200;1500]

for n in N
    test_random_dense_quadratic_on_sphere(n)
end
