function test_random_dense_quadratic_on_sphere()
using DynamicPolynomials
using SpectralPOP


function test(n::Int64)


    println("***Problem setting***")


    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n] #variables

    
    function generate_random_poly(v::Vector{Monomial{true}})
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end

    R=1.0
    
    function generate_objective_and_constraint(x)
        v=reverse(monomials(x,0:2))
        f=generate_random_poly(v) #objective function




        h=[R-sum(x.^2)] #sphere constraints

        l=length(h)

        println("Number of equality constraints: l=",l)
        println("====================")
        return f,h
    end
    
    f,h=generate_objective_and_constraint(x)
    
    k=1 # relaxed order

    println("Relaxed order: k=",k)


    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=300
        g=Vector{Polynomial{true,Float64}}([])

        opt_val = SpectralPOP.SumofSquares_POP(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5) #Limited memory bundle method

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

    opt_val,opt_sol = SpectralPOP.CTP_POP(x,f,h,k,R;method="SketchyCGAL",EigAlg="Normal",tol=1e-4)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


N=[50;75;100;125;150;175;200;250;300;350;400;500;700;900;1200;1500]

for n in N
    test(n)
end
end