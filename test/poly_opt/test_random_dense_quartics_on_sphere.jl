function test_random_dense_quartics_on_sphere()


function test(n::Int64)


    println("***Problem setting***")


    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables
    R=1.0
    
    
    function generate_random_poly(v::Vector{Monomial{true}})
        c=2*rand(Float64,length(v)).-1
        return c'*v
    end

    
    function generate_objective_and_constraints(x)
        # random quadratic objective function f
        v=reverse(monomials(x,0:4))
        f=generate_random_poly(v)


        # unit sphere constraint

        h=[R-sum(x.^2)] #type of coefficients of each polynomial must be float

        l_h=length(h)

        println("Number of equality constraints: l_h=",l_h)
        println("====================")
        return f,h
    end
    
    f,h=generate_objective_and_constraints(x)
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
    test(n)
end
end