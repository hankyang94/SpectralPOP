using DynamicPolynomials
using SpectralPOP


function test_random_dense_QCQP_on_ball(n::Int64)

    println("***Problem setting***")

    l_g=ceil(Int32, n/8)
    l_h=ceil(Int32, n/8)

    println("Number of variable: n=",n)
    println("====================")

    @polyvar x[1:n]# variables



    # random quadratic objective function f
    v=reverse(monomials(x,0:2))
    c=2*rand(Float64,length(v)).-1
    f=c'*v


    # unit sphere constraint
    R=1.0
    g=[R-sum(x.^2)] #type of coefficients of each polynomial must be float



    # random quadratic equality constraints
    randx=2*rand(n).-1# create a feasible solution
    randx=randx./sqrt(sum(randx.^2))
    randx=randx.*rand(1)

    for j in 1:l_g-1
        a=2*rand(Float64,length(v[2:end])).-1
        push!(g,a'*v[2:end])
        g[end]-=g[end](x => randx)#make constraints feasible
        g[end]+=rand(1)[1]
    end
    #g=Polynomial{true,Float64}[]
    l_g=length(g)
    println("Number of inequality constraints: l_g=",l_g)
    println("====================")

    h=Polynomial{true,Float64}[]
    for j in 1:l_h
        a=2*rand(Float64,length(v[2:end])).-1
        push!(h,a'*v[2:end])
        h[end]-=h[end](x => randx) #make constraints feasible
    end
    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")

    k=2

    println("Relaxed order: k=",k)

    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()


    if n<=30

        opt_val = SpectralPOP.SumofSquares_POP(x,f,g,h,k) # SumOfSquares.jl + Mosek

        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()
    end

    opt_val,opt_sol = SpectralPOP.CTP_POP_on_Ball(x,f,g,h,k,R)

    if n<=20
        println()
        println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        println()

        opt_val,opt_sol = SpectralPOP.BTP_POP(x,f,g,h,k,R,EigAlg="Mix",tol=1e-5,showNormGrad=false)
    end
    println()
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
    println()

end


N=[10;15;20;25;30;35]

for n in N
    test_random_dense_QCQP_on_ball(n)
end
