println("***POP with single inequality (ball) constraint***")

n=10

println("Number of variable: n=",n)
println("====================")

@polyvar x[1:n]# variables

function generate_random_poly(v)
    c=2*rand(Float64,length(v)).-1
    return c'*v
end
# random quadratic objective function f

R=1.0
function generate_objective_and_constraints(x)
    l_h=ceil(Int32, n/4)
    
    v=reverse(monomials(x,0:2))


    f=generate_random_poly(v)


    # random quadratic equality constraints
    randx=2*rand(n).-1# create a feasible solution
    randx=randx./sqrt(sum(randx.^2))
    randx=randx.*rand(1)

    
    g=[R-sum(x.^2)]

    println("Number of inequality constraints: l_g=",1)
    println("====================")

    h=Polynomial{true,Float64}[]
    for j in 1:l_h
        push!(h,generate_random_poly(v[2:end]))
        h[end]-=h[end](x => randx) #make constraints feasible
    end
    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")
    return f,g,h
end

f,g,h=generate_objective_and_constraints(x)

k=2

println("Relaxed order: k=",k)



opt_val,opt_sol = SpectralPOP.CTP_POP_UniqueBallIneq(x,f,h,k,R;method="LMBM",EigAlg="Arpack",tol=1e-5,showNormGrad=true)


println()
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println()



println("***POP on a ball***")



println("Number of variable: n=",n)
println("====================")

@polyvar x[1:n]# variables

R=1.0
function generate_objective_and_constraints2(x)
    l_g=ceil(Int32, n/8)
    l_h=ceil(Int32, n/8)
    
    # random quadratic objective function f
    v=reverse(monomials(x,0:2))

    f=generate_random_poly(v)


    # unit sphere constraint

    g=[R-sum(x.^2)] #type of coefficients of each polynomial must be float



    # random quadratic equality constraints
    randx=2*rand(n).-1# create a feasible solution
    randx=randx./sqrt(sum(randx.^2))
    randx=randx.*rand(1)

    for j in 1:l_g-1
        push!(g,generate_random_poly(v[2:end]))
        g[end]-=g[end](x => randx)#make constraints feasible
        g[end]+=rand(1)[1]
    end
    #g=Polynomial{true,Float64}[]
    l_g=length(g)
    println("Number of inequality constraints: l_g=",l_g)
    println("====================")

    h=Polynomial{true,Float64}[]
    for j in 1:l_h
        push!(h,generate_random_poly(v[2:end]))
        h[end]-=h[end](x => randx) #make constraints feasible
    end
    l_h=length(h)

    println("Number of equality constraints: l_h=",l_h)
    println("====================")
    return f,g,h
end

f,g,h=generate_objective_and_constraints2(x)

k=2

println("Relaxed order: k=",k)

println()
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println()


opt_val,opt_sol = SpectralPOP.CTP_POP_on_Ball(x,f,g,h,k,R,showNormGrad=true)
