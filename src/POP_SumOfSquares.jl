function SumofSquares_POP(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
@time begin
    println("**SumOfSquares+Mosek:")
    n=length(x)
    l_g=length(g)
    l_h=length(h)



    ceil_g=[@inbounds ceil(Int64,maxdegree(g[i])/2) for i in 1:l_g]
    sk_g=[@inbounds binomial(k-ceil_g[i]+n,n) for i in 1:l_g]

    ceil_h=[@inbounds ceil(Int64,maxdegree(h[i])/2) for i in 1:l_h]
    s2k_h=[@inbounds binomial(2*(k-ceil_h[i])+n,n) for i in 1:l_h]


    model = SOSModel(optimizer_with_attributes(Mosek.Optimizer, QUIET=true))
    @variable(model, lambda)
    @objective(model, Max, lambda)

    wsos=f-lambda
    sigma_monos = reverse(monomials(x, 0:k))


    @variable(model, sigma0, SOSPoly(sigma_monos))
    wsos-=sigma0
    for i in 1:l_g
        sigma= @variable(model,[1:1], SOSPoly(sigma_monos[1:sk_g[i]]))
            wsos-=sigma[1]*g[i]
    end

    psi_monos = reverse(monomials(x, 0:2*k))
    for i in 1:l_h
        psi=@variable(model, [1:1],Poly(psi_monos[1:s2k_h[i]]))
            wsos-=psi[1]*h[i]
    end


    @constraint(model, wsos==0)
    optimize!(model)
    println(termination_status(model))
    opt_val=objective_value(model)
    println("opt_val=",opt_val)
    end
    return opt_val


end


function SumofSquares_POP_WithExtraction(x::Vector{PolyVar{true}},f::Polynomial{true},g::Vector{Polynomial{true,Float64}},h::Vector{Polynomial{true,Float64}},k::Int64)
@time begin
    println("**SumOfSquares+Mosek:")
    n=length(x)
    l_g=length(g)
    l_h=length(h)



    ceil_g=[@inbounds ceil(Int64,maxdegree(g[i])/2) for i in 1:l_g]
    sk_g=[@inbounds binomial(k-ceil_g[i]+n,n) for i in 1:l_g]

    ceil_h=[@inbounds ceil(Int64,maxdegree(h[i])/2) for i in 1:l_h]
    s2k_h=[@inbounds binomial(2*(k-ceil_h[i])+n,n) for i in 1:l_h]


    model = SOSModel(optimizer_with_attributes(Mosek.Optimizer, QUIET=true))
    @variable(model, lambda)
    @objective(model, Max, lambda)

    wsos=f-lambda
    sigma_monos = reverse(monomials(x, 0:k))

    for i in 1:l_g
        sigma= @variable(model,[1:1], SOSPoly(sigma_monos[1:sk_g[i]]))
            wsos-=sigma[1]*g[i]
    end

    psi_monos = reverse(monomials(x, 0:2*k))
    for i in 1:l_h
        psi=@variable(model, [1:1],Poly(psi_monos[1:s2k_h[i]]))
            wsos-=psi[1]*h[i]
    end



    c = @constraint(model, wsos  in SOSCone(),maxdegree = k)
    optimize!(model)
    println(termination_status(model))
    opt_val=objective_value(model)
    println("opt_val=",opt_val)

    #println("Gram matrix=",value(gram).Q)
    #println("Moment matrix=",moment_matrix(dual(c), gram_monos))
    ν = moment_matrix(c)
    ranktol = 1e-3

    extract=extractatoms(ν, ranktol)
    println(extract)
    if extract!=nothing
        opt_sol =extract.atoms[1].center[:,1]
    else
        opt_sol=Vector{Float64}([])
    end

    end
    return opt_val,opt_sol


end
