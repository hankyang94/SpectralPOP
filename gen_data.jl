using DynamicPolynomials, SparseArrays, LinearAlgebra, JuMP, SumOfSquares, MAT
include("./src/basicfuncs.jl")
include("./src/POP_CTP.jl")

# contype = "sphere";
contype = "binary";
n       = 20;
alpha   = 4;
@polyvar x[1:n]; #variables
v       = reverse(monomials(x,0:alpha));
c       = randn(Float64,length(v));
f       = c'*v;
g       = differentiate(f, x);
R       = 1.0;
if contype == "sphere"
    h   = [R-sum(x.^2)];
elseif contype == "binary"
    h   = x.^2 - ones(n);
end
kappa   = 2;
N,l_h,vv,sk,m,a0,a,Ib,Vb,invInde,invP = ConvertStandardSDP(x,f,h,kappa)

invInde = Float64.(invInde);
filename = string("/Users/hankyang/Documents/MATLAB/sdp-solver-certifiable-perception/RANDPOP/data/",contype,"/testpop.mat");
matwrite(filename, Dict(
	"d"  => n,
	"vv" => vv,
    "n" => sk,
    "m"  => m,
    "C"  => a0,
    "At"  => a,
    "Ib"  => Ib,
    "Vb"  => Vb,
    "invInde" => invInde,
    "invP" => invP,
    "f" => f,
    "g" => g,
    "contype" => contype);)
