# UNIFAC-excess
UNIFAC described in terms of excess gibbs free energy

a lot of thermodynamic libraries implement UNIFAC, and in particular, the residual part, using the following definition:

```julia
Ge = R*T*sum(x[i]*γ[i] for i in components)
γ[i] = γ(comb)*γ(res)
#we focus on the residual part. the combinatorial part is a GC-averaged UNIQUAC
γ(res)[i] = sum(ν[i][k]*(Γ[k] - Γ[i][k]) for k in groups)
ln(Γ[k]) = Qk(1 - log(sum(Θ[k]*Ψ[m,k] for k in groups) - sum(Θ[m]*Ψ[k,m]/sum(Θ[n]*Ψ[n,m] for n in groups) for k in groups)
Ψ[i,j] = f(T) #depends on the formulation.
Θ[m] = Q[m]X[m]/sum(Q[m]*X[m] for m in groups)
X[m] = sum(x[i]*ν[i,m] for i in components)/X̄
X̄ = sum(sum(x[i]*ν[i,k] for i in components) for k in groups)
```
The evaluation of the excess gibbs (or helmholtz, because of pressure independence) free energy is needed, especially on EoS/Gᴱ mixing rules, but in the UNIFAC, the formulation presented above has more complexity than necessary. requiring evaluation (and allocation) of vectors containing activity coefficients. in particular, to avoid calculating Ψ[k,m] n times,it's necessary to allocate the entire `Ψ` matrix of group interaction parameters (O(N^2) memory). Compare this with the expression of excess gibbs free energy for UNIQUAC:

```julia
Ge(uniquac,T,z) = G_comb(uniquac,T,z) +  G_res(uniquac,T,z)
G_comb(uniquac,T,z) = sum(x[i]*log(ϕ[i]/x[i])) + (z/2)*q[i]*x[i]*log(θ[i]/ϕ[i]) 
G_res(uniquac,T,z) =  -sum(x[i]*q[i]*log(sum(θ[j]*τ[j,i]))
z = 5
ϕ[i] = r[i]*x[i]/sum(r[i]*x[i] for i in components)
θ[i] = q[i]*x[i]/sum(q[i]*x[i] for i in components)
```
That is just a quadratic accumulation loop, and if written carefully, does not allocate any temporary vectors.

UNIFAC, can indeed be rewritten in terms of excess gibbs energy.
```
Ge(UNIFAC,T,z) = Ge_comb(UNIFAC,T,z) + Ge_res(UNIFAC,T,z)
```
the combinatiorial part is exactly the same as UNIQUAC. but with GC-averaged `r` and `q`
```
r[i] = sum(ν[i,k]*R[k] for k in groups)
q[i] = sum(ν[i,k]*Q[k] for k in groups)
```

`r` and `q` can be preallocated, they dont change with any combination of (T,x).

the residual part can be rewritten as:

```
Ge_res(UNIFAC,T,z) = -X̄*sum(X[k]*Q[k]*log(sum(Θ[m]*Ψ[m,k] for m in groups) for k in groups)
```
in this particular form, `Ψ[m,k]` is just used once per combination of `(m,k)` indices, that allows us to skip the allocation of the `Ψ` and calculate each index on the fly. the complexity of the operation depends on the decision of allocating or not the vector of segment fractions, `X`. if X is allocated (`O(k)` memory), then the operation `Ge_res(UNIFAC,T,z)` depends quadratically on the number of groups.
