import ..TightBinding: Hops

"""
    Creates the hops in triangular lattice (possibly tripartite unitcell).
    Used for effective triangular lattice band models when approximating isolated bands
    in twisted bilayer graphene.
"""

function gethaldanelike(lat; t1, t2, cellrange=3, kwargs...)
    hops0 = Hops(lat, (r1,r2=0)->t_haldanelike(r1,r2;t=t1); cellrange=cellrange)
    addhops!(hops0, lat, (r1,r2=0)->t_haldanelike(r1,r2;d0=√3,t=t2); cellrange=cellrange)

    hops1 = kron(hops0, [1.0 0.0; 0.0 0.0])
    hops2 = kron(Hops(Dict(r=>conj(M) for (r,M) in hops0)), [0.0 0.0; 0.0 1.0]) # time reversal symmetric partner sector

    addhops!(hops1, hops2)

    hops1
end



function t_haldanelike(r1,r2=0.0; t=1.0, d0=1, tol=1e-2)
    θ = -0.1 * π

    δr = r1 .- r2
    δr = [cos(θ) -sin(θ); sin(θ) cos(θ)] * δr[1:2]
    d = norm(δr)

    if d0-tol < d < d0+tol
        if mod(angle(δr[1] + δr[2]*1im), 2π/3)-π/3 >= 0
            return t
        else
            return conj(t)
        end
    else
        return 0
    end
end