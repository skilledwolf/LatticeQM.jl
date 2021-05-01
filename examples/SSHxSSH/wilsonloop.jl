using Plots
using LatticeQM

using LatticeQM.Algebra: σ0, σ1, σ2, σ3
using LatticeQM.Structure.Lattices: Lattice, addbasis!, addorbital!, addorbitals!, addextra!

# Benalcazar, Bernevig, Hughes, Science 357, 61–66 (2017)
function getSSH2D(t1X,t2X, t1Y,t2Y) #t1X,δtX, t1Y,δtY

    lat = Lattice() # 0D lattice
    addbasis!(lat, [1,0]) # 1D lattice
    addbasis!(lat, [0,1]) # 2D lattice
    # addbasis!(lat, [0,0,1], :finite) # 2D lattice with z-coordinates
    addextra!(lat, "sublattice") # non-spatial coordinate
    addorbital!(lat, [0,   0,   2])
    addorbital!(lat, [1/2, 0,   4])
    addorbital!(lat, [0,   1/2, 3])
    addorbital!(lat, [1/2, 1/2, 1])

    lat.specialpoints = LatticeQM.Structure.Geometries.kdict_sq

    h = DenseHops()
    h[[0,0]] = zeros(Complex, 4,4) # intra cell hops
    h[[0,0]][3,1] = t1X
    h[[0,0]][2,3] = -t1Y
    h[[0,0]][4,2] = t1X
    h[[0,0]][1,4] = t1Y
    h[[0,0]] += h[[0,0]]'

    h[[1,0]] = zeros(Complex, 4,4) # inter cell hops in X
    h[[1,0]][3,1] = t2X
    h[[1,0]][2,4] = t2X
    h[[-1,0]] = h[[1,0]]'

    h[[0,1]] = zeros(Complex, 4,4) # inter cell hops in Y
    h[[0,1]][2,3] = -t2Y
    h[[0,1]][4,1] = t2Y
    h[[0,-1]] = h[[0,1]]'

    
    lat, h
end

lat, h = getSSH2D(0.5,1.0,0.5,1.0)

pol1, U1, pol2, U2 = Spectrum.NestedWilson2D(h, 100, 100,1:2)

# en1 = mod.(en1, 1.0)
p = scatter(pol1[:,1])
scatter!(p, pol1[:,2])
mkpath("output"); savefig(p, "output/pol1.pdf")

p = scatter(pol2[:,1])
scatter!(p, pol2[:,2])
mkpath("output"); savefig(p, "output/pol2.pdf")