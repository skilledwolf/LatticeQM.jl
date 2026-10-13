import ..Paths

function cropcircle(lat::Lattice)
	d = latticedim(lat)
	D = spacedim(lat)
	@assert d==2 "Only implemented for 2d at the moment"

	points = allpositions(lat)
	A = getA(lat)

	lconst = minimum(norm, eachcol(A))

	f = p-> norm(p) < lconst/2

	points = points[:,f.(eachcol(points[1:d,:]))]
	N = size(points, 2)

	return Lattice(Matrix(1.0*I, D, D),
		0,
	    points[1:D,:],
	    points[D+1:end,:],
	    # row-index order — Dict key order would permute labels against rows
	    extralabels=extralabelsinorder(lat),
	    specialpoints=Paths.LabeledPoints(
	        ["γ"],
	        [zeros(0)],
	        ["\$\\gamma\$"],
	        ["γ"]
	    )
	)
end

function circleregion(lat::Lattice, r=30.0)
	d = latticedim(lat)
	D = spacedim(lat)
	@assert d==2 "Only implemented for 2d at the moment"
	
	I = ceil(Int,2r/minimum(norm, eachcol(getA(lat))))

	slat = superlattice(lat, [[I, 0] [0, I]])
	foldPC!(slat)

	cropcircle(slat)
end

# circleregion(lat::Lattice, r=30.0) = fillregion(lat, p-> norm(p) < r)

"""
	fillregion(lat::Lattice, f::Function)

Takes the d-dimensional lattice and tiles the d-dimensional region defined by function f.
Returns a 0-dimensional lattice (single unit cell with no periodicities).

f(p)=true if point p belongs to the region and f(p)=false otherwise

The region must be finite and include the origin. The sampled cells containing
an orbital or their origin in the region must be connected through adjacent
cells of the reduced lattice basis. Disconnected regions should be filled
separately.
"""
function fillregion(lat::Lattice, f::Function)
	d = latticedim(lat)
	D = spacedim(lat)
	@assert d==2 "Only implemented for 2d at the moment"

	unitcell = allpositions(lat)

	A = getA(lat)

	# Include both spanning directions even for unequal lattice constants. A
	# shortest-neighbor shell alone may only span one direction.
	_, U = _lattice_reduce(A)
	steps = [Tuple(U * [i, j]) for i in -1:1 for j in -1:1 if i != 0 || j != 0]
	@assert f(zeros(d)) "The origin must lie in the region."
	queue = [(0, 0)]
	seen = Set(queue)
	points = Vector{Float64}[]
	index = 1
	while index <= length(queue)
		cell = queue[index]
		index += 1
		shift = A * collect(cell)
		occupied = f(shift[1:d])
		for orbital in eachcol(unitcell)
			point = copy(orbital)
			point[1:D] .+= shift
			if f(point[1:d])
				push!(points, point)
				occupied = true
			end
		end
		occupied || continue
		for step in steps
			neighbor = cell .+ step
			if !(neighbor in seen)
				push!(seen, neighbor)
				push!(queue, neighbor)
			end
		end
	end
	points = isempty(points) ? zeros(size(unitcell, 1), 0) : reduce(hcat, points)

	# Return a zero-dimensional lattice object. The unitcell contains the positions of orbitals in the region
	# Original orbital labels (such as sublattice or layer-index) are preserved.
	return Lattice(Matrix(1.0*I, D, D),
		0,
	    points[1:D,:],
	    points[D+1:end,:],
	    # row-index order — Dict key order would permute labels against rows
	    extralabels=extralabelsinorder(lat),
	    specialpoints=Paths.LabeledPoints(
	        ["γ"],
	        [zeros(0)],
	        ["\$\\gamma\$"],
	        ["γ"]
	    )
	)
end

