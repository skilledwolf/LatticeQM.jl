function minsum(x::AbstractVector,y::AbstractVector)
	if length(y) > length(x)
		minsum(y,x)
	else
		N = minimum([length(x),length(y)])
		return [x[1:N] .+ y[1:N]; x[N+1:end]]
	end
end

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
	    extralabels=collect(keys(lat.extralabels)), 
	    specialpoints=LabeledPoints(
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

Note: f must define a finite region that includes the origin.
"""
function fillregion(lat::Lattice, f::Function)
	d = latticedim(lat)
	D = spacedim(lat)
	@assert d==2 "Only implemented for 2d at the moment"

	unitcell = allpositions(lat)

	A = getA(lat)

	δA = hcat(getneighborcells(lat; halfspace=false, innerpoints=true, excludeorigin=true)...)


	function getpoints()

		@assert f(zeros(d)) "The origin must lie in the region."
		candidates = Vector{Int64}[zeros(d)]
		confirmed = Vector{Int64}[]
		dismissed = Vector{Int64}[]

		while length(candidates) > 0
			for p in candidates

				if f(A * p) # the candidate is confirmed
					found = true
					append!(confirmed, [p]) # ... added to the list

					# Generate new candidates
					for δx in eachcol(δA) # ... and it's (unknown) neighbors proposed for the next step
						pnew = p + δx

						if !in(pnew,confirmed) && !in(pnew,dismissed) && !in(pnew,candidates)
							append!(candidates, [pnew])
						end
					end
				else
					append!(dismissed, [p])
				end

				deleteat!(candidates, findall(x->x==p, candidates)) # ... remove from candidates
			end
		end

		confirmed
	end

	# points = A*hcat(getpoints()...)

	# mysum(M,p) = (M[1:D,:].+p[1:D]; M)

	# points = hcat([mysum(points,p0) for p0 = eachcol(unitcell)]...)
	# points = points[:,f.(eachcol(points[1:d,:]))]
	# N = size(points, 2)

	points0 = getpoints()
	points = Vector{Float64}[]
	for p in points0
		for p0 = eachcol(unitcell)
			p1 = minsum(p0, A *p)

			if f(p1[1:d])
				append!(points, [p1])
			end
		end
	end
	N = length(points)
	points = hcat(points...)

	# points = vcat(points, hcat(fill(extracoordinates(lat), N)...))
	# extralabels = [["x$i" for i=1:d]; collect(keys(lat.extralabels))[sortperm(collect(values(lat.extralabels)))]]

	# Return a zero-dimensional lattice object. The unitcell contains the positions of orbitals in the region
	# Original orbital labels (such as sublattice or layer-index) are preserved.
	return Lattice(Matrix(1.0*I, D, D),
		0,
	    points[1:D,:],
	    points[D+1:end,:], 
	    extralabels=collect(keys(lat.extralabels)), 
	    specialpoints=LabeledPoints(
	        ["γ"],
	        [zeros(0)],
	        ["\$\\gamma\$"],
	        ["γ"]
	    )
	)
end



# function fillregion(lat::Lattice, f::Function)
# """
# Takes the d-dimensional lattice and tiles the d-dimensional region defined function f.
# Returns a 0-dimensional lattice (single unit cell with no periodicities).

# f(p)=true if point p belongs to the region and f(p)=false otherwise

# Note: f must define a finite region that includes the origin.
# """

# 	d = latticedim(lat)
# 	@assert d==2 "Only implemented for 2d at the moment"

# 	points = Vector{Float64}[] # vector of vectors

# 	unitcell = allpositions(lat)
# 	A = getA(lat)





# 	function addpoints!(points, args...)
# 		found = false

# 		if length(args)==d # inner-most recursion level
# 			x0 = A * [args...]
# 			for x in eachcol(unitcell)
# 				p = minsum(x, x0)
# 				if f(p[1:d])
# 					append!(points, [p])
# 					found = true
# 				end 
# 			end

# 			return found

# 		elseif length(args)<d # still recurse deeper!

# 			for i = Iterators.countfrom(0)
# 				if !addpoints!(points, args..., i)
# 					break
# 				else
# 					found = true
# 				end
# 			end

# 			for i = Iterators.countfrom(-1,-1)
# 				if !addpoints!(points, args..., i)
# 					break
# 				else 
# 					found = true
# 				end
# 			end

# 			# return found

# 		else
# 			error("Woah, something went seriously wrong here.")
# 		end

# 		found
# 	end

# 	addpoints!(points)

# 	N = length(points)

# 	if N == 0
# 		error("The region defined your region functions does not contain any points")
# 	end

# 	points = hcat(points...)
# 	# points = vcat(points, hcat(fill(extracoordinates(lat), N)...))
# 	extralabels = [["x$i" for i=1:d]; collect(keys(lat.extralabels))[sortperm(collect(values(lat.extralabels)))]]

# 	# Return a zero-dimensional lattice object. The unitcell contains the positions of orbitals in the region
# 	# Original orbital labels (such as sublattice or layer-index) are preserved.
# 	return Lattice(zeros(0,0), 
# 	    zeros(0,N), 
# 	    points, 
# 	    extralabels=extralabels, 
# 	    specialpoints=LabeledPoints(
# 	        ["γ"],
# 	        [zeros(0)],
# 	        ["\$\\gamma\$"],
# 	        ["γ"]
# 	    )
# 	)
# end