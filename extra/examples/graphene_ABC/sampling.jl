
function getcustomsampling(; N::Int=5, k::Int=1, l::Int=20)

    xs = LinRange(0,1,3*N+1)[1:end-1]
    ys = LinRange(0,1,3*N+1)[1:end-1]

    function isin1(p)

        ((xs[N-k] < p[1] <= xs[N+k]) && (ys[2*N-k] < p[2] <= ys[2*N+k]))
    end

    function isin2(p)

        ((xs[2*N-k] < p[1] <= xs[2*N+k]) && (ys[N-k] < p[2] <= ys[N+k]))
    end

    X1 = hcat(([x,y] for x in xs for y in ys if !isin1([x,y])&&!isin2([x,y]))...)
    weights1 = 1/(3*N)^2*ones(size(X1,2))

    xs2 = LinRange(xs[N-k+1], xs[N+k+1], l+1)[1:end-1]
    ys2 = LinRange(ys[2*N-k+1], ys[2*N+k+1], l+1)[1:end-1]
    X2 = hcat(([x,y] for x in xs2 for y in ys2)...)
    weights2 = ((xs[N+k]-xs[N-k])/l)^2*ones(size(X2,2))

    # @assert (xs[N+k]-xs[N-k])≈2*k/(3*N)

    xs3 = LinRange(xs[2*N-k+1], xs[2*N+k+1], l+1)[1:end-1]
    ys3 = LinRange(ys[N-k+1], ys[N+k+1], l+1)[1:end-1]
    X3 = hcat(([x,y] for x in xs3 for y in ys3)...)
    weights3 = ((xs[N+k]-xs[N-k])/l)^2*ones(size(X3,2))

    Xall = [X1 X2 X3]
    weights = [weights1; weights2; weights3]

    Xall, weights
end