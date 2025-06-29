function RKStep(f,t0,x0,h,A,b)

    # A Butcher tabelaux
    # b weights
    # h step size
    # x0 init value (column)
    # t0 init time
    # f function handle of RHS accepting f(float,column vector)

    n = length(x0)
    s=size(A,1)
    c=A*ones(s)


    # each stage is a column vector
    # all those are placed next to each other hKold, hKnew

    tol = 1e-14

    hKold=h*zeros(n,s);
    fcast(tt,xx) = hcat(f.(eachcol(tt),eachcol(xx))...)
    hKnew=h*fcast(t0.+h*c', x0 .+ hKold*A')

    while norm(hKnew-hKold,2) > tol
        hKold = hKnew
        hKnew = h*fcast(t0.+h*c', x0 .+ hKold*A')
    end

    return x0 .+ hKnew*b

end


function RKStep(f,t0,x0,h,A,b,N::Int)

    # A Butcher tabelaux
    # b weights
    # h step size
    # x0 init value (column)
    # t0 init time
    # f function handle of RHS accepting f(float,column vector)
    # N number of steps

    Step(x) = RKStep(f,t0,x,h,A,b)

    # pre-allocation
    x = zeros(length(x0),N+1)
    x[:,1] = x0
    
    # iteration
    @showprogress for k=1:N
           x[:,k+1] = Step(x[:,k])
    end

    return x
end