function cx_matrix(velocity::Float64)
    nlevs = Int32(6) #Maximum is 6, but can be lower
    σ = Array{Float64,2}(undef,nlevs,nlevs)
    ccall((:cx_matrix, "./cx_helper.so"), Cvoid, (Ptr{Float64},Ref{Int32},Ptr{Float64}), Ref(velocity),Ref(nlevs),σ)
    return σ
end