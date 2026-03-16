function cx_matrix(velocity::Float64)
    nlevs = Int32(6) #Maximum is 6, but can be lower
    σ = Array{Float64,2}(undef,nlevs,nlevs)
    ccall((:cx_matrix, "./cx_helper.so"), Cvoid, (Ptr{Float64},Ref{Int32},Ptr{Float64}), Ref(velocity),Ref(nlevs),σ)
    return σ
end

function cx_rates(neutral_density::Array{Float64,1},neutral_velocity::Array{Float64,1},ion_velocities::Array{Float64,2})
    nlevs = Int64(6)
    @assert size(ion_velocities,1) == 3
    @assert size(neutral_density,1) == nlevs
    @assert size(neutral_velocity,1) == 3
    num_ions = size(ion_velocities,2)
    rates = Array{Float64,2}(undef,nlevs,num_ions)
    ccall((:cx_rates,"./cx_helper.so"),Cvoid,(Ptr{Float64},Ptr{Float64},Ptr{Float64},Ptr{Int64},Ptr{Int64},Ptr{Float64}),neutral_density,neutral_velocity,ion_velocities,Ref(num_ions),Ref(nlevs),rates)
    return rates
end