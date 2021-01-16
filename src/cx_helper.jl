
function cx_matrix(velocity::Float64)
    #convert velocity 
    σ = Array{Float64,2}(undef, 6,6)
    ccall((:__cx_helper_MOD_cx_matrix, "./cx_helper.so"), Cvoid, (Ptr{Float64},Ptr{Float64}), Ref(velocity),σ)
    return σ
end