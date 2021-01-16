neutral_density = [1.0,0.,0.,0.,0.,0.]
neutral_velocity = [0.,0.,0.]
start = 6
stop = 9
length = 50
function logspace(start,stop,length=50)
    exp10.(range(start, stop=stop, length=length))
end
ion_velocities = Array{Float64,2}(undef,3,length)
ion_velocities[2:end,:] .= 0
ion_velocities[1,:] = logspace(start,stop,length)
cx_rate_results = cx_rates(neutral_density,neutral_velocity,ion_velocities)

