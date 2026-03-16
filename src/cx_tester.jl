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

using Plots
plotly()

plot(ion_velocities[1,:], cx_rate_results[1,:],ticks = :native)
#Need to transpose to get alignment of shapes for automatically plotting multiple lines, one for each energy level
plot(ion_velocities[1,:], transpose(cx_rate_results),ticks = :native)
cxmin = minimum(cx_rate_results[cx_rate_results.>0.0])
plot(ion_velocities[1,:], transpose(cx_rate_results),xaxis=:log,yaxis=:log,ylim=[cxmin,maximum(cx_rate_results)])

