program fidatest

use fidanet

call setinputs(ai,ab,6)
call settables()

call setstates([5.6d8,4.5d7,6.2d6,4.3d5,2.1d5,1.6d4],states)

call testcol(1,2.3d0,1.2d-4,states,dens,photons)

print*, states
print*, dens
print*, photons

end program fidatest
