program grid_min_max
implicit none
integer :: num, ios, steps
real*8 :: grid_min1, grid_min2, grid_max1, grid_max2
real*8 :: m, cv1, cv2, dumm

!open(1,file='prob_4D.dat')
open(1,file='COLVAR_PB')
steps=0
read(1,*,iostat=ios) m,cv1,cv2,dumm,dumm,dumm
if (ios.ne.0) stop 'error reading colvar file'
grid_min1=cv1
grid_max1=cv1
grid_min2=cv2
grid_max2=cv2

rloop : DO
     read(1,*,iostat=ios)m,cv1,cv2,dumm,dumm,dumm
        if(ios.ne.0) exit rloop
        steps=steps+1
     grid_min1=MIN(cv1,grid_min1)
     grid_max1=MAX(cv1,grid_max1)
     grid_min2=MIN(cv2,grid_min2)
     grid_max2=MAX(cv2,grid_max2)
end do rloop
print *, 'steps =', steps
print *,' grid_min1 =', grid_min1
print *,' grid_min2 =', grid_min2
print *,' grid_max1 =', grid_max1
print *,' grid_max2 =', grid_max2

end program
