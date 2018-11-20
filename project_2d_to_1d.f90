program projected_free_energy 
IMPLICIT NONE

REAL*8, ALLOCATABLE :: p(:,:)
REAL*8 :: dum, pp1, pp2,  norm, beta_cv
REAL*8 :: s1, s2, grid_width
INTEGER*8 :: i,j, nx, ny
REAL*8, PARAMETER :: kb=8.314472e-3
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0

 
open(1,file='input')
open(2, file='PROB_2D', STATUS='OLD')

read(1,*) nx, ny, beta_cv
read(1,*) grid_width

ALLOCATE(p(nx,ny))

grid_width = (3.14-(-3.14))/float(nx-1)

print *, 'grid width =', grid_width

norm=0.d0
do i=1,nx
  do j=1,ny
       read(2,*) dum,dum,p(i,j)
       norm=norm+p(i,j)
  end do
 read(2,*)
end do

print *, 'norm =', norm
norm = norm*grid_width
print *, 'norm =', norm

open(10,file='PROB_1D',status='replace')

  do i=1,nx
  s1= -3.14 + float(i-1)*grid_width
  pp1=0.d0
      do j=1,ny
      pp1=pp1+p(i,j)
      end do
  pp1 = pp1/norm
  write(10,*) s1, pp1
  end do


print *, 'done'
close(10)
end program

