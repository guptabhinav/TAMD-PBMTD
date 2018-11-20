program projected_free_energy 
IMPLICIT NONE

REAL*8, ALLOCATABLE :: p(:,:,:,:)
REAL*8 :: dum, pp, norm, grid2, grid3
INTEGER*8 :: i,j,k,l, nx, ny, nz, nd

 
open(1,file='input')
open(2, file='PROB_4D',form='unformatted',status='old')

read(1,*) nx, ny, nz, nd
read(1,*) grid2, grid3

ALLOCATE(p(nx,ny,nz,nd))

norm=0.d0
p(nx,ny,nz,nd)=0.d0

do i=1,nx
  do j=1,ny
    do k=1,nz
      do l=1,nd
      read(2) p(i,j,k,l)
      norm = norm + p(i,j,k,l)
      end do
    end do
  end do
end do

print *, 'norm =', norm
norm = norm*grid2*grid3
print *, 'norm*area =', norm

open(10,file='PROB_2D',status='replace')
do j=1,ny
   do k=1,nz
      pp=0.d0
      do i=1,nx
        do l=1,nd
        pp=pp+p(i,j,k,l)
        end do
      end do
      write(10,*) j, k, pp/norm 
    end do
    write(10,*)
end do

end program

