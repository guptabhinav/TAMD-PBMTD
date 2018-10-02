! calculating 2D free energy from 4D probability distribution
program free_energy
implicit none
REAL*8, ALLOCATABLE :: prob(:,:), grid_min(:), grid_max(:), grid_width(:)
REAL*8 :: dum1, dum2, norm, beta_cv
REAL*8 :: grid_min1,grid_min2,grid_min3,grid_min4,grid_max1,grid_max2,grid_max3,grid_max4
REAL*8 :: t_cv, cv1, cv2, cv3, cv4, s1, s2
INTEGER :: i, ns ,ns_cv, steps, bin1, bin2, i_s1, i_s2 
INTEGER, ALLOCATABLE ::  nbin(:)
REAL*8, PARAMETER :: kb=8.314472e-3                                                                  ! kJ/(mol*K)
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0

open(1,file='input')

read(1,*) ns
read(1,*) ns_cv
read(1,*) t_cv

open(2,file='prob_4D.dat')

call step_count(2,steps)
print *, 'steps form prob_4D.dat file =', steps

allocate(grid_min(ns_cv))
allocate(grid_max(ns_cv))
allocate(grid_width(ns_cv))
allocate(nbin(ns_cv))

beta_cv = 1.d0/(kb*t_cv)
print *, 'beta cv =', beta_cv

read(1,*) nbin(1:ns_cv)
read(1,*) grid_min(1:ns_cv)
read(1,*) grid_max(1:ns_cv)

! =======================getting minima and maxima from the generated data ===========================

!CALL grid_min_max(2,grid_min1,grid_min2,grid_min3,grid_min4,grid_max1,grid_max2,grid_max3,grid_max4)
!grid_min(1) = grid_min1
!grid_max(1) = grid_max1
!grid_min(2) = grid_min2
!grid_max(2) = grid_max2
!grid_min(3) = grid_min3
!grid_max(3) = grid_max3
!grid_min(4) = grid_min4
!grid_max(4) = grid_max4
!======================================================================================================
print*, 'grid_min1, grid_min2, grid_min3 and grid_min4 =', grid_min(1:ns_cv)
print*, 'grid_max1, grid_max2, grid_max3 and grid_max4 =', grid_max(1:ns_cv)
print *, 'bins used=', nbin(1:ns_cv)


grid_width(1:ns_cv) = (grid_max(1:ns_cv) - grid_min(1:ns_cv))/float(nbin(1:ns_cv)-1)

print *, 'grid widths =', grid_width(1:ns_cv)

allocate(prob(nbin(1),nbin(2)))
prob=0.d0
norm=0.d0

!allocate(cv(ns_cv,steps))

rewind(2)
do i=1,steps
          read(2,*)  cv1, cv2, cv3, cv4, dum2
          bin1 = nint((cv1-grid_min(1))/grid_width(1))+1
          bin2 = nint((cv2-grid_min(2))/grid_width(2))+1
               if(bin1.gt.0.and.bin2.gt.0.and.bin1.le.nbin(1).and.bin2.le.nbin(2)) then
               prob(bin1,bin2) = prob(bin1,bin2) + dum2
               norm = norm + dum2
               end if
end do

print *, 'norm =', norm

do i=1,ns
norm =norm*grid_width(i)
end do

print *, 'norm =', norm
open(10,file='probability_2D.dat',status='unknown')
open(20,file='free_energy.dat',status='unknown')

  DO i_s1=1,nbin(1)
    s1= grid_min(1) + float(i_s1-1)*grid_width(1)
       DO i_s2=1,nbin(2)
       s2= grid_min(2) + float(i_s2-1)*grid_width(2)
       prob(i_s1,i_s2) = prob(i_s1,i_s2)*(1.d0/norm)

       write(10,'(3F16.8)') s1, s2,  prob(i_s1,i_s2) 
       WRITE(20,'(3F16.8)') s1, s2, -(1.d0/beta_cv)*DLOG(MAX(prob(i_s1,i_s2),1e-16))*kj_to_kcal

       END DO
     write(20,*)
   END DO
print *, 'free_energy.dat => S_1, S_2, prob_density, free_energy'

end program

subroutine step_count(file_number,steps)
integer :: file_number, steps, ios,i
steps=0
do
 read(file_number,*,iostat=ios)
 if(ios.ne.0) exit
  steps=steps+1
end do
end subroutine


SUBROUTINE grid_min_max(num,grid_min1,grid_min2,grid_min3,grid_min4,grid_max1,grid_max2,grid_max3,grid_max4)
implicit none
integer :: num, ios
real*8 :: grid_min1, grid_min2, grid_min3, grid_min4
real*8 :: grid_max1, grid_max2, grid_max3, grid_max4
real*8 :: cv1, cv2, cv3, cv4, dumm

rewind(num)
read(num,*,iostat=ios) cv1, cv2, cv3, cv4, dumm
if (ios.ne.0) stop 'error reading colvar file'
grid_min1=cv1
grid_max1=cv1
grid_min2=cv2
grid_max2=cv2
grid_min3=cv3
grid_max3=cv3
grid_min4=cv4
grid_max4=cv4

rloop : DO
     read(num,*,iostat=ios) cv1, cv2, cv3, cv4, dumm
        if(ios.ne.0) exit rloop
     grid_min1=MIN(cv1,grid_min1)
     grid_max1=MAX(cv1,grid_max1)
     grid_min2=MIN(cv2,grid_min2)
     grid_max2=MAX(cv2,grid_max2)
     grid_min3=MIN(cv3,grid_min3)
     grid_max3=MAX(cv3,grid_max3)
     grid_min4=MIN(cv4,grid_min4)
     grid_max4=MAX(cv4,grid_max4)
end do rloop


end subroutine
