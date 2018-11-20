! PBMTD along phi1 and phi2 and TAMD along phi1,phi2,psi1,psi2
PROGRAM reweighting
implicit none
real*8, allocatable :: hh(:), ss(:), grid_max(:), grid_min(:), grid_width(:), dum(:), ds2(:), dt(:), diff_s2(:)
real*8, allocatable :: v(:,:), hill(:,:), ht(:,:), width(:,:),  grid(:,:), ct(:), pb_bias(:), num(:), den(:)
real*8, allocatable :: prob(:,:,:,:), cv(:,:),vbias(:)
real*8 :: kt,ktb,bias_fact,dum4,dum5,dum3,norm,t_cv,t_sys,deltaT,alpha_sys,alpha_cv,beta_sys,beta_cv,addition,dummy                 
real*8 :: gridmin1,gridmin2, griddif1, griddif2, rweight, s1, s2, s3, s4, dum1,dum2, sys_temp                                                       
real*8 :: dumm1, dumm2
integer, allocatable :: nbin(:)
integer :: i_mtd,mtd_steps,ns,ios,is,ig,i,j,i_md,md_steps,mtd_max,dum10
integer :: w_cv, w_hill, k, t_min, t_max, index1, index2, index3, index4, i_s1, i_s2,i_s3, i_s4
REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)
REAL*8, PARAMETER :: kb=8.314472e-3                                                                  ! kJ/(mol*K)
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0

open(1,file='input')

read(1,*) ns
read(1,*) t_cv, t_sys, deltaT                                   ! CVs_temp, system_temperature, bias_factor
read(1,*) w_cv, w_hill                                ! frequency of printing cvmdck file, frequency of hill added 
read(1,*) t_min, t_max


open(2,file='HILLS_phi1')
open(3,file='HILLS_phi3')
open(9,file='HILLS_phi4')
open(4,file='COLVAR')
open(7,file='vbias.dat')
open(8,file='ct.dat')

open(50,file='PROB_4D',form='unformatted' )


CALL step_count(2,mtd_steps)
print *, 'mtd steps from count =', mtd_steps
rewind(2)

CALL step_count(4,md_steps)
print *, 'md steps from count =', md_steps
rewind(4)


allocate(width(ns,mtd_steps))
allocate(ht(ns,mtd_steps))
allocate(hill(ns,mtd_steps))
allocate(hh(ns))
allocate(ss(ns))
allocate(grid_max(ns))
allocate(grid_min(ns))
allocate(grid_width(ns))
allocate(nbin(ns))
allocate(ct(mtd_steps))
allocate(ds2(ns))
allocate(diff_s2(ns))
allocate(dt(ns))
allocate(num(ns))
allocate(den(ns))


read(1,*) grid_max(1:ns)
read(1,*) grid_min(1:ns)
read(1,*) nbin(1:ns)
read(1,*) mtd_steps

print *, 'bins used=', nbin(1:ns)

grid_width(1:ns) = (grid_max(1:ns) - grid_min(1:ns))/float(nbin(1:ns)-1)

print *, 'grid widths =', grid_width(1:ns)


do i=1,mtd_steps
  read(2,*)   dum2, hill(1,i), width(1,i), ht(1,i), dum4
  IF( hill(1,i) .gt. pi) hill(1,i) = hill(1,i) - 2.d0*pi
  IF( hill(1,i) .lt. -pi) hill(1,i) = hill(1,i) + 2.d0*pi
end do

do i=1,mtd_steps
  read(3,*)   dum3, hill(2,i), width(2,i), ht(2,i), dum5 
  IF( hill(2,i) .gt. pi) hill(2,i) = hill(2,i) - 2.d0*pi
  IF( hill(2,i) .lt. -pi) hill(2,i) = hill(2,i) + 2.d0*pi
end do


do i=1,mtd_steps
  read(9,*)   dum3, hill(3,i), width(3,i), ht(3,i), dum5
  IF( hill(3,i) .gt. pi) hill(3,i) = hill(3,i) - 2.d0*pi
  IF( hill(3,i) .lt. -pi) hill(3,i) = hill(3,i) + 2.d0*pi
end do


print *, 'mtd steps used =', mtd_steps
!md_steps = mtd_steps*w_hill/w_cv +1
!md_steps= mtd_steps*w_hill/w_cv

print *, 'md steps used =', md_steps

allocate(pb_bias(md_steps))
allocate(cv(ns,md_steps))

do k = 1, md_steps
     read(4,*) dum3, cv(1,k), dum3, cv(2,k) ,dum3, cv(3,k), dum3, cv(4,k), dum3
              if (cv(1,k) .gt. pi)  cv(1,k) = cv(1,k) - 2.d0*pi
              if (cv(1,k) .lt. -pi) cv(1,k) = cv(1,k) + 2.d0*pi
              if (cv(2,k) .gt. pi)  cv(2,k) = cv(2,k) - 2.d0*pi
              if (cv(2,k) .lt. -pi) cv(2,k) = cv(2,k) + 2.d0*pi
              if (cv(3,k) .gt. pi)  cv(3,k) = cv(3,k) - 2.d0*pi
              if (cv(3,k) .lt. -pi) cv(3,k) = cv(3,k) + 2.d0*pi
              if (cv(4,k) .gt. pi)  cv(4,k) = cv(4,k) - 2.d0*pi
              if (cv(4,k) .lt. -pi) cv(4,k) = cv(4,k) + 2.d0*pi
end do



alpha_cv = (t_cv+deltaT)/deltaT
alpha_sys=(t_sys+deltaT)/deltaT
beta_cv =  1.d0/(kb*t_cv)
beta_sys = 1.d0/(kb*t_sys)

print *, 'alpha_cv = (t_cv+deltaT)/deltaT = ', alpha_cv
print *, 'alpha_sys = (t_sys+deltaT)/deltaT = ', alpha_sys
print *, 'beta_cv = 1/(kb*t_cv) =  ', beta_cv
print *, 'beta_sys = 1/(kb*t_sys) =  ', beta_sys

allocate (grid(ns,nbin(1)))

do i=1,nbin(1)
 grid(1,i) = grid_min(1) + float(i-1)*grid_width(1)
end do

do i=1,nbin(2)
 grid(2,i) = grid_min(2) + float(i-1)*grid_width(2)
end do

do j=1,nbin(3)
 grid(3,j) = grid_min(3) + float(j-1)*grid_width(3)
end do

do j=1,nbin(4)
 grid(4,j) = grid_min(4) + float(j-1)*grid_width(4)
end do

allocate(v(ns,nbin(1)))


allocate(vbias(md_steps))

print *, 'reading ct factor.'
do i=1,mtd_steps
 read(8,'(3I10,F16.6)') dum10, dum10, dum10, ct(i)
end do

print *, 'reading vbias.'
do k=1,md_steps
 read(7,'(I10,F16.6)') dum10, vbias(k)
end do



print *, 'calculating probability.'
allocate(prob(nbin(1),nbin(2),nbin(3),nbin(4)))
norm=0.0
prob=0.0
print *, 'calculating unbiased probabilty density => free energy.dat'

DO i_md=1,md_steps
        IF((i_md.GE.t_min).AND.(i_md.LT.t_max))THEN                                      ! t_min = 1 not 0
          index1 = nint((cv(1,i_md)-grid_min(1))/grid_width(1))+1            ! pbmtd applied
          index2 = nint((cv(2,i_md)-grid_min(2))/grid_width(2))+1            ! US is applied
          index3 = nint((cv(3,i_md)-grid_min(3))/grid_width(3))+1            ! pbmtd applied
          index4 = nint((cv(4,i_md)-grid_min(4))/grid_width(4))+1            ! pbmtd applied
        
           if(index1.gt.0.and.index2.gt.0.and.index3.gt.0.and.index1.le.nbin(1).and.index2.le.&
              nbin(2).and.index3.le.nbin(3).and.index4.gt.0.and.index4.le.nbin(4))then
              i_mtd =  (i_md-1)*w_cv/w_hill 

!             i_mtd= i_md*w_cv/w_hill
!             print *, 'i_md, i_mtd =', i_md, i_mtd 

              dummy = vbias(i_mtd) + ct(i_mtd)             
              rweight = dexp(beta_cv*dummy)
              prob(index1,index2,index3,index4) = prob(index1,index2,index3,index4) + rweight
              norm=norm+rweight
           end if
        
        END IF
END DO



 print *,' norm =', norm
 do is=1,ns
 norm=norm*grid_width(is)
 end do                         
 print *, 'norm per unit area=', norm

    DO i_s1=1,nbin(1)
       DO i_s2=1,nbin(2)
          DO i_s3=1,nbin(3)
              DO i_s4=1,nbin(4)

               prob(i_s1,i_s2,i_s3,i_s4) = prob(i_s1,i_s2,i_s3,i_s4)*(1.d0/norm)                                    
              
               WRITE(50) prob(i_s1,i_s2,i_s3,i_s4)

              END DO
          END DO
       END DO
    END DO

print *, 'free_energy.dat => S_1, S_2,S_3, S_4, prob_density'
close (50)

print *, 'done.'
close (1)
close (2)
close (3)


deallocate(cv)
deallocate(width)
deallocate(ht)
deallocate(hill)
deallocate(hh)
deallocate(ss)
deallocate(nbin)
deallocate(grid_max)
deallocate(grid_min)
deallocate(grid_width)
deallocate(pb_bias)
deallocate(ct)
deallocate(ds2)
deallocate(diff_s2)
deallocate(grid)
deallocate(num)
deallocate(den)
deallocate(prob)

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
