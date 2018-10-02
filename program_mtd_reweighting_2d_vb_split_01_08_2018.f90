! PBMTD along phi1 and phi2 and TAMD along phi1,phi2,psi1,psi2
PROGRAM reweighting
implicit none
real*8, allocatable :: hh(:), ss(:), grid_max(:), grid_min(:), grid_width(:), dum(:), ds2(:), dt(:), diff_s2(:)
real*8, allocatable :: v(:,:), hill(:,:), ht(:,:), width(:,:),  grid(:,:), ct(:), pb_bias(:), num(:), den(:)
real*8, allocatable :: prob(:,:), cv(:,:),vbias(:)
real*8 :: kt,ktb,bias_fact,dum4,dum5,dum3,norm,t_cv,t_sys,deltaT,alpha_sys,alpha_cv,beta_sys,beta_cv,addition,dummy                 
real*8 :: gridmin1,gridmin2, griddif1, griddif2, rweight, s1, s2, dum1,dum2, sys_temp                                                       
real*8 :: dumm1, dumm2
integer, allocatable :: nbin(:)
integer :: i_mtd,mtd_steps,ns,ios,is,ig,i,j,i_md,md_steps,mtd_max,dum10
integer :: w_cv, w_hill, k, t_min, t_max, index1, index2, i_s1, i_s2
REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)
REAL*8, PARAMETER :: kb=8.314472e-3                                                                  ! kJ/(mol*K)
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0

open(1,file='input')

read(1,*) ns
read(1,*) t_cv, t_sys, deltaT                                   ! CVs_temp, system_temperature, bias_factor
read(1,*) w_cv, w_hill                                ! frequency of printing cvmdck file, frequency of hill added 
read(1,*) t_min, t_max


open(2,file='HILLS_phi1')
open(3,file='HILLS_phi2')
open(4,file='COLVAR_PB')

open(20,file='ct.dat',status='unknown')
open(30,file='vbias.dat',status='unknown' )
open(50,file='free_energy.dat',status='unknown' )


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
!read(1,*) mtd_steps

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


print *, 'mtd steps used =', mtd_steps
md_steps = mtd_steps*w_hill/w_cv +1
print *, 'md steps used =', md_steps

allocate(pb_bias(md_steps))
allocate(cv(ns,md_steps))

do k = 1, md_steps 
  read(4,*)   dum3, cv(1,k), cv(2,k), dum3, dum3, pb_bias(k)
              if (cv(1,k) .gt. pi)  cv(1,k) = cv(1,k) - 2.d0*pi 
              if (cv(1,k) .lt. -pi) cv(1,k) = cv(1,k) + 2.d0*pi 
              if (cv(2,k) .gt. pi)  cv(2,k) = cv(2,k) - 2.d0*pi
              if (cv(2,k) .lt. -pi) cv(2,k) = cv(2,k) + 2.d0*pi
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

do j=1,nbin(2)
grid(2,j) = grid_min(2) + float(j-1)*grid_width(2)
end do

allocate(v(ns,nbin(1)))

print *, 'calculating ct factor.'

      v=0.d0
      DO i_mtd=1,mtd_steps
         ds2(1:ns)=width(1:ns,i_mtd)*width(1:ns,i_mtd)
         ss(1:ns)=hill(1:ns,i_mtd)
         hh(1:ns)=ht(1:ns,i_mtd)/alpha_sys                  
         num=0.D0
         den=0.D0
         DO is=1,ns

            DO ig=1,nbin(is)
            diff_s2(is)=grid(is,ig)-ss(is)

               if (diff_s2(is) .gt. pi ) diff_s2(is) =diff_s2(is) - 2.d0*pi  
               if (diff_s2(is) .lt.-pi ) diff_s2(is) =diff_s2(is) + 2.d0*pi 

            diff_s2(is)=diff_s2(is)*diff_s2(is)*0.5D0
            v(is,ig) = v(is,ig) + hh(is)*dexp(-diff_s2(is)/ds2(is))               
            num(is) = num(is) + dexp(beta_cv*(alpha_sys-1.d0)*v(is,ig))                 
            den(is) = den(is) + dexp(beta_cv*alpha_sys*v(is,ig))                                 
            END DO
         END DO

       addition=0.d0
       do is=1,ns
       addition = addition + num(is)/den(is)
       end do       
       
     ct(i_mtd)= (1.d0/beta_cv)*dlog(addition)

     write(20,'(F16.6)')  ct(i_mtd)
     END DO

allocate(vbias(md_steps))
allocate(dum(ns))
close (20)
print *, 'done calculating ct. => ct.dat'

     DO i_md=1,md_steps
        mtd_max=  (i_md-1)*w_cv/w_hill                
        ss(1:ns)=cv(1:ns,i_md)
         dum=0.d0
        DO i_mtd=1,mtd_max
           ds2(1:ns)=width(1:ns,i_mtd)*width(1:ns,i_mtd)
           hh(1:ns)=ht(1:ns,i_mtd)/alpha_sys                                          
           diff_s2(1:ns)=ss(1:ns)-hill(1:ns,i_mtd)
            
              if (diff_s2(is) .gt. pi ) diff_s2(is) =diff_s2(is) - 2.d0*pi
              if (diff_s2(is) .lt.-pi ) diff_s2(is) =diff_s2(is) + 2.d0*pi
            
           diff_s2(1:ns)=diff_s2(1:ns)*diff_s2(1:ns)*0.5D0
           dum(1:ns)=dum(1:ns)+hh(1:ns)*DEXP(-diff_s2(1:ns)/ds2(1:ns))
        END DO
         vbias(i_md)= -(1.d0/beta_cv)*(dlog( dexp(-beta_cv*dum(1)) + dexp(-beta_cv*dum(2)) )- dlog(dfloat(ns)))

         write(30,'(I10,F16.6)') i_md, vbias(i_md)
!                                  vbias(i_md)+(1.d0/beta_cv)*dlog( dexp(-beta_cv*0.d0) + dexp(-beta_cv*0.d0) )

      END DO
close (30)
print *, 'done calculating pb_bias => vbias.dat'

allocate(prob(nbin(1),nbin(2)))
norm=0.0
prob=0.0
print *, 'calculating unbiased probabilty density => free energy.dat'

DO i_md=1,md_steps
        IF((i_md.GE.t_min).AND.(i_md.LT.t_max))THEN                                      ! t_min = 1 not 0
          index1 = nint((cv(1,i_md)-grid_min(1))/grid_width(1))+1
          index2 = nint((cv(2,i_md)-grid_min(2))/grid_width(2))+1
          if(index1.gt.0.and.index2.gt.0.and.index1.le.nbin(1).and.index2.le.nbin(2))then
             i_mtd =  (i_md-1)*w_cv/w_hill 
             dummy = vbias(i_mtd) + ct(i_mtd)             
             rweight = dexp(beta_cv*dummy)
             prob(index1,index2) = prob(index1,index2) + rweight
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
    s1= grid_min(1) + float(i_s1-1)*grid_width(1)
       DO i_s2=1,nbin(2)
       s2= grid_min(2) + float(i_s2-1)*grid_width(2)
       prob(i_s1,i_s2) = prob(i_s1,i_s2)*(1.d0/norm)                                    
       WRITE(50,'(4F16.8)') s1,s2,prob(i_s1,i_s2), -(1.d0/beta_cv)*dlog(max(prob(i_s1,i_s2),1e-16))*kj_to_kcal     !in kcal/mol.
       END DO
    write(50,*)
    END DO
print *, 'free_energy.dat => S_1, S_2, prob_density, free_energy'
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
