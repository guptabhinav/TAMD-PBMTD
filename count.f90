program addition
real*8 :: dum,prob,sum

integer :: ios

open(1,file='probability_2D.dat')

sum=0.d0
 
do
 read(1,*,iostat=ios)dum,dum,prob
 if(ios.ne.0) exit
 sum=sum+prob
end do

 print *, ' probability sum =', sum
end
