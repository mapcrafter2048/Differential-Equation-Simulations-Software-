Program model
Implicit none 
integer, parameter::tran=20000,n=1000,nite=30000   !!!change for number of oscillators
real,dimension(n) :: omega,theta,dth,tho, dth_o
real:: pi,t,r1,rx1,ry1,shi1,r,k1,alpha
real:: A,var,rx2,ry2,r_2,r2,shi2
real, parameter :: h=0.05d0,k2=8.0d0,eta=1.0d0
integer :: i,j,it, k_loop
pi=4.0*atan(1.0)

open(unit = 1,file="/home/priyanka/Downloads/final_model/rt.txt")

!=================frequency===========
alpha=1
do i=1,n
   !omega(i)=-1+1*rand()
   omega(i)=alpha*tan((i*pi)/float(n)-((n+1)*pi)/float(2*n))
 ! omega(i)=alpha*tan((pi/2)*(2*i-n-1)/(n+1))
 !var=rand();
 !omega(i)=tan(pi*(var-0.5))
  ! write(*,*) omega(i)
enddo 
A=n*eta 
do i=1,n
!if (i.le.A) then
  !theta(i)=0
  theta(i)=-pi+2*pi*rand()
!else
 !theta(i)=2*pi
!endif    
enddo

!=======transient loop===============
k1=0
do k_loop=1,80
r=0.0d0;
r_2=0.0d0;
t=0.0d0;

do it=1,nite
!======order parameters=======================   
  rx1=0.0
  ry1=0.0
  rx2=0.0
  ry2=0.0
	do i=1,n
		rx1= rx1+ COS(theta(i))
		ry1= ry1+ SIN(theta(i))
	        rx2= rx2+ COS(2*theta(i))
		ry2= ry2+ SIN(2*theta(i))
	enddo
	r1=SQRT(rx1*rx1 + ry1*ry1)/float(n)
	r2=SQRT(rx2*rx2 + ry2*ry2)/float(n)
	shi1 =  atan(ry1/rx1)
	shi2 =  atan(ry2/rx2)
	if (it>tran) then
	r=r+r1
	r_2=r_2+r2
	endif
	!write(14,*) t,theta(1),theta(2)   
	!write(16,*)t 
	
        call derivs(t,theta, dth, k1,k2,omega,r1,shi1,r2,shi2)
        call rk4(theta,dth,n,t,h,tho,derivs,k1,k2, omega,r1,shi1,r2,shi2)
        do i=1,n
             theta(i)=MOD(tho(i),2*pi)
        enddo
        do i=1,n
            dth_o(i)=dth(i)
        enddo    
   t=t+h
!   write(1,*) t,r1
!	write(*,*) t,r1
enddo
   r=r/real(nite-tran)
  r_2=r_2/real(nite-tran)

 !  do i=1,n
  ! write(1,*) omega(i), theta(i)
  ! write(*,*) omega(i),  theta(i)
  ! enddo
  write(1,*) real(k1), r ,r_2
  write(*,*) real(k1), r,r_2
  k1=k1+0.5
enddo 
   
stop       
end program model



!==============code for equations======================
Subroutine derivs(t,theta, dth, k1,k2, omega,r1,shi1,r2,shi2)
implicit none
integer, parameter :: n=1000                            !!!change for number of oscillators
real,dimension(n):: dth1,dth2,dth,omega,theta
real:: t,r1,shi1,k1,k2,r2,shi2
integer:: i
real, parameter :: alpha1=0.0d0, beta1=0.0d0                  !!!!!!!!!!!!!change here for alpha and beta!!!!!!!!!!    
                                                !!!!!!next calculate r value here

   do i=1,n
        dth1(i)=(k1*(r1**(alpha1+1))*sin(shi1-theta(i)))
        dth2(i)=(k2*(r2*r1**(beta1+1))*sin(shi2-shi1-theta(i)))
        dth(i)=omega(i)+dth1(i)+dth2(i)
   enddo
 return
 end subroutine
 
 
 
!=========Rk4 code============================= 
SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs,k1,k2,omega,r1,shi1,r2,shi2)
Implicit none
INTEGER:: n
REAL,intent(in)::h,x,dydx(n),y(n), omega(n)
real,intent(out)::yout(n)
EXTERNAL derivs
integer,PARAMETER::NMAX=100000 
real::dym(NMAX),dyt(NMAX),yt(NMAX), k1, k2,r1, shi1,r2,shi2
real:: xh, hh, h6
integer:: i
hh=h*0.5d0
h6=h/6.d0
xh=x+hh
do i=1,n
  yt(i)=y(i)+hh*dydx(i)
enddo 
call derivs(xh,yt,dyt, k1, k2,omega,r1,shi1,r2,shi2)
do i=1,n
  yt(i)=y(i)+hh*dyt(i)
enddo 
call derivs(xh,yt,dym, k1,k2, omega,r1,shi1,r2,shi2)
do i=1,n
 yt(i)=y(i)+h*dym(i)
 dym(i)=dyt(i)+dym(i)
enddo 
call derivs(x+h,yt,dyt, k1, k2,omega,r1,shi1,r2,shi2)
do i=1,n
 yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
enddo 
END subroutine rk4






