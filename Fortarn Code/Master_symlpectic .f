c	symplectic all to all, all to all multilayer
c	For all to all, based on mean-field eq.
c	for any query: ask sarika

c	nstep: total time, itrans: initial transient: change as per need
c	h: rk4 step size
c	second layer has alpha. change lambda2 (43rd line) as per as need, and Dx (43rd line) multiplexing strength

	integer ndim,ilambda,nlines1,nlines2
	parameter(ndim=1000)

	real omega1(ndim),omega2(ndim),omega_mean1,omega_mean2            !! parameters to for intrinsic frequency 
	real theta1in(ndim),theta2in(ndim),alpha,alpha2
	real theta1out(ndim),theta2out(ndim),dydt1(ndim),dydt2(ndim)
	integer a2(ndim,ndim), a1(ndim,ndim)
	real h,t,r1_final,r2_final,shi1,shi2,lambda1_step	

	real lambda1,lambda2,Dx,coup1,coup2,add1,add2,lambda1_max        !!parameters for coupling strength
	real rx1,rx2,ry1,ry2,r1,r2,x
	real rho1,rho2,phi1,phi2,lambda3,lambda1_min
	real rhox1,rhox2,rhoy1,rhoy2,rho1_final,rho2_final

	real pi
	external rk4,derivs
	data h/0.01/, nstep/20000/, t/0.0/
	itrans=10000
10	format(5x(f10.6))

	pi=4.0*atan(1.0)

c  ==== !!!Chnage these parameters according to differential equation=================
	val=1
	a=1.0
	b=0.0
	lambda2=8.0
	lambda3=0.0
	lambda1_step=0.2
	lambda1_max = 3.0
	lambda1_min= 1.0

c=============== Lorenzian freq============
	alpha=1.0
	 do i=1,ndim
             omega1(i)=alpha*tan( (i*pi)/float(ndim)
     *      -  (float(ndim+1)*pi)/float(2*ndim) )
             omega2(i)=alpha*tan( (i*pi)/float(ndim)
     *      -  (float(ndim+1)*pi)/float(2*ndim) )
c   random lorenzian distribution----
		!x=rand()
	        !if(x .le. 0.0001)then
	         !  x=0.1
	        !endif
c	    omega1(i)=tan( pi*(rand()-0.5d0) )
c	    omega2(i)=tan( pi*(rand()-0.5d0) )
c ============================================
c              write(80,*)i,omega1(i)
c             write(*,*)i,x,omega1(i)
              omega_mean1 = omega_mean2 + omega1(i) 
              omega_mean2 = omega_mean2 + omega2(i) 
           enddo
	
	 omega_mean1 = omega_mean1/float(ndim)
	 omega_mean2 = omega_mean2/float(ndim)
	! write(*,*)"omega1_mean",omega_mean1
	! write(*,*)"omega2_mean",omega_mean2

c ====homogenous freq distribution=============
	   !do i=1,ndim
	    ! omega1(i)=-0.5+1.0*rand()
	     !omega2(i)=-0.5+1.0*rand()
           !enddo
c==========================================

	
	do i=1,ndim
	! val=rand()
	! if (val>0.8) then
            ! theta1in(i)=2*pi                           !!initial conditions for theta in backward direction
             theta1in(i)=-pi+2*pi*rand()                 !!initial conditions for theta in forward direction
         ! else
          !   theta1in(i)=2.0*pi
          !endif
             theta2in(i)=2*pi
        enddo
       
	do i=1,ndim
             !write(*,*)theta1in(i)
c             theta2in(i)=2*pi
        enddo
c=========================

c	min and max values of lambda 1 (layer without alpha)
	       !lambda1= lambda1_min -  lambda1_step 
	       lambda1=lambda1_min                              !!Change here for forward and backward 
		do while(lambda1 .le. lambda1_max)
	         lambda1=lambda1+lambda1_step
		!do while(lambda1 .le. lambda1_max)
	          !lambda1=lambda1+lambda1_step

	!do i=1,ndim
         !    theta1in(i)=2*pi*rand()
          !   theta2in(i)=2*pi*rand()
        !enddo
	   r1=0.0
	   r2=0.0
	   r1_final=0.0
	   r2_final=0.0
	   rho1=0.0
	   rho2=0.0
	   rho1_final=0.0
	   rho2_final=0.0

c	=======evolving in time==============
	do it=1,nstep
c	====Calculations for order parameters r1 and r2 for one layer and rho1 and rho2 for another layer==============
	rx1=0
	rx2=0
	ry1=0.0
	ry2=0.0
	rhox1=0
	rhox2=0
	rhoy1=0.0
	rhoy2=0.0
	do i=1,ndim
		rx1= rx1+ COS(theta1in(i))
		ry1= ry1+ SIN(theta1in(i))
		rx2= rx2+ COS(2.0*theta1in(i))
		ry2= ry2+ SIN(2.0*theta1in(i))
		rhox1= rhox1+ COS(theta2in(i))
		rhoy1= rhoy1+ SIN(theta2in(i))
		rhox2= rhox2+ COS(2.0*theta2in(i))
		rhoy2= rhoy2+ SIN(2.0*theta2in(i))
	enddo
	r1 = SQRT(rx1*rx1 + ry1*ry1)/float(ndim)
	r2 = SQRT(rx2*rx2 + ry2*ry2)/float(ndim)
	rho1 = SQRT(rhox1*rhox1 + rhoy1*rhoy1)/float(ndim)
	rho2 = SQRT(rhox2*rhox2 + rhoy2*rhoy2)/float(ndim)
	shi1 =  atan(ry1/rx1)
	shi2 =  atan(ry2/rx2)
	phi1 =  atan(rhoy1/rhox1)
	phi2 =  atan(rhoy2/rhox2)

	call rk4(alpha2,lambda2,dx,r1,r2,ndim,omega1,omega2,lambda1,
     *	 theta1in,theta2in,dydt1,dydt2,t,h,theta1out,theta2out,shi1,shi2,
     *	  rho1,rho2,phi1,phi2,lambda3) 
	 do i=1,ndim
	 	theta1in(i)=MOD(theta1out(i),2.0*pi)
	 	theta2in(i)=MOD(theta2out(i),2.0*pi)
	 enddo
	 

	if(it > itrans)then
	r1_final = r1 + r1_final
	r2_final = r2 + r2_final
	rho1_final = rho1 + rho1_final
	rho2_final = rho2 + rho2_final
	endif
	enddo  !nstep

	r1_final= r1_final/float(nstep-itrans)
	r2_final= r2_final/float(nstep-itrans)
	rho1_final= rho1_final/float(nstep-itrans)
	rho2_final= rho2_final/float(nstep-itrans)

	write(*,*)lambda1,r1_final,r2_final
	write(21,*)lambda1,r1_final,r2_final

	enddo !ilambda1
	

	!enddo !lambda2
	!enddo !dx
c	do i=1,ndim
c             write(11,*)theta1in(i)
c             write(12,*)theta2in(i)
c        enddo
	
	
	

	return
	end
c======================  rk-4   =============================
	Subroutine rk4(alpha2,lambda2,Dx,r1,r2,n,omega1,omega2,
     *        lambda1,y1,y2,dydt1,dydt2,t,h,yout1,yout2,shi1,shi2,
     *		rho1,rho2,phi1,phi2,lambda3)

	integer nmax,n
	real h,t,dydt1(n),y1(n),yout1(n)
	real dydt2(n),y2(n),yout2(n)
	parameter (nmax=1000)
	integer i
	real h6,hh,th,dym1(nmax),dyt1(nmax),yt1(nmax)
	real dym2(nmax),dyt2(nmax),yt2(nmax)
	real omega1(n),omega2(n),lambda1,lambda2,Dx,alpha2
	integer a2(n,n),a1(n,n)
	real r1,r2,shi1,shi2,lambda3
	real rho1,rho2,phi1,phi2
	external derivs

	hh = h*0.50
	h6 = h/6.0
	th = t + hh
	call derivs(alpha2,lambda2,Dx,r1,r2,n,t,omega1,omega2,
     *  lambda1,y1,y2,dydt1,dydt2,shi1,shi2,rho1,rho2,phi1,
     *    phi2,lambda3)
	  do 11 i=1,n
	     yt1(i) = y1(i) + hh*dydt1(i)
	     yt2(i) = y2(i) + hh*dydt2(i)
11        continue
	call derivs(alpha2,lambda2,Dx,r1,r2,n,th,omega1,omega2,
     *	 lambda1,yt1,yt2,dyt1,dyt2,shi1,shi2,rho1,rho2,phi1,phi2,
     *	lambda3)
	   do 12 i=1,n
	      yt1(i) = y1(i) + hh*dyt1(i)
	      yt2(i) = y2(i) + hh*dyt2(i)
12         continue
	call derivs(alpha2,lambda2,Dx,r1,r2,n,th,omega1,omega2,
     *	   lambda1,yt1,yt2,dym1,dym2,shi1,shi2,rho1,rho2,phi1,phi2,
     *	lambda3)
	  do 13 i =1,n
	    yt1(i) = y1(i) + h*dym1(i)
	    dym1(i) = dyt1(i) + dym1(i)
	    yt2(i) = y2(i) + h*dym2(i)
	    dym2(i) = dyt2(i) + dym2(i)
13        continue
	call derivs(alpha2,lambda2,Dx,r1,r2,n,t+h,omega1,omega2,
     *	lambda1,
     *	yt1,yt2,dyt1, dyt2,shi1,shi2,rho1,rho2,phi1,phi2,lambda3)
	  do 14 i =1,n
	     yout1(i) = y1(i) + h6*(dydt1(i)+dyt1(i)+2.0*dym1(i))
	     yout2(i) = y2(i) + h6*(dydt2(i)+dyt2(i)+2.0*dym2(i))
14	  continue
	return
	end

!!!!!!!!!!!!!     Differential Equation !!!!!!!!!!!!!!!!!!!!!!!
	subroutine derivs(alpha2,lambda2,Dx,r1,r2,ndim,t,omega1,
     *       omega2,lambda1,theta1in,theta2in,dydt1,dydt2,shi1,shi2,
     *		rho1,rho2,phi1,phi2,lambda3)
	integer ndim
	real omega1(ndim),omega2(ndim),theta1in(ndim),theta2in(ndim)
	real dydt1(ndim),dydt2(ndim),coup1,coup2,lambda2,Dx,lambda1,alpha2
	integer a2(ndim,ndim),a1(ndim,ndim),degree
	real r1,r2,shi1,shi2,sigma1,sigma2,sigma3
	real rho1,rho2,phi1,phi2,lambda3

c for all to LL
	sigma1=0.0
	sigma2=0.0
	sigma3=0.0


	Do i =1,ndim

	
	dydt1(i) = omega1(i) + lambda1*(r1)*sin(shi1-theta1in(i)) +
     *  	lambda2*r1*r2*sin(shi2-shi1-theta1in(i)) +
     *		lambda3*r1*r1*rho1*sin(shi1-theta1in(i))

	dydt2(i) = omega2(i) + sigma1*rho1*sin(phi1-theta2in(i)) +
     *  	sigma2*rho1*rho2*sin(phi2-phi1-theta2in(i)) +
     *		sigma3*rho1*rho1*rho1*sin(phi1-theta2in(i))

	enddo
c ========
	return
	end
