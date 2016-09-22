C	Subroutines in this file include:

c-----------------------------------------------------------------------
!	subroutine clebgord(n,n_,m,Cn10)
	!Calculates the Clebsch-Gordan coefficients as required in the
	!orientational averaging of the T_matrix (equation 3.27
	!Mishchenko, 1991)

	!subroutine legendrecoeff(nmax,pn)
	!calculate the Legendre expansion coefficients pn for the Perfect
	!Davis-Greenstein alignment of elongated particles up to n = nmax
	!as required in orientationally averaging the T-matrix  (equation 3.27
	!Mishchenko, 1991)

!	subroutine avgTmatrix(NMAX)
	
	!Performs the orientation averaging of equation 3.27 in Mishchenko
	!1991
	
	
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	subroutine clebgord(n,n_,m,NMAX,Cn10)
	!Calculates the Clebsch-Gordan coefficients as required in the
	!orientational averaging of the T_matrix (equation 3.27
	!Mishchenko, 1991) only works for n>=n_.
	!For n<n_ clebgord(n,n_,m,Cn10) = clebgord(n_,n,m,Cn10)
	!
	!-------------------------------------------------------------
	implicit none
	integer n,n_,m
	integer n1,ndummy,NMAX
	real Cn10(0:2*NMAX)
	real Cnn_m,a,b,c
	!-------------------------------------------------------------
	!if n < n_, apply recipricity relation
	!if (n<n_) then
	!	ndummy=n
	!	n=n_
	!	n_=ndummy
	!endif
	
	!calculate Cnn_m (B4)
	Cnn_m=(-1.0)**(n_+m)*(2.0*n_+1.0)**(-0.5)
	!print *, "Cnn_m = ", Cnn_m
	do ndummy=n_+1,n
	  Cnn_m=Cnn_m*((ndummy+m)*(ndummy-m)*(2.0*ndummy-2.0*n_+1.0)/
     &	ndummy/(2.0*ndummy+1.0)/(ndummy-n_))**0.5
		!print *, "Cnn_m = ", Cnn_m
	end do
	!calculate first 2 elements
	Cn10(n-n_)=Cnn_m
	n1=n-n_+1
	Cn10(n1)=(4.0*(2.0*n1+1.0)*(2.0*n1-1.0)/(n-n_+n1)/(-n+n_+
     &	n1)/(n+n_-n1+1.0)/(n+n_+n1+1.0))**0.5*
     &	(m*Cn10(n1-1)-0.0)
     	!print *, "Cn10(n1) = ", Cn10(n1),n1
	!now the rest
	do n1=n-n_+2,n+n_
		a=(4.0*(2.0*n1+1.0)*(2.0*n1-1.0)/(n-n_+n1)/(-n+n_+
     &		n1)/(n+n_-n1+1.0)/(n+n_+n1+1.0))**0.5
     		!print *,"a = ",a
		b=m*Cn10(n1-1)
		!print *,"b = ",b
		c=((-n+n_+n1-1.0)*(n-n_+n1-1.0)*(n+n_-n1+2.0)*(n+n_+n1)/4.0/
     &	 (2.0*n1-3.0)/(2.0*n1-1.0))**0.5*Cn10(n1-2)
     		!print *,"c = ",c
		Cn10(n1)=a*(b-c)
      	!print *, "Cn10(n1) = ", Cn10(n1)
	
	end do
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	subroutine legendrecoeff(nmax,pn)
	!calculate the Legendre expansion coefficients pn for the Perfect
	!Davis-Greenstein alignment of elongated particles up to n = nmax
	!as required in orientationally averaging the T-matrix  (equation 3.27
	!Mishchenko, 1991)
	!----------------------------------------------------------------
	implicit none
	integer nmax
	integer i
	real pn(0:nmax)
	
	pn(0)=1.0
	do i = 1,nmax
		if (mod(i,2)==1) then
			pn(i)=0
		else
			pn(i)=-(i-1.0)/i*pn(i-2)
		end if
	end do
	
	return
	end	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
C********************************************************************
	subroutine avgTmatrix(NMAX)
	
	!Performs the orientation averaging of equation 3.27 in Mishchenko
	!1991
	
	INCLUDE 'amplq.par.f'
	
	integer NMAX
	integer m,m_,n,n_,n1,m1,i,j
	REAL*4
     &     RT11(NPN6,NPN4,NPN4),RT12(NPN6,NPN4,NPN4),
     &     RT21(NPN6,NPN4,NPN4),RT22(NPN6,NPN4,NPN4),
     &     IT11(NPN6,NPN4,NPN4),IT12(NPN6,NPN4,NPN4),
     &     IT21(NPN6,NPN4,NPN4),IT22(NPN6,NPN4,NPN4)
	
	REAL*4
     &     newRT11(NPN6,NPN4,NPN4),newRT12(NPN6,NPN4,NPN4),
     &     newRT21(NPN6,NPN4,NPN4),newRT22(NPN6,NPN4,NPN4),
     &     newIT11(NPN6,NPN4,NPN4),newIT12(NPN6,NPN4,NPN4),
     &     newIT21(NPN6,NPN4,NPN4),newIT22(NPN6,NPN4,NPN4)
	real::Cn10(0:2*NMAX),Cn10sum2(0:2*NMAX),Cn10sum2m0(0:2*NMAX)
	real::pn(0:2*NMAX)
	real Resum1(2,2),Imsum1(2,2),Resum2(2,2),Imsum2(2,2),ReTij(2,2),
     &	ImTij(2,2)
	
	COMMON /TMAT/ RT11,RT12,RT21,RT22,IT11,IT12,IT21,IT22
	
	!get vector of pns
	call legendrecoeff(2*NMAX,pn)
	
	DO m_=0,NMAX 
	  m=m_+1!!An annoying feature.  m_ is m in the equations
	  	!!m are the indices for RTIJ etc
	  DO n=Max(m_,1),NMAX
	  	DO n_=Max(m_,1),NMAX
		  !calculate new Tmatrix componenets T^ij_mnn1nn2
		  !first summation in eq. 3.27
		  Resum1(1,1)=0
		  Resum1(1,2)=0
		  Resum1(2,1)=0
		  Resum1(2,2)=0
		  Imsum1(1,1)=0
		  Imsum1(1,2)=0
		  Imsum1(2,1)=0
		  Imsum1(2,2)=0
c		  allocate (Cn10(abs(n-n_):n+n_))
c		  allocate (Cn10sum2(abs(n-n_):n+n_))
c		  allocate (Cn10sum2m0(abs(n-n_):n+n_))
		  
		  if (n<n_) then
		  	!use reciprocity
		  	call clebgord(n_,n,m_,NMAX,Cn10)
		      call clebgord(n_,n,0,NMAX,Cn10sum2m0)
		  else
		  	call clebgord(n,n_,m_,NMAX,Cn10)
		  	call clebgord(n,n_,0,NMAX,Cn10sum2m0)
		  end if
		  do n1=abs(n-n_),n+n_
			!second summation in eq. 3.27
			!put m1=0 calculation here
			ReTij(1,1)=RT11(1,n,n_)
			ReTij(1,2)=RT12(1,n,n_)
			ReTij(2,1)=RT21(1,n,n_)
			ReTij(2,2)=RT22(1,n,n_)
			ImTij(1,1)=IT11(1,n,n_)
			ImTij(1,2)=IT12(1,n,n_)
			ImTij(2,1)=IT21(1,n,n_)
			ImTij(2,2)=IT22(1,n,n_)
			do i=1,2
			   do j=1,2
			   Resum2(i,j)=0.5*Cn10sum2m0(n1)*ReTij(i,j)
			   Imsum2(i,j)=0.5*Cn10sum2m0(n1)*ImTij(i,j)
			   end do
			end do
			do m1=1,min(n,n_)!m1 is the m1 in equations
			  ReTij(1,1)=RT11(m1+1,n,n_)
			  ReTij(1,2)=RT12(m1+1,n,n_)
			  ReTij(2,1)=RT21(m1+1,n,n_)
			  ReTij(2,2)=RT22(m1+1,n,n_)
			  ImTij(1,1)=IT11(m1+1,n,n_)
			  ImTij(1,2)=IT12(m1+1,n,n_)
			  ImTij(2,1)=IT21(m1+1,n,n_)
			  ImTij(2,2)=IT22(m1+1,n,n_)
			  !calculate C^n10_nm1n_-m
			  if (n<n_) then
		  		!use reciprocity
				call clebgord(n_,n,m1,NMAX,Cn10sum2)
			  else
				call clebgord(n,n_,m1,NMAX,Cn10sum2)
			  end if
			  do i=1,2
			   do j=1,2
			  	Resum2(i,j)=Resum2(i,j)+(-1.0)**m1*Cn10sum2(n1)*
     &				 ReTij(i,j)
			  	Imsum2(i,j)=Imsum2(i,j)+(-1.0)**m1*Cn10sum2(n1)*
     &				ImTij(i,j)
     			   end do
			  end do
			end do
			do i=1,2
			  do j=1,2
			   Resum1(i,j)=Resum1(i,j)+(1.0+(-1.0)**(n+n_+n1+i+j))
     &			    *pn(n1)*Cn10(n1)*Resum2(i,j)
     			   Imsum1(i,j)=Imsum1(i,j)+(1.0+(-1.0)**(n+n_+n1+i+j))
     &			    *pn(n1)*Cn10(n1)*Imsum2(i,j)
     			    
     			  end do
			end do
		  end do
		  newRT11(m,n,n_)=(-1.0)**m_*Resum1(1,1)
		  newRT12(m,n,n_)=(-1.0)**m_*Resum1(1,2)
		  newRT21(m,n,n_)=(-1.0)**m_*Resum1(2,1)
		  newRT22(m,n,n_)=(-1.0)**m_*Resum1(2,2)
		  newIT11(m,n,n_)=(-1.0)**m_*Imsum1(1,1)
		  newIT12(m,n,n_)=(-1.0)**m_*Imsum1(1,2)
		  newIT21(m,n,n_)=(-1.0)**m_*Imsum1(2,1)
		  newIT22(m,n,n_)=(-1.0)**m_*Imsum1(2,2)
c		  deallocate (Cn10)
c		  deallocate (Cn10sum2)
c		  deallocate (Cn10sum2m0)
		end do	
	  end do
	end do
	DO m_=0,NMAX 
	  m=m_+1!!An annoying feature.  m_ is m in the equations
	  	!!m are the indices for RTIJ etc
	  DO n=Max(m_,1),NMAX
	  	DO n_=Max(m_,1),NMAX
		RT11(m,n,n_)=newRT11(m,n,n_)
		RT12(m,n,n_)=newRT12(m,n,n_)
		RT21(m,n,n_)=newRT21(m,n,n_)
		RT22(m,n,n_)=newRT22(m,n,n_)
		IT11(m,n,n_)=newIT11(m,n,n_)
		IT12(m,n,n_)=newIT12(m,n,n_)
		IT21(m,n,n_)=newIT21(m,n,n_)
		IT22(m,n,n_)=newIT22(m,n,n_)
		end do
	  end do
	end do  
	return
	end		
