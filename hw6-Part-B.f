	program OverRelaxation
c..Taha Yaşar Demir / 1881978
c..CE580 - HomeWork #6 PART-B
	parameter(mx=101)
	common/mk_grid/   dx,dy,x(mx),y(mx),N
	common/analytic/  u(mx,mx)
	common/numeric/   s(mx,mx),w,dw,r(mx,mx)
	common/error/     e,old_e

	open(11,file="error.dat")

	print*, "Enter Omega Value"
	read*, w

	call grid
	call init
	do k=1,1000 ! solution iteration
		call solution
		call error_cal
		call output(k)
	enddo

	close(11)

	stop
	end
c-----------------------------------------------------------------------
	subroutine grid
	parameter(mx=101)
	common/mk_grid/   dx,dy,x(mx),y(mx),N
	common/analytic/  u(mx,mx)
	common/numeric/   s(mx,mx),w,dw,r(mx,mx)
	common/error/     e,old_e

	N    = 51
	dx   = 1./(N-1)
	dy   = 1./(N-1)
	x(1) = 0.
	y(1) = 0.

	do i=2,N
		x(i) = x(i-1) + dx
		y(i) = y(i-1) + dy
	enddo

	return
	end
c-----------------------------------------------------------------------
	subroutine init
	parameter(mx=101)
	common/mk_grid/   dx,dy,x(mx),y(mx),N
	common/analytic/  u(mx,mx)
	common/numeric/   s(mx,mx),w,dw,r(mx,mx)
	common/error/     e,old_e

	do i=1,N
		do j=1,N
			s(i,j) = 0. ! Initialize internal points
			u(i,j) = x(i)**2 - y(j)**2 ! Analytical solution
		enddo 
	enddo

	! boundary conditions
	do i=1,N
		s(i,1) = u(i,1)
		s(i,N) = u(i,N) 
	enddo
	do j=1,N
		s(1,j) = u(1,j)
		s(N,j) = u(N,j)
	enddo

	return
	end
c-----------------------------------------------------------------------
	subroutine solution
	parameter(mx=101)
	common/mk_grid/   dx,dy,x(mx),y(mx),N
	common/analytic/  u(mx,mx)
	common/numeric/   s(mx,mx),w,dw,r(mx,mx)
	common/error/     e,old_e

	do i=2,N-1
		do j=2,N-1
			r(i,j) = 0.25*(s(i+1,j) + s(i-1,j)+
     +				s(i,j+1) + s(i,j-1) - 4*s(i,j))
			s(i,j) = s(i,j) + w*r(i,j)
		enddo
	enddo

	return
	end
c-----------------------------------------------------------------------
	subroutine error_cal
	parameter(mx=101)
	common/mk_grid/   dx,dy,x(mx),y(mx),N
	common/analytic/  u(mx,mx)
	common/numeric/   s(mx,mx),w,dw,r(mx,mx)
	common/error/     e,old_e
	
	e = 0. ! e : Error

	do i=2,N-1
		do j=2,N-1
 			e = e + abs(u(i,j)-s(i,j))/((N-2)*(N-2)) 
		enddo
	enddo

	return
	end
c-----------------------------------------------------------------------
	subroutine output(iter)
	parameter(mx=101)
	common/mk_grid/   dx,dy,x(mx),y(mx),N
	common/analytic/  u(mx,mx)
	common/numeric/   s(mx,mx),w,dw,r(mx,mx)
	common/error/     e,old_e
	integer iter

	write(11,*) iter,e


	return
	end