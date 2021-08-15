subroutine update_solution
    use global
    implicit real(8)(a-h,o-z)
    dimension resave(4)
    allocatable :: qcon(:,:,:)
    save qcon
    if(.not.allocated(qcon)) then
      allocate( qcon(4,in-1,jn-1) )
      do 11 j=1,jn-1
      do 11 i=1,in-1
        qcon(1,i,j) = rou(i,j)
        qcon(2,i,j) = rou(i,j)*u(i,j)
        qcon(3,i,j) = rou(i,j)*v(i,j)
        qcon(4,i,j) = p(i,j)/gama1 + 0.5d0*rou(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
  11      continue
      open(70,file='monitor.dat')
    end if
    resm1 = 0.d0
    resave = 0.d0
    imax = 0
    jmax = 0
    do 21 j=1,jn-1
    do 21 i=1,in-1
  ! calculate the local time step:
      tx1 = 0.5d0*( x1(i,j) + x1(i+1,j) )
      ty1 = 0.5d0*( y1(i,j) + y1(i+1,j) )
      tz1 = 0.5d0*( z1(i,j) + z1(i+1,j) )
      tx2 = 0.5d0*( x2(i,j) + x2(i,j+1) )
      ty2 = 0.5d0*( y2(i,j) + y2(i,j+1) )
      tz2 = 0.5d0*( z2(i,j) + z2(i,j+1) )
      vnmci = tx1*u(i,j) + ty1*v(i,j)
      vnmcj = tx2*u(i,j) + ty2*v(i,j)
      sonic = sqrt(gama*p(i,j)/rou(i,j))
      chveL = abs(vnmci) + abs(vnmcj) + sonic*(tz1+tz2)
      dtcon = cfL/chveL
  ! update the flow variables:
      qcon(1:4,i,j) = qcon(1:4,i,j) - dtcon*flux(1:4,i,j)
      rou(i,j) = qcon(1,i,j)
      u(i,j) = qcon(2,i,j)/qcon(1,i,j)
      v(i,j) = qcon(3,i,j)/qcon(1,i,j)
      p(i,j) = gama1*( qcon(4,i,j) - 0.5d0*rou(i,j)*(u(i,j)*u(i,j)+v(i,j)*v(i,j)) )
  
      tres1 = abs(flux(1,i,j))
      if(tres1>resm1) then
        resm1 = tres1
        imax = i
        jmax = j
      end if
      resave(1:4) = resave(1:4) + abs(flux(1:4,i,j))
  21    continue
    resave = resave/real(in-1)/real(jn-1)
  
    if(mod(n,10)==0) write(70,7001) n, imax, jmax, resave(1:4)
    if(mod(n,10)==0) write(* ,7001) n, imax, jmax, resave(1:4)
  7001  format(3(i6),10(1x,e11.4))
    if(n==m) close(70)
  
    return
  end