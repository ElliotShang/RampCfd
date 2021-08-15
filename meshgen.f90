subroutine meshgen
    use global
    implicit real(8)(a-h,o-z)
  
  ! generate the mesh
    xinl=0.d0
    xout=3.d0
    do i=1,in
      do j=1,jn
      xx = xinl + real(i-1.)/real(in-1.)*(xout-xinl)
      if(xx<=0.6d0) then
        yi = 0.0d0+real(j-1.0d0)/real(jn-1.)*1.0d0
      else if(xx>=1.2d0) then
        yi = 0.20d0+real(j-1.0d0)/real(jn-1.)*(1.0d0-0.2d0)
      else
        ywall = ((xx-0.6d0)/0.6)*0.2
        yi = ywall+real(j-1.0d0)/real(jn-1.0d0)*(1.0d0-ywall)
      end if
      x(i,j) = xx
      y(i,j) = yi
    end do
    end do
  
    open(71,file='pltgrid.dat')
    write(71,*) 'title="noname"'
    write(71,*) 'variables="x","y"'
    write(71,*) 'zone t="noname", i=',in,',j=',jn,', f=point'
    do j=1,jn
      do i=1,in
        write(71,*) x(i,j), y(i,j)
      end do
    end do
    close(71)
                
    do 21 j=1,jn-1
    do 21 i=1,in
      x1(i,j) = y(i,j+1)-y(i,j)
      y1(i,j) = x(i,j)-x(i,j+1)
      z1(i,j) = sqrt(x1(i,j)*x1(i,j)+y1(i,j)*y1(i,j))
  21    continue
  
    do 22 j=1,jn
    do 22 i=1,in-1
      x2(i,j) = y(i,j)-y(i+1,j)
      y2(i,j) = x(i+1,j)-x(i,j)
      z2(i,j) = sqrt(x2(i,j)*x2(i,j)+y2(i,j)*y2(i,j))
  22    continue
  
    do 23 j=1,jn-1
    do 23 i=1,in-1
      ii = i + 1
      jj = j + 1
      area(i,j) = 0.5d0*abs((x(i,j)-x(ii,j))*y(ii,jj)+(x(ii,j)-x(ii,jj))*y(i,j)+(x(ii,jj)- x(i,j))*y(ii,j))
      area(i,j) = area(i,j)+0.5d0*abs((x(i,j)-x(ii,jj))*y(i,jj)+(x(ii,jj)-x(i,jj))*y(i,j)+(x(i,jj)-x(i,j))*y(ii,jj))
  23    continue
  
    return
end