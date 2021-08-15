subroutine initalizeflow
    use global
    implicit real(8)(a-h,o-z)
  
    dd = d1
    pp = p1
    uu = vo*sqrt(gama*p1/d1)
    do 11 j=1,jn-1
    do 11 i=1,in-1
      rou(i,j) = dd
      u(i,j) = uu
      v(i,j) = 0.d0
      p(i,j) = pp
  11  continue
  
    return
end