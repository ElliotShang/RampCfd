subroutine updateboundary
    use global
    implicit real(8)(a-h,o-z)
  
  ! inlet
    temp=1.d0+gama3*vo**2
    p0=p1*temp**gama2
    ts0=p1/d1/r
    tt0=ts0*temp
    u0=vo*sqrt(gama*r*ts0)
    do 11 j=1,jn-1
      rou(0,j)=d1
      u(0,j)=u0
      v(0,j)=0.d0
      p(0,j)=p1
  11   continue
  ! outlet:
    do 21 j=1,jn-1
    vmchb=u(in-1,j)/sqrt(gama*p(in-1,j)/rou(in-1,j))
    p(in,j) = p(in-1,j)
    rou(in,j) = rou(in-1,j)
    u(in,j) = u(in-1,j)
    v(in,j) = v(in-1,j)
  21  continue
  
  ! the upper side 
    do 31 i=1,in-1
      vnormal=u(i,jn-1)*x2(i,jn)+v(i,jn-1)*y2(i,jn)
      vtemp=2.d0*vnormal/z2(i,jn)/z2(i,jn)
      rou(i,jn)=rou(i,jn-1)
      u(i,jn)=u(i,jn-1)-vtemp*x2(i,jn)
      v(i,jn)=v(i,jn-1)-vtemp*y2(i,jn)
      p(i,jn)=p(i,jn-1)
  31    continue
  
  ! the lower side 
    do 41 i=1,in-1
      vnormal=u(i,1)*x2(i,1)+v(i,1)*y2(i,1)
      vtemp=2.d0*vnormal/z2(i,1)/z2(i,1)
      rou(i,0)=rou(i,1)
      u(i,0)=u(i,1)-vtemp*x2(i,1)
      v(i,0)=v(i,1)-vtemp*y2(i,1)
      p(i,0)=p(i,1)
  41    continue
  
    return
end