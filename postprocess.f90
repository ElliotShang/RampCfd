subroutine postprocess
    use global
    implicit real(8)(a-h,o-z)
    open(71,file='flowvar.dat')
    write(71,*) 'title="noname"'
    write(71,*) 'variables="x","y","rou","u","v","spre","tpre","ttem","mach"'
    write(71,*) 'zone t="noname", i=',in+1,', j=',jn+1,', f=point'
    do  j=1,jn
    do  i=1,in
      stem = p(i,j)/(r*rou(i,j))
      vmac = sqrt((u(i,j)*u(i,j)+v(i,j)*v(i,j))/(gama*p(i,j)/rou(i,j)))
      temp = 1.d0 + gama3*vmac*vmac
      ttem = stem*temp
      tpre = p(i,j)*temp**gama2
      write(71,7101) x(i,j),y(i,j),rou(i,j),u(i,j),v(i,j),p(i,j),tpre,ttem,vmac
    end do
    end do
    close(71)
  7101  format(20(1x,f18.8))
  
    return
  end