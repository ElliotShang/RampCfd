subroutine updateresidual
    use global
    implicit real(8)(a-h,o-z)
    real(8)::lambda1,lambda2,lambda3
    real(8),allocatable :: fx(:,:,:),fy(:,:,:)
    allocate (fx(4,in,jn-1),fy(4,in-1,jn))
    error=0.1d0
  ! x 方向
    do i=1,in
      do j=1,jn-1
        roul=rou(i-1,j)
        rour=rou(i,j)
        ul=u(i-1,j)
        ur=u(i,j)
        vL=v(i-1,j)
        vr=v(i,j)
        pL=p(i-1,j)
        pr=p(i,j)
       
        unl=ul*x1(i,j)/z1(i,j)+vL*y1(i,j)/z1(i,j)
        unr=ur*x1(i,j)/z1(i,j)+vr*y1(i,j)/z1(i,j)
        hL=gama2*pL/roul+0.5d0*(ul**2.d0+vL**2.d0)
        hr=gama2*pr/rour+0.5d0*(ur**2.d0+vr**2.d0)
        roubar=sqrt(roul*rour)
        ubar=(ul*sqrt(roul)+ur*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        unbar=(unl*sqrt(roul)+unr*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        vBar=(vL*sqrt(roul)+vr*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        hBar=(hL*sqrt(roul)+hr*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        c=sqrt((gama-1.d0)*(hBar-0.5d0*(ubar**2.d0+vBar**2.d0))) 
        drou=rour-roul
        dp=pr-pL
        du=ur-ul
        du=ur-ul
        dv=vr-vL
        dun=unr-unl
        lambda1=abs(ubar)
        lambda2=abs(ubar-c)
        lambda3=abs(ubar+c)
  
  ! harten修正
        if(abs(lambda1).lt.error) then
            lambda1=(lambda1**2.d0+error**2.d0)/(2.d0*error)
        end if
        if(abs(lambda2).lt.error) then
            lambda2=(lambda2**2.d0+error**2.d0)/(2.d0*error)
        end if     
        if(abs(lambda3).lt.error) then
            lambda3=(lambda3**2.d0+error**2.d0)/(2.d0*error)
        end if     
        
        beta1=z1(i,j)*abs(lambda1)*(drou-dp/(c**2.d0)) 
        beta2=z1(i,j)*abs(lambda3)*(dp+roubar*c*dun)/(2.d0*c**2.d0)       
        beta3=z1(i,j)*abs(lambda2)*(dp-roubar*c*dun)/(2.d0*c**2.d0)
        beta4=beta1+beta2+beta3
        beta5=c*(beta2-beta3)
        beta6=z1(i,j)*abs(lambda1)*roubar*(du-(x1(i,j)/z1(i,j))*dun) 
        beta7=z1(i,j)*abs(lambda1)*roubar*(dv-(y1(i,j)/z1(i,j))*dun)
  ! x方向通量计算            
        tempL1=z1(i,j)*roul*unl
        tempL2=z1(i,j)*roul*ul*unl+pL*x1(i,j)
        tempL3=z1(i,j)*roul*vL*unl+pL*y1(i,j)
        tempL4=z1(i,j)*roul*hL*unl
        tempr1=z1(i,j)*rour*unr
        tempr2=z1(i,j)*rour*ur*unr+pr*x1(i,j)
        tempr3=z1(i,j)*rour*vr*unr+pr*y1(i,j)
        tempr4=z1(i,j)*rour*hr*unr
        fx(1,i,j)=0.5d0*(tempL1+tempr1-beta4)
        fx(2,i,j)=0.5d0*(tempL2+tempr2-ubar*beta4-(x1(i,j)/z1(i,j))*beta5-beta6)
        fx(3,i,j)=0.5d0*(tempL3+tempr3-vBar*beta4-(y1(i,j)/z1(i,j))*beta5-beta7)
        fx(4,i,j)=0.5d0*(tempL4+tempr4-hBar*beta4-unbar*beta5-ubar*beta6-vBar*beta7+c*c*beta1/(gama-1.d0))
     end do
    end do           
  
  ! y 方向
    do j=1,jn
      do i=1,in-1
        roul=rou(i,j-1)
        rour=rou(i,j)
        ul=u(i,j-1)
        ur=u(i,j)
        vL=v(i,j-1)
        vr=v(i,j)
        pL=p(i,j-1)
        pr=p(i,j)
       
        unl=ul*x2(i,j)/z2(i,j)+vL*y2(i,j)/z2(i,j)
        unr=ur*x2(i,j)/z2(i,j)+vr*y2(i,j)/z2(i,j)
        hL=gama2*pL/roul+0.5d0*(ul**2.d0+vL**2.d0)
        hr=gama2*pr/rour+0.5d0*(ur**2.d0+vr**2.d0)
        roubar=sqrt(roul*rour)
        ubar=(ul*sqrt(roul)+ur*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        unbar=(unl*sqrt(roul)+unr*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        vBar=(vL*sqrt(roul)+vr*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        hBar=(hL*sqrt(roul)+hr*sqrt(rour))/(sqrt(roul)+sqrt(rour))
        c=sqrt((gama-1.d0)*(hBar-0.5d0*(ubar**2.d0+vBar**2.d0)))
        drou=rour-roul
        dp=pr-pL
        du=ur-ul
        du=ur-ul
        dv=vr-vL
        dun=unr-unl
        lambda1=abs(ubar)
        lambda2=abs(ubar-c)
        lambda3=abs(ubar+c)
  
  ! harten修正
        if(abs(lambda1).lt.error) then
            lambda1=(lambda1**2.d0+error**2.d0)/(2.d0*error)
        end if
        if(abs(lambda2).lt.error) then
            lambda2=(lambda2**2.d0+error**2.d0)/(2.d0*error)
        end if     
        if(abs(lambda3).lt.error) then
            lambda3=(lambda3**2.d0+error**2.d0)/(2.d0*error)
        end if     
        
        beta1=z2(i,j)*abs(lambda1)*(drou-dp/(c**2.d0)) 
        beta2=z2(i,j)*abs(lambda3)*(dp+roubar*c*dun)/(2.d0*c**2.d0)       
        beta3=z2(i,j)*abs(lambda2)*(dp-roubar*c*dun)/(2.d0*c**2.d0)
        beta4=beta1+beta2+beta3
        beta5=c*(beta2-beta3)
        beta6=z2(i,j)*abs(lambda1)*roubar*(du-(x2(i,j)/z2(i,j))*dun) 
        beta7=z2(i,j)*abs(lambda1)*roubar*(dv-(y2(i,j)/z2(i,j))*dun)
  ! y方向通量计算            
        tempL1=z2(i,j)*roul*unl
        tempL2=z2(i,j)*roul*unl*ul+pL*x2(i,j)
        tempL3=z2(i,j)*roul*vL*unl+pL*y2(i,j)
        tempL4=z2(i,j)*roul*hL*unl
        tempr1=z2(i,j)*rour*unr
        tempr2=z2(i,j)*rour*ur*unr+pr*x2(i,j)
        tempr3=z2(i,j)*rour*vr*unr+pr*y2(i,j)
        tempr4=z2(i,j)*rour*hr*unr
        fy(1,i,j)=0.5d0*(tempL1+tempr1-beta4)
        fy(2,i,j)=0.5d0*(tempL2+tempr2-ubar*beta4-x2(i,j)/z2(i,j)*beta5-beta6)
        fy(3,i,j)=0.5d0*(tempL3+tempr3-vBar*beta4-y2(i,j)/z2(i,j)*beta5-beta7)
        fy(4,i,j)=0.5d0*(tempL4+tempr4-hBar*beta4-unbar*beta5-ubar*beta6-vBar*beta7+c*c*beta1/(gama-1.d0))
     end do
    end do
  ! flux
    do K=1,4
      do i=1,in-1
          do j=1,jn-1
              flux(K,i,j)=fx(K,i+1,j)-fx(K,i,j)+fy(K,i,j+1)-fy(K,i,j) 
          end do
      end do
    end do 
    deallocate(fx,fy)               
  
    return
  end