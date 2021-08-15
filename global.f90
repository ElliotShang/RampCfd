module global
    integer in,jn,n,m
    real(8) p1,d1,vo,cfl,pouts,gama,r,gama1,gama2,gama3
    real(8), allocatable :: x(:,:), y(:,:),x1(:,:),y1(:,:),z1(:,:),x2(:,:),y2(:,:),z2(:,:),area(:,:)
    real(8), allocatable :: rou(:,:), u(:,:), v(:,:), p(:,:), flux(:,:,:) 
end module
  
  