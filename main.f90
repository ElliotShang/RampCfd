include 'global.f90'
include 'Initalizeflow.f90'
include 'meshgen.f90'
include 'boundarycon.f90'
include 'updateresidual.f90'
include 'update_solution.f90'
include 'postprocess.f90'
program main
    use global 
    implicit real(8)(a-h,o-z)     
    in = 400
    jn = 100 
    p1 = 31836.315d0
    d1 = 0.535815d0
    vo = 2.0d0
    m = 80000
    cfL = 0.3d0
    allocate(x(in,jn),y(in,jn))
    allocate(x1(in,jn-1),y1(in,jn-1),z1(in,jn-1),x2(in-1,jn),y2(in-1,jn),z2(in-1,jn),area(in-1,jn-1))
    allocate(rou(0:in,0:jn),u(0:in,0:jn),v(0:in,0:jn),p(0:in,0:jn),flux(4,in-1,jn-1))
    r  = 287.06d0
    gama  = 1.4d0
    gama1 = gama - 1.d0
    gama2 = gama/gama1
    gama3 = gama1/2.d0
    call meshgen
    call initalizeflow
    do n=1,m
      call updateboundary
      call updateresidual
      call update_solution
    end do
    call updateboundary
    call postprocess

    stop
end


