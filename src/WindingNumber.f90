program WindingNumber
  implicit none
  integer          :: n, ios, i, m, j
  real             :: ip, sa, sb, wn, ov
  real, allocatable:: lon_k(:), lat_k(:), lon_t(:), lat_t(:), theta(:)
  real, parameter  :: pi=3.14159265

  ! Input closed polygon
  open(10, file='data/Polygon.txt', status='old', action='read')
  n = 0
  do
     read(10, *, iostat=ios)
     if (ios == 0) then
        n = n + 1
     else if (ios < 0) then
        exit
     end if
  end do
  rewind(10)
  allocate(lon_k(1:n), lat_k(1:n))
  do i = 1, n
     read(10, *) lon_k(i), lat_k(i)
  end do
  close(10)
  
  ! Input target points you wish to judge if these wre inside or outside of polygon
  open(20, file='data/TargetPoints.txt', status='old', action='read')
  m = 0
  do
     read(20, *, iostat=ios)
     if (ios == 0) then
        m = m + 1
     else if (ios < 0) then
        exit
     end if
  end do
  rewind(20)
  allocate(lon_t(1:m), lat_t(1:m))
  do j = 1, m
     read(20, *) lat_t(j), lon_t(j)
  end do
  close(20)
  
  ! Target point is inside polygon if wn = 1.000 
  open(30, file='work/WindingNumber.txt', status='replace', action='write')
  do j = 1, m
     allocate(theta(1:n-1))
     do i = 1, n-1
        ip = (lon_t(j)-lon_k(i))*(lon_t(j)-lon_k(i+1))+(lat_t(j)-lat_k(i))*(lat_t(j)-lat_k(i+1))
        sa = sqrt((lon_t(j)-lon_k(i))  **2+(lat_t(j)-lat_k(i))  **2)
        sb = sqrt((lon_t(j)-lon_k(i+1))**2+(lat_t(j)-lat_k(i+1))**2)
        theta(i) = acos(ip/(sa*sb))

        ! Judge if spanning angle between two pointis is clockwise or
        ! counter-clockwise by sign of outer vector
        ! If sign of outer vector is negative, then flip sign of spanning angle
        ov = (lon_t(j)-lon_k(i))*(lat_t(j)-lat_k(i+1))-(lon_t(j)-lon_k(i+1))*(lat_t(j)-lat_k(i))
        if (ov < 0) then
           theta(i) = theta(i) * (-1.0)
        end if
        
        if (isnan(theta(i))) then
           theta(i) = 0.000
        end if

     end do
     wn = (1.0/(2.0*pi))*sum(theta(1:n-1))
     write(30, '(3f10.3)') lon_t(j), lat_t(j), wn
     deallocate(theta)
  end do
  close(30)  
end program WindingNumber
