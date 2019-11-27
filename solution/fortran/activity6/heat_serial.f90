
#define N           27
#define MAX_ITER    4000
#define TOL         1e-4
#define MAX_TEMP    100.0

module grid
    real, dimension(:,:), allocatable ::T_old,T_new
end module grid

program heat_serial
    use grid

    implicit none
    integer :: rows,cols
    integer :: i,j,iter
    real :: dT
    rows = N+2
    cols = N+2
    ! allocate memory
    allocate(T_old(rows,cols))
    allocate(T_new(rows,cols))

    call initialize(rows,cols)
    iter = 0
    dT = MAX_TEMP
    do while ( (dT > TOL) .AND. (iter <= MAX_ITER) )
        do i = 2,rows-1
            do j = 2,cols-1
               T_new(i,j) = 0.25 * ( T_old(i-1,j) + T_old(i+1,j) + T_old(i,j-1) + T_old(i,j+1) )
            end do
        end do

        dT = 0.0;
        do i = 2,rows-1
            do j = 2,cols-1
               dT = max(abs(T_new(i,j) - T_old(i,j)) , dT)
               T_old(i,j) = T_new(i,j)
            end do
        end do
        iter=iter+1
    end do

    !writing the output to the file is optional and is there for validation.
    call write_grid(rows,cols,iter)

    if ((iter - 1) == MAX_ITER) then
        print *, "Reached maximum iterations ",iter,". Error = ", dT
    else
        print *, "Converged in ",iter," iterations with and error of ",dT
    end if

end program heat_serial


! Initialize the grids with appropirate boundary and initial conditions
subroutine initialize(rows,cols)
    use grid

    implicit none
    integer :: rows,cols
    integer :: i,j

    do i=1,rows
        do j=1,cols
                T_old(i,j) = 0.75 * MAX_TEMP
                T_new(i,j) = 0.75 * MAX_TEMP
        end do
    end do

    do j=1,cols
        T_old(rows,j) = MAX_TEMP
        T_new(rows,j) = MAX_TEMP
    end do
    do i=1,rows
        T_old(i,1) = MAX_TEMP
        T_new(i,1) = MAX_TEMP
    end do
end subroutine initialize

! This is optional. Write the 2D grid in a file
subroutine write_grid(rows,cols,iter)
    use grid
    implicit none
    integer :: rows,cols,iter
    integer :: i,j
    character(len=32) :: fname
    write(fname,'(A,I0,A)') 'output_serial_t',iter,'.txt'

    open (unit = 10, file = fname, status='replace')
    do i=2,rows-1
       write(10,*) (T_old(i,j), j=2,cols-1)
    end do
    close(10)
end subroutine write_grid
