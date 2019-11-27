
! 2D decomposition of the domain a 2D grid representing a
! metal plate.
! N :          Number of points such that NxN square grid will be created as a domain
! MAX_ITER:    Maximum number of iterations to run
! TOL:         Minimum delta T such that below this value no change is deemed and the simulation needs to stop
! MAX_TEMP:    Temperature on the boundary

#define N           27
#define MAX_ITER    4000
#define TOL         1e-4
#define MAX_TEMP    100.0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating a module for globalizing the allocatable arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module grid
    real (kind=4), dimension(:,:), allocatable ::T_old,T_new
end module grid


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Creating a module for globalizing mpi related data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module my_mpi_vars
    use mpi

  ! variables to store virtual topology information
    ! Global rank and number of MPI processes
    integer ::  grank, gprocs
    integer ::  globalN, localN
    ! Local rank and number of MPI processes
    integer lrank, lprocs
    integer :: ndims = 2                                    ! Since we are decomposing grid in 2D
    integer, dimension(2) :: dims                           ! placeholder for dimensions
    integer, dimension(2) :: mycoords                       ! placeholder for Cartesian coordinates
    integer, dimension(2) :: period = (/0,0/)               ! placeholder for periodicity information
    integer :: reorder = 0

    integer :: x_axis = 1, y_axis = 0                       ! convinience variables for referencing axes
    integer :: right, left, up, down                        ! placeholder for ranks of nearest neighbours

    integer :: cartcomm                                     ! Declaration of new communicator
    integer :: ierr                                         ! error handler for MPI calls

    integer, dimension(8)  :: req                           ! placeholder for MPI Isend and Irecv requests
    integer, dimension(MPI_STATUS_SIZE,8)  :: status        ! placeholder for maintaining status of the MPI requests. Size defined in mpif.h
end module my_mpi_vars



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program heat_mpi_2d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use grid
    use mpi
    use my_mpi_vars

    implicit none


    ! Local rows and columns along with ghost region
    integer rows, cols;


    integer :: i,j,iter,rem
    integer :: flag
    real (kind=4) :: dT,dT_local


    ! Initialize MPI
    call MPI_Init(ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD,grank,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,gprocs,ierr)

    call is_grid_decomposible(flag)
    if ( flag .ne. 0) then
        if (grank .eq. 0) then
            write(0,'(A)') "Grid decomposition will fail."
            write(0,'(A,I0,A)') "nprocs = ",gprocs,"     --- Number of MPI Processes"
            write(0,'(A,I0,A)') "N      = ",N,"          --- Number of grid points in each dimension"
            write(0,'(A)') "nprocs should be a perfect square (e.g. 1,4,9,16...)"
            write(0,'(A)')"Also, N should be exactly divisible by sqrt(nprocs)"
            flush(0)
            call MPI_Abort(MPI_COMM_WORLD,911,ierr);
        end if
    end if
 ! It may be a good idea to create a Cartesian Communicator
 ! Here are some hints:
 ! 1.   Get a good guess on dimensions
 ! 2.   Create a Cartesian communicator and call it cartcomm
 ! 3.   Query rank and its coordinates in new communicator
 ! 4.   Query who are the neighbours
 !


 ! Domain decomposition: decide work for each rank by figuring out localN
 ! hint : globalN/number of MPI ranks


    rows = localN + 2
    cols = localN + 2

    ! allocate memory
    allocate(T_old(rows,cols))
    allocate(T_new(rows,cols))

    call initialize(rows,cols)
    iter = 0
    dT = MAX_TEMP

    do while ( (dT > TOL) .AND. (iter <= MAX_ITER) )

        !   update halo region
        call halo_update(rows,cols)

        !   Evaluate temperature in inner domain
        do i = 2,rows-1
            do j = 2,cols-1
               T_new(i,j) = 0.25 * ( T_old(i-1,j) + T_old(i+1,j) + T_old(i,j-1) + T_old(i,j+1) )
            end do
        end do

        ! Figure out the maximum local error
        dT = 0.0
        dT_local = 0.0
        do i = 2,rows-1
            do j = 2,cols-1
               dT_local = max(abs(T_new(i,j) - T_old(i,j)) , dT_local)
               T_old(i,j) = T_new(i,j)
            end do
        end do

       ! Communicate to everyone the maximum of local errors from all ranks
       ! hint: use MPI_Allreduce with  MPI_MAX operation

        iter=iter+1
    end do

    !writing the output to the file is optional and is there for validation.
    call write_grid(rows,cols,iter)
    call MPI_Barrier(cartcomm,ierr)
    !call print_grid(rows,cols,iter)

    if (lrank == 0) then
        if ((iter - 1) == MAX_ITER) then
            print *, "Reached maximum iterations ",iter,". Error = ", dT
        else
            print *, lrank,": Converged in ",iter," iterations with and error of ",dT
        end if
    end if

call MPI_FINALIZE(ierr)
end program heat_mpi_2d


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize the grids with appropirate boundary and initial conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize(rows,cols)
    use grid
    use my_mpi_vars
    implicit none
    integer :: rows,cols
    integer :: i,j

    do i=1,rows
        do j=1,cols
                T_old(i,j) = 0.75 * MAX_TEMP
                T_new(i,j) = 0.75 * MAX_TEMP
        end do
    end do

    if (mycoords(2) == 0) then
        do i=1,rows
            T_old(i,1) = MAX_TEMP
            T_new(i,1) = MAX_TEMP
        end do
    end if
    if (mycoords(1) == dims(1)-1) then
        do j=1,cols
            T_old(rows,j) = MAX_TEMP
            T_new(rows,j) = MAX_TEMP
        end do
    end if
end subroutine initialize



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine for Halo exchange
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine halo_update(rows,cols)
    use grid
    use my_mpi_vars
    implicit none
    integer :: rows,cols

     ! Send and receive ghost regions to the 4 neighbors
     ! Hint:
     ! 1. Create temporary arrays for send and receive buffers for each side
     ! 2. load buffers
     ! 3. communicate synchronous or asynchronous, its your preference
     ! 4. unload data in receive buffer to the respective rows and columns of 2D array T.

 end subroutine halo_update



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! write grid
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_grid(rows,cols,iter)
    use grid
    use my_mpi_vars
    integer :: i,j,iter,rows,cols,root
    integer :: x,y,indx                            ! some scratch variables
    integer :: numelems_local,recv_buf_size
    ! following allocations are significant on to root rank
    integer :: recv_nelems(lprocs)                 ! number of elements in a vector to receive from each rank
    integer :: disp(lprocs)                        ! array to hold offsets where data from each rank will go in the receive buffer

    real (kind=4), dimension(:), allocatable :: sendbuf,recvbuf      ! send buffer to send the data to root and recvbuffer is only signifcant to root
    character(len=256)   ::  fname                         ! output filename

    root=0
    numelems_local = localN * localN
    ! tell Writer how many elements each rank is going to send
    call MPI_Gather(numelems_local, 1, MPI_INT, recv_nelems, 1, MPI_INT, root,cartcomm,ierr)

    ! disp is an array of offsets or displacements
    ! i.e. start index where to receive vector from each rank (useful only on master)
    if (lrank == root) then
        disp(1) = 0
        do i = 2,lprocs
            disp(i) = disp(i-1)+ recv_nelems(i)
        end do
    end if

    ! Now prepare for gathering the actual data
    allocate(sendbuf(numelems_local))
    indx=1
    do i = 2, rows-1
       do j = 2, cols-1
            sendbuf(indx) = T_old(i,j)
            indx=indx+1
       end do
    end do
    recv_buf_size=0
    do i = 1,lprocs
        recv_buf_size = recv_buf_size + recv_nelems(i)
    end do

    if (lrank == root) then
        allocate(recvbuf(recv_buf_size))
    end if

    call MPI_Gatherv(sendbuf, numelems_local, MPI_FLOAT, recvbuf, recv_nelems,disp, MPI_FLOAT, root, cartcomm,ierr);

    if (lrank == root) then
        write(fname,'(A,I0,A)') 'output_t',iter,'.txt'
        open (unit = 10, file = fname, status='replace')
        indx=0
        do y=1,recv_buf_size,dims(2)*numelems_local
            do i=1,numelems_local,localN
                do x=1,dims(1)*numelems_local,numelems_local
                    do j=1,localN
                        indx = (y-1) + (i-1) + (x-1) + j
                        write(10,fmt='((F10.5))',advance='no') recvbuf(indx)
                    end do
                end do
                write(10,*)
            end do
        end do
        close(10)
    end if
end subroutine write_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! general validity check
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine is_grid_decomposible(flag)

    use my_mpi_vars
    implicit none
    ! condition 1 : N should be exactly divisible by sqrt(nprocs)
    ! condition 2 : nprocs should be a perfect square i.e. sqrt(nprocs) is a whole number
    integer :: cond_1=0,cond_2=0
    integer :: proc_1D,flag

    proc_1D = sqrt(real(gprocs))
    ! Check first condition
    if ( mod(N,proc_1D) .ne. 0) then
            cond_1 = 1
    end if
    ! Check second condition
    if ( (gprocs/proc_1D) .ne. proc_1D ) then
            cond_2 = 1
    end if
    if ( (cond_1==0)  .and. (cond_2==0) ) then
        flag=0
    else
        flag=1
    end if

end subroutine is_grid_decomposible
