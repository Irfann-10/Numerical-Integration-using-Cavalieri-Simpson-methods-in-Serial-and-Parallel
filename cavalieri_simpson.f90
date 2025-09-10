program cavalieri_simpson_parallel_serial
  use mpi
  implicit none
  
  integer :: rank, size, ierr
  integer :: n, i, local_n, local_start, local_end
  real :: a, b, h, x, local_integral, parallel_integral, serial_integral
  double precision :: start_time, end_time, parallel_time, serial_time
  
  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  
  ! Number of intervals
  n = 1000000
  
  ! Integration limits
  a = 1.0
  b = 5.0
  
  ! Calculate local number of intervals for each process
  local_n = n / size
  
  ! Calculate local start and end points for each process
  local_start = rank * local_n
  local_end = (rank + 1) * local_n - 1
  
  ! Calculate step size
  h = (b - a) / real(n)
  
  ! Initialize local integral value
  local_integral = 0.0
  
  ! Start the timer for parallel version
  start_time = MPI_Wtime()
  
  ! Perform local integration
  do i = local_start, local_end
    x = a + i * h
    
    ! Calculate the function value at x
    local_integral = local_integral + f(x) + 4.0 * f(x + h/2.0) + f(x + h)
  end do
  
  ! Multiply by h/6 to complete the local integration
  local_integral = h/6.0 * local_integral
  
  ! Reduce local integrals to obtain the global integral
  call MPI_Reduce(local_integral, parallel_integral, 1, MPI_REAL, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
  ! End the timer for parallel version
  end_time = MPI_Wtime()
  
  ! Calculate the parallel execution time
  parallel_time = end_time - start_time
  
  ! Print parallel result and execution time from rank 0
  if (rank == 0) then
    write(*,*) "Parallel Approximate integral: ", parallel_integral
    write(*,*) "Parallel Execution Time: ", parallel_time
  end if

 ! Perform the integration in serial on process 0
  if (rank == 0) then
    ! Start the timer for serial version
    start_time = MPI_Wtime()
    
    ! Perform the integration in serial
    do i = 0, n-1
      x = a + i * h
    
      ! Calculate the function value at x
      serial_integral = serial_integral + f(x) + 4.0 * f(x + h/2.0) + f(x + h)
    end do
  
    ! Multiply by h/6 to complete the integration
    serial_integral = h/6.0 * serial_integral
    
    ! End the timer for serial version
    end_time = MPI_Wtime()
    
    ! Calculate the serial execution time
    serial_time = end_time - start_time
    
    ! Print serial result and execution time on process 0
    write(*,*) "Serial Approximate integral: ", serial_integral
    write(*,*) "Serial Execution Time: ", serial_time
  end if
  
  ! Finalize MPI
  call MPI_Finalize(ierr)
  
contains

  ! Function to be integrated (modify it according to your needs)
  real function f(x)
    real, intent(in) :: x
    f = (6*x**4)-(7*x**3)+4*x  ! Example function: 6x^4-7x^3+4x
  end function f
  
end program cavalieri_simpson_parallel_serial
