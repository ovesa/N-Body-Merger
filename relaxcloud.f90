module definitions ! this module simply lists all of my definitions
    real*8, allocatable, dimension(:,:) :: x(:,:), v(:,:), grav(:,:)

    integer :: i, j, k
    integer :: N, G, seed, t_final, t_dynamical

    real*8 :: xx,yy,zz,summa, m, rad,xx2,yy2,zz2,summa2
    real*8 ::  t, time_step, dt, dump_time, eps = 0.001,eps2 = 0.001

    real*8 :: kinetic, potential,total_energy, sum_total
    real*8, PARAMETER :: pi =4.*atan(1.0)
    character*1000 :: filename, filename_energy, file_velocity
    character*100 :: outfile, output_energy, output_total, outfile2



end module definitions


Program project
    use definitions
    use omp_lib

    Real*8   :: sigma, mean
   
  
    
    ! Enter the number of particles 
    write(*,*) "Enter N:"
    read(*,*) N

    allocate(grav(3,N),v(3,N),x(3,N))

    sigma = .5 
    mean  = 0.
    G = 1. !gravitational constant. Set to 1 to be unitless
    m = 1./N !mass
   

    t = 0.00
    t_dynamical =  1 !sqrt((3.*pi)/(16.*G*1.))
    t_final = t_dynamical*1
    dt =  0.001
   
    


    call init_random_seed !calls a random seed number. Used to give random initial positions & velocities to particles
    rad = 0.5 !radius of the cloud
  


    ! gives initial positions & velocities using a Guassian Distribution
    do i=1,N
        80 x(1,i)= sigma*sqrt(-2*log(rand(seed)))*cos(2*pi*rand(seed)) + mean
        x(2,i) = sigma*sqrt(-2*log(rand(seed)))*cos(2*pi*rand(seed)) + mean
        x(3,i) = sigma*sqrt(-2*log(rand(seed)))*cos(2*pi*rand(seed)) + mean

        v(1,i) = sigma*sqrt(-2*log(rand(seed)))*cos(2*pi*rand(seed)) + mean
        v(2,i) = sigma*sqrt(-2*log(rand(seed)))*cos(2*pi*rand(seed)) + mean
        v(3,i) = sigma*sqrt(-2*log(rand(seed)))*cos(2*pi*rand(seed)) + mean

        sum_total = abs((x(1,i)**2+x(2,i)**2+x(3,i)**2))
        if (sum_total > rad) goto 80
       
    end do


    
    call energy
    total_energy = kinetic+potential 

  
   ! output for initial positions
    outfile= '/home/users/ovesa/ASTRO506/particles/try/data_t00000000.txt'
    open(unit=92,file=outfile) 
    write(92,*) "# t = ",t
    write(92,*) "x   y   z   vx   vy   vz"
    do i=1,N
       write(92,'(3f12.5,3x,3es12.4,3x,es12.4)') x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
    enddo
    close(92)

    ! output for initial energies
    output_energy='/home/users/ovesa/ASTRO506/particles/try/energy.txt'
    open(unit=93,file=output_energy)
    write(93,*) " t            KE            PE            E            KE/PE"
    write(93,'(5f12.6)') t,kinetic,potential,total_energy,kinetic/potential


  

    do
        ! KDK leapfrog
        v = v + grav*(dt/2)
        x = x + v*(dt)
        t =t +dt
        call acceleration
        v = v + grav*(dt/2)

        
        if (mod(int(t/dt),100) .eq. 0) then
            write(*,*)  t
            call energy
            total_energy = kinetic+potential 
            write(93,'(5f12.6)') t,kinetic,potential,total_energy,kinetic/potential

            write(outfile,'(a,I8.8,a)') '/home/users/ovesa/ASTRO506/particles/try/data_t',int(t),'.txt'
            open(unit=92,file=TRIM(outfile))
            write(92,*) "# t = ", t
            write(92,*) "x   y   z   vx   vy   vz"
            do i=1,N
                write(92,'(3f12.5,3x,3es12.4,3x,es12.4)') x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
            enddo
            close(92)
        endif
        
      
     !  if (t .eq. t_final) then
     !   write(97,*) "# t = ",t
     !   write(97,*) "x   y   z   vx   vy   vz"
     !   write(97,'(3f12.5,3x,3es12.4,3x,es12.4)') x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
     !  endif


       
        if (t> t_final) exit
        end do 


   
    close(unit=93)
    close(unit=94)


end program project





subroutine acceleration
    use definitions  
    grav = 0.0

    call omp_set_num_threads(25)
     !$OMP Parallel Do Default(Shared)  PRIVATE(j,xx,yy,zz,summa)
     Do i=1,N 
        Do j=1,N
            xx = (x(1,i) - x(1,j))**2       
            yy = (x(2,i) - x(2,j))**2
            zz = (x(3,i) - x(3,j))**2
            summa = xx+yy+zz
            grav(:,i) = grav(:,i) - m*(x(:,i)-x(:,j))/(sqrt(summa+ eps**2)**3)
        end do
     end do
    !$OMP END PARALLEL DO

     
 end subroutine acceleration





subroutine energy
    use definitions

kinetic = 0.0
potential = 0.0

! kinetic energy = -1/2*sum(mi*v^2)
do i=1,N 
    kinetic = kinetic + (0.5)*(m*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2))
end do
! potential energu = -1/2* sum(Gmm/(ri-rj))

do i = 1,N
    do j = i,N
            xx= (x(1,i) - x(1,j))**2       
            yy = (x(2,i) - x(2,j))**2
            zz = (x(3,i) - x(3,j))**2
            summa = xx+yy+zz
            potential = potential - m*m*(G)/(sqrt(summa+ eps**2))
        end do
end do


end subroutine energy

subroutine init_random_seed

    INTEGER :: i, k, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = k)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, k) /)
    CALL RANDOM_SEED(PUT = seed)

    !DEALLOCATE(seed)
end
