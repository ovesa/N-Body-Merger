module definitions ! this module simply lists all of my definitions
    integer :: i, j ,counter=0, k
    integer, PARAMETER :: N=40000
    real*8, dimension(3,N) :: x, v, grav
    real*8,  dimension(6,N/2) :: datak
    !real, allocatable, dimension(:) :: m(:)
    integer :: G, iseed, seed, t_final,t_dynamical
    real*8 :: ra_val, rp_val, ecc2, mu, a, xx, yy, zz, summa, b, randf, m
    real*8 :: period, t, time_step, dt,dump_time, is, rad, eps = 0.001 

    real*8 :: kinetic, potential,total_energy, sum_total
    real*8, PARAMETER :: pi =4.*atan(1.0)
    character*1000 :: filename, filename_energy, file_velocity
    character(100) :: outfile, output_energy, output_total, open_file
    CHARACTER(LEN=*), PARAMETER  :: FMT2 = '(a,i8.8,a)'






end module definitions



Program project
    use definitions
    use omp_lib

    Real*8   :: sigma, mean
    Real*4 :: D
   
    
    
     
   !allocate(grav(3,N),v(3,N),x(3,N), datak(6,N/2))

    open_file = '/home/users/ovesa/ASTRO506/particles/data2/data_t00000872.txt'
    !data_t00000443.txt
    open(unit=99,file=open_file)
    ! Enter the number of particles 
  
    do i = 1, N/2
            read(99,*) (datak(j,i), j=1,6)
    enddo

    do i = 1,N/2
        x(1,i) = datak(1,i) 
        x(2,i)= datak(2,i)
        x(3,i) = datak(3,i)
        v(1,i) = datak(4,i) + 0.2
        v(2,i) = datak(5,i)  + 0.1
        v(3,i) = datak(6,i) + .02

        x(1,i+(N/2)) = datak(1,i)  + 150
        x(2,i+(N/2)) = datak(2,i)  + 20
        x(3,i+(N/2)) = datak(3,i)
        v(1,i+(N/2)) = datak(4,i) -  0.2
        v(2,i+(N/2)) = datak(5,i) - 0.1
        v(3,i+(N/2)) = datak(6,i) - 0.02

    enddo

    sigma = .5 
    mean  = 0.
    G = 1.
    m = 1./N
   

    t = 0.00
    t_dynamical =  500
    !sqrt((3.*pi)/(16.*G*1.))
    t_final = t_dynamical ! amount to get one complete orbit using time step of 0.001
    !time_step = 1000.
    !dt = (t_final/10.)/(time_step)
    dt =  .1 !.1
    
    
    is = INT((t/dt) + 0.5) ! used to not write to file every timestep
    dump_time = 100



    
    rad = 0.5

  
    !iseed = 11
    !iseed = 11ç∂ç    seed = 1123

    !call init_random_seed2
    !call srand(seed)
    ! outputs
    !filename = '/Users/oanavesa/Desktop/GradSchool/FirstYear/ASTRO610/Project/positions_test.txt'
    !file_velocity = '/Users/oanavesa/Desktop/GradSchool/FirstYear/ASTRO610/Project/velocity_test.txt'
    !filename_energy = '/Users/oanavesa/Desktop/GradSchool/FirstYear/ASTRO610/Project/energy_test.txt'
    ! opening up the properfiles
    !OPEN(unit = 89, file = filename, status='replace')
    !OPEN(unit = 90, file = filename_energy, status='replace')
    !OPEN(unit = 91, file = file_velocity, status='replace')

    !write (*,*) grav


    !call init_random_seed
 
   call acceleration
   call energy
   total_energy = kinetic+potential 

   
   ! output for initial positions
   write(*,*) t
   outfile= '/home/users/ovesa/ASTRO506/meger_data_files/data_t0000000.txt'
   open(unit=92,file=outfile) 
   write(92,*) "# t = ",t
   write(92,*) "x   y   z   vx   vy   vz"
   do i=1,N
       write(92,'(3f12.5,3x,3es12.4,3x,es12.4)') x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
   enddo
   close(92)

     ! output for imitial energies
   output_energy='/home/users/ovesa/ASTRO506/meger_data_files/energy.txt'

   open(unit=93,file=output_energy)
   write(93,*) " t            KE            PE            E            KE/PE"
   write(93,'(5f12.6)') t,kinetic,potential,total_energy,kinetic/potential

   
   
  
    do
        !call srand(iseed)
        !do i = 1,N
        !10 x(:,N) = (/2.*rad*(rand(iseed)-0.5),2.*rad*(rand(iseed)-0.5),0.0/)
        !if (x(1,N)**2+x(2,N)**2+x(3,N)**2 > rad) goto 10

     !  call initial_positions
        !call write_to_file

        ! KDK leapfrog
        !write(*,*) x
        v = v + grav*(dt/2)
        x = x + v*(dt)
        t =t +dt
        call acceleration
        v = v + grav*(dt/2)

        
        
        if (mod(int(t/dt),100) .eq. 0) then
            write(*,*) t
            call energy
            total_energy = kinetic+potential 
            write(93,'(5f12.6)') t,kinetic,potential,total_energy,kinetic/potential

            !write(94,'(3f12.5,3x,3es12.4,3x,es12.4)') x(1,:),x(2,:),x(3,:),v(1,:),v(2,:),v(3,:)
            !write(*,*) t, x(:,1)

            write(outfile,'(a,I8.8,a)') '/home/users/ovesa/ASTRO506/meger_data_files/data_t',int(t),'.txt'
            !write(outfile,'(a,I8.8,a)') '/Users/oanavesa/Desktop/GradSchool/FirstYear/ASTRO610/Project/Data3/data_t',int(t),'.txt'
            open(unit=92,file=TRIM(outfile))
            write(92,*) "# t = ", t
            write(92,*) "x   y   z   vx   vy   vz"
            do i=1,N
                write(92,'(3f12.5,3x,3es12.4,3x,es12.4)') x(1,i),x(2,i),x(3,i),v(1,i),v(2,i),v(3,i)
            enddo
            close(92)
       endif
        
        !if (mod(int(t/dt),10) .eq. 0) write(*,*) t, x(:,1)
        !if (mod(counter,100) .eq. 0) write(90,*) t, potential, kinetic, total_energy
        !if (mod(counter,100) .eq. 0) write(91,*) t, v

       


        !write(*,*) t, x, potential
        !write(89,*) t, x
        !write(91,*) t, v
        !write(90,*) t, potential, kinetic, total_energy
        !end do
        if (t> t_final) exit
        end do ! t> 1000* with counter 1000


   
    close(unit=93)
    close(unit=94)    
    close(99)



end program project





subroutine acceleration
    use definitions  
    grav = 0.0

CALL OMP_SET_NUM_THREADS(40)
    !$OMP Parallel Do Default(Shared) PRIVATE(j,xx,yy,zz, summa)
     Do i=1,N 
        Do j=1,N
            xx = (x(1,i) - x(1,j))**2       
            yy = (x(2,i) - x(2,j))**2
            zz = (x(3,i) - x(3,j))**2
            summa = xx+yy+zz
            grav(:,i) = grav(:,i) - m*(x(:,i)-x(:,j))/(sqrt(summa + eps**2)**3)
        end do
     end do
    !$OMP END PARALLEL DO

     
 end subroutine acceleration








subroutine energy
    use definitions

kinetic = 0
potential = 0 

! kinetic energy = -1/2*sum(mi*v^2)
do i=1,N 
    kinetic = kinetic + (0.5)*(m*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2))
end do
! potential energu = -1/2* sum(Gmm/(ri-rj))

do i = 1,N
    do j = i,N
            xx = (x(1,i) - x(1,j))**2       
            yy = (x(2,i) - x(2,j))**2
            zz = (x(3,i) - x(3,j))**2
            summa = xx+yy+zz
            potential = potential - m*m*(G)/(sqrt(summa + eps**2))
        end do


end do

end subroutine energy


