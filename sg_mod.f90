module signal_generator
    use ifport
    
    type struct
        integer*4 rxCount(2)            ! number of recievers (X, Y)
        integer*4 txCount(3)            ! number of transmitters (X, Y, Z)
        integer*4 windowHeight          ! number of sample
        integer*4 windowCount
        integer*4 txCoord(3)            ! coordinates of real transmitter (X, Y, Z) (km)
        real*4 X(2)                     ! dimensions of cube (km)
        real*4 Y(2)
        real*4 Z(2)        
        real*4 dT                       ! time of discretization (s)
        real*4 alpha
        real*4 omega
        real*4 noiseCoef
        real*4 startTime                ! time of micro-explosion (s)
        real*4 speed                    ! speed of wave (km/s)
    end type struct  
!------------------------------------------------------------------------------
!-----CONTAINS-----------------------------------------------------------------
!------------------------------------------------------------------------------     
    contains    
!------------------------------------------------------------------------------
!-----GENERATE DATA------------------------------------------------------------
!------------------------------------------------------------------------------
    subroutine SG_generateData(inpData, output)
		!----------------------------------------------------------------------
        type (struct) inpData        
        real*4 output(inpData%windowHeight, inpData%rxCount(1) * inpData%rxCount(2))
        !----------------------------------------------------------------------          
        integer*4 timeArray(3)
        real*4 rand_number
	    real*4 distance(inpData%rxCount(1), inpData%rxCount(2))              ! distance between supposing transmitter and every reciever        
        real*4 time
        real*4 amp
        real*4 step(2)                
        step(1) = (inpData%X(2) - inpData%X(1)) / (inpData%rxCount(1) - 1)   ! distance between recievers by X and Y axes
        step(2) = (inpData%Y(2) - inpData%Y(1)) / (inpData%rxCount(2) - 1)
		!----------------------------------------------------------------------- 
        call itime(timeArray)                                                ! get the current time
        j = rand(timeArray(1)+timeArray(2)+timeArray(3))
        do j = 1, inpData%rxCount(2)            
            do i = 1, inpData%rxCount(1)
                distance(i, j) = sqrt((inpData%X(1) + (i-1) * step(1) - inpData%txCoord(1))**2 + (inpData%Y(1) + (j-1) * step(2) - inpData%txCoord(2))**2 + (inpData%txCoord(3))**2)
                time = inpData%startTime + distance(i, j) / inpData%speed
                amp = inpData%txCoord(3) / ((distance(i, j))**2)                               
                do k = 1, inpData%windowHeight
                    output(k, i + (j-1) * inpData%rxCount(1)) = amp * exp(-inpData%alpha * inpData%omega * (inpData%dT * (k-1) - time)**2) * cos(4.0 * acos(0.0) * inpData%omega * (inpData%dT * (k-1) - time))                   
                    rand_number = rand(0)                    
                    output(k, i + (j-1) * inpData%rxCount(1)) = output(k, i + (j-1) * inpData%rxCount(1)) +((rand_number * 2) - 1) * inpData%noiseCoef * amp ! -1 <= noise < 1    
					!output(k, i + (j-1) * inpData%rxCount(1)) = ((rand_number * 2) - 1) * inpData%noiseCoef * amp ! Only  noise                   
                end do
            end do
        end do        
    end subroutine SG_generateData
end module signal_generator
