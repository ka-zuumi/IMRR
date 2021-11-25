module velocityScaling
use PARAMETERS
use ANALYSIS
implicit none

contains




!subroutine scaleVelocities(coordsIN,&
!                           momentaIN)
subroutine scaleVelocities(coordsIN,&
                           momentaIN,&
                           momentaOther)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
real(dp),dimension(3*Natoms),intent(in) :: momentaOther
integer,dimension(1,Natoms) :: clusters1
integer,dimension(2,Natoms) :: clusters2

logical :: successfulScaling

integer :: rxnProgress

integer :: i

!if (noVelocityScaling_flag) return

select case(velocityScalingID)

! ID = 0
! No scaling
case(0)
    return

! ID = 1
! This was the OLD scheme
case(1)

    call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                      clusters1,clusters2)
    
    call updateRxnProgress1(coordsIN,clusters2,rxnProgress)
    
    ! This is algorithm A
    if (rxnProgress == 1) then
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm B (analytical derivative)
    else if (.true.) then
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm B (numerical derivative)
    else
        call scaleVelocities3(coordsIN,momentaIN,clusters2)
    end if

! ID = 2
! This is the NEW (and more robust!) scheme
case(2)

    call updateRxnProgress2(coordsIN,momentaIN)
    
    ! This is algorithm B
    if (rxnStage == 0) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm A
    else if (rxnStage < 3) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm B
    else
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    end if

! ID = 3
! This is the NEWER (and even more robust!) scheme
case(3)

    call updateRxnProgress2(coordsIN,momentaIN)
    
    ! This is algorithm B*
    if (rxnStage == 0) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2plus(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm A*
    else if (rxnStage < 3) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1plus(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm B*
    else
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2plus(coordsIN,momentaIN,clusters2)
    end if

! ID = 4
! This is the NEWER (and even more robust!) scheme
case(4)

    call updateRxnProgress2(coordsIN,momentaIN)

    if (Naccept == 0) then
        resetPCMcounter = resetPCMcounter + 1
        if (resetPCMcounter == resetPCMcounter_max) then
            resetPCMcounter = 0
            call resetPCM(momentaIN)
            print *, ""
            print *, " Resetting the PCM to zero!"
            print *, ""
        end if
    end if
    
    ! This is algorithm *B
    if (rxnStage == 0) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm *A
    else if (rxnStage < 3) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm *B
    else
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    end if

! ID = 5
! This is the new temporary scheme to test if
! algorithm A works
case(5)
    
    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities1(coordsIN,momentaIN,clusters2)

! ID = 6
! This is the new temporary scheme to test if
! algorithm *A works
case(6)

    if (Naccept == 0) then
        resetPCMcounter = resetPCMcounter + 1
        if (resetPCMcounter == resetPCMcounter_max) then
            resetPCMcounter = 0
            call resetPCM(momentaIN)
            print *, ""
            print *, " Resetting the PCM to zero!"
            print *, ""
        end if
    end if
 
    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities1(coordsIN,momentaIN,clusters2)

! ID = 7
! AKA the BAB-glue scheme
case(7)

    if (Naccept == 0) then
        call updateRxnProgress3(coordsIN,momentaIN)
        resetPCMcounter = resetPCMcounter + 1
        if (resetPCMcounter == resetPCMcounter_max) then
            resetPCMcounter = 0
            call resetPCM(momentaIN)
            print *, ""
            print *, " Resetting the PCM to zero!"
            print *, ""
        end if
    end if
    
    ! This is algorithm *B
    if ((rxnStage == 0).or.(rxnStage==4)) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm *A
    else if (rxnStage == 2) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    end if

! ID = 8
! AKA the adjusted BAB-glue scheme
! (it uses a different updateRxnProgress)
case(8)

    if (Naccept == 0) then
        call updateRxnProgress4(coordsIN,momentaIN)
        resetPCMcounter = resetPCMcounter + 1
        if (resetPCMcounter == resetPCMcounter_max) then
            resetPCMcounter = 0
            call resetPCM(momentaIN)
            print *, ""
            print *, " Resetting the PCM to zero!"
            print *, ""
        end if
    end if
    
    ! This is algorithm *B
    if ((rxnStage == 0).or.(rxnStage==4)) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm *A
    else if (rxnStage == 2) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    end if

! ID = 9
! AKA the BAB-glue scheme
! but now only scaling after interpolation
case(9)

    call updateRxnProgress3(coordsIN,momentaIN)
    resetPCMcounter = resetPCMcounter + 1
    if (resetPCMcounter == resetPCMcounter_max) then
        resetPCMcounter = 0
        call resetPCM(momentaIN)
        print *, ""
        print *, " Resetting the PCM to zero!"
        print *, ""
    end if

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is algorithm *B
    if ((rxnStage == 0).or.(rxnStage==4)) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm *A
    else if (rxnStage == 2) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    end if

! ID = 10
! AKA the adjusted BAB-glue scheme
! (it uses a different updateRxnProgress)
! but now only scaling after interpolation
case(10)

    call updateRxnProgress4(coordsIN,momentaIN)
    resetPCMcounter = resetPCMcounter + 1
    if (resetPCMcounter == resetPCMcounter_max) then
        resetPCMcounter = 0
        call resetPCM(momentaIN)
        print *, ""
        print *, " Resetting the PCM to zero!"
        print *, ""
    end if

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is algorithm *B
    if ((rxnStage == 0).or.(rxnStage==4)) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm *A
    else if (rxnStage == 2) then
    
        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0
    
        call scaleVelocities1(coordsIN,momentaIN,clusters2)
    
    end if

! ID = 11
! The sleek new algorithm C which is for a
! single cluster and no PCM shifting
! (iterative algorithm C)
case(11)

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities4(coordsIN,momentaIN,clusters2)

! ID = 12
! The sleek new algorithm C which is for a
! single cluster and no PCM shifting
! (analytical algorithm C)
case(12)

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities5(coordsIN,momentaIN,clusters2)

! ID = 13
! The sleek new algorithm C which is for a
! single cluster and no PCM shifting
! (numerical-analytical algorithm C)
case(13)

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities6(coordsIN,momentaIN,clusters2)

! ID = 14
! The brute-force algorithm D which subtracts out
! the center of mass translational momentum and
! then scales to the total energy of the system
case(14)

!!  if (.not.(prevNaccept_was_nonzero)) return

    ! Does not require any sort of clustering
    ! (because it is one-cluster by default)
    
    call scaleVelocities7(coordsIN,momentaIN,clusters2)

! ID = 15
! A repeat of ID = 13 (algorithm C) but it
! rescales whenever possible
case(15)

    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities6(coordsIN,momentaIN,clusters2)

! ID = 16
! A repeat of ID = 13 (algorithm C) but it
! rescales whenever possible but starting
! from a slightly messed-up gradient
case(16)

    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    ! Mess up the momenta slightly
    do i = 1, 3*Natoms
        momentaIN(i) = momentaIN(i) * &
              (((modulo(i,2)-0.5)*0.4) + 1.0)
!             (((rand()-0.5)*0.4) + 1.0)
    end do

    call scaleVelocities6(coordsIN,momentaIN,clusters2,&
                          successfulScaling)

    if (.not.(successfulScaling)) then
    do i = 1, 3*Natoms
        momentaIN(i) = momentaIN(i) / &
              (((modulo(i,2)-0.5)*0.4) + 1.0)
    end do
    end if

! ID = 17
! The sleek new MODIFIED algorithm C which is for a
! single cluster and no PCM shifting
! (numerical-analytical algorithm modC)
case(17)

    ! This is my brute-force way of making
    ! it do a one-cluster analysis instead
    ! of changing the subroutine
    clusters2(1,:) = 1
    clusters2(2,:) = 0

    call scaleVelocities9(coordsIN,momentaIN,clusters2)

! ID = 20
case(20)

    ! Here, we must update the rxnStage OUTSIDE of
    ! the velocity scaling because we don't want
    ! to interfere with a rewind

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is algorithm C from ID = 13
    if (rxnStage == 0) then

        ! This is my brute-force way of making
        ! it do a one-cluster analysis instead
        ! of changing the subroutine
        clusters2(1,:) = 1
        clusters2(2,:) = 0

        call scaleVelocities6(coordsIN,momentaIN,clusters2)
    
    ! This is algorithm B
    elseif (rxnStage == 2) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    end if

! ID = 21
! This is pure algorithm D
case(21)

    call scaleVelocities8(coordsIN,momentaIN,momentaOther)

! ID = 22
! This is pure algorithm D but we mess up the momenta
! if there was an unsuccessful scaling so that VENUS
! will be forced to rewind
case(22)

    call scaleVelocities8(coordsIN,momentaIN,momentaOther,&
                          successfulScaling)

    if (.not.(successfulScaling)) then
    do i = 1, 3*Natoms
        momentaIN(i) = sqrt(2*C1*(2*H_baseline)/&
                            sum(1.0d0 / masses(1:NATOMS)))
    end do
    end if

! ID = 30
case(30)

    ! Here, we must update the rxnStage OUTSIDE of
    ! the velocity scaling because we don't want
    ! to interfere with a rewind

!!  if (.not.(prevNaccept_was_nonzero)) return
    
    ! This is algorithm D from ID = 21
    if (rxnStage == 0) then

        call scaleVelocities8(coordsIN,momentaIN,momentaOther)
    
    ! This is algorithm B
    elseif (rxnStage == 2) then
        call clusterAtoms(reshape(coordsIN, (/ 3, NATOMS /)),&
                          clusters1,clusters2)
        call scaleVelocities2(coordsIN,momentaIN,clusters2)
    
    end if


! ID = unknown
! No scaling
case default
    return
    
end select

return

end subroutine scaleVelocities

subroutine resetPCM(momentaIN)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(inout) :: momentaIN

real(dp),dimension(3) :: momentumTCM

integer :: i, startIndex, endIndex

momentumTCM = 0.0d0
startIndex = 1
endIndex = 3
do i = 1, Natoms
    momentumTCM = momentumTCM + &
        momentaIN(startIndex:endIndex)
    startIndex = startIndex + 3
    endIndex = endIndex + 3
end do

momentumTCM = momentumTCM / Natoms

startIndex = 1
endIndex = 3
do i = 1, Natoms
    momentaIN(startIndex:endIndex) = &
        momentaIN(startIndex:endIndex) -&
        momentumTCM
    startIndex = startIndex + 3
    endIndex = endIndex + 3
end do

return

end subroutine resetPCM

subroutine resetPCMweighted(momentaIN)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(inout) :: momentaIN

real(dp),dimension(3) :: velocityTCM

integer :: i, startIndex, endIndex

velocityTCM = 0.0d0
startIndex = 1
endIndex = 3
do i = 1, Natoms
    velocityTCM = velocityTCM + &
        momentaIN(startIndex:endIndex)
    startIndex = startIndex + 3
    endIndex = endIndex + 3
end do
velocityTCM = velocityTCM / sum(masses(1:Natoms))

startIndex = 1
endIndex = 3
do i = 1, Natoms
    momentaIN(startIndex:endIndex) = &
        momentaIN(startIndex:endIndex) -&
        masses(i)*velocityTCM
    startIndex = startIndex + 3
    endIndex = endIndex + 3
end do

return

end subroutine resetPCMweighted







subroutine updateRxnProgress1(coordsIN,clusters,rxnProgress)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
integer,dimension(2,Natoms),intent(in) :: clusters
integer,intent(out) :: rxnProgress

real(dp) :: distij
real(dp) :: distij_threshold = 25.0
integer :: i, j, i3, j3

! To monitor the reaction progress, just look
! at the shortest distance between any two
! atoms in the two clusters

! By default we assume the reaction is in
! stage zero
rxnProgress = 0
do i = 1, Natoms
    if (clusters(1,i) > 0) then
        i3 = 3*i
        do j = 1, Natoms
            if (clusters(2,j) > 0) then
                j3 = 3*j
                distij = ((coordsIN(i3)-coordsIN(j3))**2)+&
                         ((coordsIN(i3-1)-coordsIN(j3-1))**2)+&
                         ((coordsIN(i3-2)-coordsIN(j3-2))**2)

                ! If it meets our various criteria then
                ! it has reached stage one
                if (distij < distij_threshold) then
                    rxnProgress = 1
                    exit
                end if
            end if
        end do
    end if
end do

return

end subroutine updateRxnProgress1

subroutine updateRxnProgress2(coordsIN,momentaIN)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(in) :: &
            coordsIN, momentaIN

real(dp) :: distBrC
!real(dp) :: distBrC_threshold1 = 36.0
!real(dp) :: distBrC_threshold2 = 25.0

! To monitor the reaction progress, just look
! at the distance between the carbon and the
! bromine atoms

distBrC = ((coordsIN(4)-coordsIN(10))**2)+&
          ((coordsIN(5)-coordsIN(11))**2)+&
          ((coordsIN(6)-coordsIN(12))**2)

! In rxnStage 0, we are using algorithm B and
! check to see whether we should switch
if ((rxnStage == 0).and.(distBrC < rxnProgress2_threshold1)) then
    print *, ""
    print *, " Switching to reaction stage 1"
    print *, ""
    rxnStage = 1

! In rxnStage 1, we are using algorithm A and
! check to see whether we have reached the
! "middle" region of the intermediate region
else if ((rxnStage == 1).and.(distBrC < rxnProgress2_threshold2)) then
    print *, ""
    print *, " Switching to reaction stage 2"
    print *, ""
    rxnStage = 2

! In rxnStage 2, we are using algorithm A and
! check to see whether we should switch
else if ((rxnStage == 2).and.(distBrC < rxnProgress2_threshold1)) then
    print *, ""
    print *, " Switching to reaction stage 3"
    print *, ""
    rxnStage = 3

    ! We also want to measure the relative
    ! translational energy for a new baseline
    ! because it most likely shifted
    call getEReWOclusters(coordsIN,momentaIN,ERe_baseline)

    ! One assumption is that this assumes that
    ! the ERe measurement is stable at this
    ! point even before scaling (?).

! In rxnStage 3, we are using algorithm B and
! do not want anymore switching
else
end if

return

end subroutine updateRxnProgress2

subroutine updateRxnProgress3(coordsIN,momentaIN)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(in) :: &
            coordsIN, momentaIN

real(dp) :: distBrC
!real(dp) :: distBrC_threshold1 = 36.0
!real(dp) :: distBrC_threshold2 = 25.0

! To monitor the reaction progress, just look
! at the distance between the carbon and the
! bromine atoms

distBrC = ((coordsIN(4)-coordsIN(10))**2)+&
          ((coordsIN(5)-coordsIN(11))**2)+&
          ((coordsIN(6)-coordsIN(12))**2)

! In rxnStage 0, we are using algorithm *B and
! check to see whether we should switch
if ((rxnStage == 0).and.(distBrC < rxnProgress3_threshold1)) then
    print *, ""
    print *, " Switching to reaction stage 1"
    print *, ""
    rxnStage = 1

    ! When we switch, we turn reject ON
    prev_reject_flag = reject_flag
    reject_flag = .true.

! In rxnStage 1, we are rejecting and
! check to see whether we have reached the
! "middle" region of the intermediate region
else if ((rxnStage == 1).and.(distBrC < rxnProgress3_threshold2)) then
    print *, ""
    print *, " Switching to reaction stage 2"
    print *, ""
    rxnStage = 2

    ! We can turn reject back OFF (if it was off)
    reject_flag = prev_reject_flag

! In rxnStage 2, we are using algorithm *A and
! check to see whether we should switch
else if ((rxnStage == 2).and.(distBrC > &
          (rxnProgress3_threshold2+rxnProgress3_threshold3)/2)) then
    print *, ""
    print *, " Switching to reaction stage 3"
    print *, ""
    rxnStage = 3

    ! When we switch, we turn reject ON
    prev_reject_flag = reject_flag
    reject_flag = .true.

! In rxnStage 3, we are rejecting and
! check to see whether we are in the
! relatively stable "far away" region
else if ((rxnStage == 3).and.(distBrC > rxnProgress3_threshold3)) then
    print *, ""
    print *, " Switching to reaction stage 4"
    print *, ""
    rxnStage = 4

    ! We can turn reject back OFF (if it was off)
    reject_flag = prev_reject_flag

    ! We also want to measure the relative
    ! translational energy for a new baseline
    ! because it most likely shifted
    call getEReWOclusters(coordsIN,momentaIN,ERe_baseline)

    ! One assumption is that this assumes that
    ! the ERe measurement is stable at this
    ! point even before scaling (?).

! In rxnStage 4, we are using algorithm *B and
! do not want anymore switching
else
end if

return

end subroutine updateRxnProgress3

subroutine updateRxnProgress4(coordsIN,momentaIN)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(in) :: &
            coordsIN, momentaIN

real(dp) :: distBrC
!real(dp) :: distBrC_threshold1 = 36.0
!real(dp) :: distBrC_threshold2 = 25.0

! To monitor the reaction progress, just look
! at the distance between the carbon and the
! bromine atoms

distBrC = ((coordsIN(4)-coordsIN(10))**2)+&
          ((coordsIN(5)-coordsIN(11))**2)+&
          ((coordsIN(6)-coordsIN(12))**2)

! In rxnStage 0, we are using algorithm *B and
! check to see whether we should switch
if (rxnStage == 0) then
    if (distBrC < rxnProgress4_threshold1) then
        print *, ""
        print *, " Switching to reaction stage 1"
        print *, ""
        rxnStage = 1

        ! When we switch, we turn reject ON
        prev_reject_flag = reject_flag
        reject_flag = .true.

        ! When we switch, we reset what I call
        ! the "glue counter" which makes sure
        ! that enough glue is applied
        rxnProgress4_thresholdCounter = 0
    end if

! In rxnStage 1, we are rejecting and
! check to see whether we have reached the
! "middle" region of the intermediate region
else if (rxnStage == 1) then

    rxnProgress4_thresholdCounter = &
        rxnProgress4_thresholdCounter + 1

    if (distBrC < rxnProgress4_threshold2) then
        print *, ""
        print *, " Switching to reaction stage 2"
        print *, ""
        rxnStage = 2

        ! We can turn reject back OFF (if it was off)
        reject_flag = prev_reject_flag

        ! When we switch, we reset what I call
        ! the "glue counter" which makes sure
        ! that enough glue is applied
        rxnProgress4_thresholdCounter = 0

    else if ((rxnProgress4_thresholdCounter > &
              rxnProgress4_thresholdCounterMax).and.&
             (distBrC > rxnProgress4_threshold1)) then
        print *, ""
        print *, " Switching to reaction stage 4"
        print *, ""
        rxnStage = 4

        ! We can turn reject back OFF (if it was off)
        reject_flag = prev_reject_flag

        ! We also want to measure the relative
        ! translational energy for a new baseline
        ! because it most likely shifted
        call getEReWOclusters(coordsIN,momentaIN,ERe_baseline)

        ! One assumption is that this assumes that
        ! the ERe measurement is stable at this
        ! point even before scaling (?).

    else
    end if

! In rxnStage 2, we are using algorithm *A and
! check to see whether we should switch
else if ((rxnStage == 2).and.(distBrC > rxnProgress4_threshold2)) then
    print *, ""
    print *, " Switching back to reaction stage 1"
    print *, ""
    rxnStage = 1

    ! When we switch, we turn reject ON
    prev_reject_flag = reject_flag
    reject_flag = .true.

    ! When we switch, we reset what I call
    ! the "glue counter" which makes sure
    ! that enough glue is applied
    rxnProgress4_thresholdCounter = 0

! In rxnStage 4, we are using algorithm *B and
! do not want anymore switching
else

    rxnProgress4_thresholdCounter = &
        rxnProgress4_thresholdCounter + 1

end if

return

end subroutine updateRxnProgress4

subroutine updateRxnProgress5(coordsIN,momentaIN)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: &
            coordsIN, momentaIN

real(dp) :: distHBrCO2

real(dp) :: wa,wb
real(dp),dimension(3) :: pcma, pcmb, qr
real(dp) :: prel,erel

integer :: i

!real(dp) :: distBrC_threshold1 = 36.0
!real(dp) :: distBrC_threshold2 = 25.0

! To monitor the reaction progress, look at:

! (1) the distance between the carbon and the
! bromine atoms

distHBrCO2 = (&
     ((coordsIN(4) - coordsIN(10))**2) + &
     ((coordsIN(5) - coordsIN(11))**2) + &
     ((coordsIN(6) - coordsIN(12))**2) )
!   (((coordsIN(1)+coordsIN(4))/2 - &
!     (coordsIN(7)+coordsIN(10)+coordsIN(13))/3 )**2) + &
!   (((coordsIN(2)+coordsIN(5))/2 - &
!     (coordsIN(8)+coordsIN(11)+coordsIN(14))/3 )**2) + &
!   (((coordsIN(3)+coordsIN(6))/2 - &
!     (coordsIN(9)+coordsIN(12)+coordsIN(15))/3 )**2) )

! (2) the relative translational energy
! between the two clusters

!     RDMASS = SWA * SWB / (SWA + SWB)
!     RCM = 0.0D0
!     DO I=1, 3
!        QR(I) = QCMA(I) - QCMB(I)
!        VR(I) = VCMA(I) - VCMB(I)
!        RCM = RCM + QR(I)**2
!     ENDDO

!     VREL = 0.0D0
!     DO I=1, 3
!        VREL = VREL + VR(I) * QR(I)
!     ENDDO

!     EREL = RDMASS * ((VREL**2) / RCM) / 2.0D0 / C1

!     Etrans = SIGN(EREL,VREL)

wa = masses(1)+masses(2)
wb = masses(3)+masses(4)+masses(5)

pcma(1) = momentaIN(1)+momentaIN(4)
pcma(2) = momentaIN(2)+momentaIN(5)
pcma(3) = momentaIN(3)+momentaIN(6)

pcmb(1) = momentaIN(7)+momentaIN(10)+momentaIN(13)
pcmb(2) = momentaIN(8)+momentaIN(11)+momentaIN(14)
pcmb(3) = momentaIN(9)+momentaIN(12)+momentaIN(15)

!qr(1) = (coordsIN(1)+coordsIN(4)) - &
!        (coordsIN(7)+coordsIN(10)+coordsIN(13))
!qr(2) = (coordsIN(2)+coordsIN(5)) - &
!        (coordsIN(8)+coordsIN(11)+coordsIN(14))
!qr(3) = (coordsIN(3)+coordsIN(6)) - &
!        (coordsIN(9)+coordsIN(12)+coordsIN(15))
qr(1) = (coordsIN(1)*masses(1)+&
         coordsIN(4)*masses(2))      - &
        (coordsIN(7)*masses(3)+&
         coordsIN(10)*masses(4)+&
         coordsIN(13)*masses(5))
qr(2) = (coordsIN(2)*masses(1)+&
         coordsIN(5)*masses(2))      - &
        (coordsIN(8)*masses(3)+&
         coordsIN(11)*masses(4)+&
         coordsIN(14)*masses(5))
qr(3) = (coordsIN(3)*masses(1)+&
         coordsIN(6)*masses(2))      - &
        (coordsIN(9)*masses(3)+&
         coordsIN(12)*masses(4)+&
         coordsIN(15)*masses(5))

prel = (((pcma(1)*wb)-(pcmb(1)*wa))*qr(1) + &
        ((pcma(2)*wb)-(pcmb(2)*wa))*qr(2) + &
        ((pcma(3)*wb)-(pcmb(3)*wa))*qr(3))

erel = sign(&
         ((prel**2) / ((wa+wb)*(wa*wb) * &
           (qr(1)*qr(1) + qr(2)*qr(2) + qr(3)*qr(3)))) / &
          (2.0d0 * C1),prel)


if (rxnStage == 0) then
    if ((distHBrCO2 > rxnProgress5_threshold1) .and. &
        (      erel > rxnProgress5_threshold2)) then
!       (      erel < rxnProgress5_threshold2)) then
        print *, ""
        print *, " Switching to reaction stage 1"
        print *, ""
        rxnStage = 1

        ! When we switch, we turn reject ON
        prev_reject_flag = reject_flag
        reject_flag = .true.

        ! When we switch, we reset what I call
        ! the "glue counter" which makes sure
        ! that enough glue is applied
        rxnProgress5_thresholdCounter = 0
    end if

elseif (rxnStage == 2) then
    if ((distHBrCO2 < rxnProgress5_threshold3) .and. &
        (      erel < 0)) then
!       (      erel < rxnProgress5_threshold2)) then
        print *, ""
        print *, " Switching to reaction stage 1"
        print *, ""
        rxnStage = 1

        ! When we switch, we turn reject ON
        prev_reject_flag = reject_flag
        reject_flag = .true.

        ! When we switch, we reset what I call
        ! the "glue counter" which makes sure
        ! that enough glue is applied
        rxnProgress5_thresholdCounter = 0
    end if

else

    ! Keep a history of the relative
    ! translational energy
    do i = 1, rxnProgress5_thresholdCounterMax-1
        erelCounter(i) = erelCounter(i+1)
    end do
    erelCounter(rxnProgress5_thresholdCounterMax) = erel

    ! Make sure the "glue counter" has elapsed
    if (rxnProgress5_thresholdCounter > &
        rxnProgress5_thresholdCounterMax) then

        ! If the molecules are coming back together,
        ! then just go back to reaction stage 0
        if (erel < 0) then

        print *, ""
        print *, " Switching to reaction stage 0"
        print *, ""
        rxnStage = 0

        reject_flag = prev_reject_flag

        ! If the molecules are dissociating and
        ! the relative translational energy is
        ! starting to flatline, then go to
        ! reaction stage 2
        elseif ((distHBrCO2 > rxnProgress5_threshold5) .and. &
                (sum((erelCounter-(sum(erelCounter)/&
                    rxnProgress5_thresholdCounterMax))**2)/&
                 rxnProgress5_thresholdCounterMax < &
                 rxnProgress5_threshold4)) then

        print *, ""
        print *, " Switching to reaction stage 2"
        rxnStage = 2

        ! We may want to measure the relative
        ! translational energy for a new baseline
        ! because it most likely shifted
        call getEReWOclusters(coordsIN,momentaIN,ERe_baseline)

        print *, "   ERe:", ERe_baseline
        print *, ""

        reject_flag = prev_reject_flag
        end if

    end if

    rxnProgress5_thresholdCounter = &
        rxnProgress5_thresholdCounter + 1

end if

return

end subroutine updateRxnProgress5





subroutine getERe(coordsIN,momentaIN,clusters,ERe)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(in) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters
real(dp),intent(out) :: ERe

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaTCM,momentaTRe,momentaRV
real(dp),dimension(3) :: &
        velocityCMRe,velocityCMRV
real(dp) :: reducedmass

real(dp) :: f1rere,f1rerv,f1rvrv

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities2(coords,momenta,clusters,&
        momentaTCM,momentaTRe,momentaRV,&
        velocityCMRe,velocityCMRV,reducedmass)

f1rere = reducedmass*sum(velocityCMRe**2)/2
f1rerv = reducedmass*dot_product(velocityCMRe,velocityCMRV)
f1rvrv = reducedmass*sum(velocityCMRV**2)/2

ERe = (f1rere + f1rerv + f1rvrv) / C1

return

end subroutine getERe

subroutine getEReWOclusters(coordsIN,momentaIN,ERe)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(in) :: momentaIN
integer,dimension(1,Natoms) :: clusters1
integer,dimension(2,Natoms) :: clusters2
real(dp),intent(out) :: ERe

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaTCM,momentaTRe,momentaRV
real(dp),dimension(3) :: &
        velocityCMRe,velocityCMRV
real(dp) :: reducedmass

real(dp) :: f1rere,f1rerv,f1rvrv

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call clusterAtoms(coords,clusters1,clusters2)

call decomposeVelocities2(coords,momenta,clusters2,&
        momentaTCM,momentaTRe,momentaRV,&
        velocityCMRe,velocityCMRV,reducedmass)

f1rere = reducedmass*sum(velocityCMRe**2)/2
f1rerv = reducedmass*dot_product(velocityCMRe,velocityCMRV)
f1rvrv = reducedmass*sum(velocityCMRV**2)/2

ERe = (f1rere + f1rerv + f1rvrv) / C1

return

end subroutine getEReWOclusters




subroutine scaleVelocities1(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaT,momentaR,momentaV

real(dp) :: scalingThreshold = 1.0d-6
real(dp) :: scaling
real(dp) :: scalingT,scalingValue
real(dp) :: scalingtt,scalingrr,scalingvv
real(dp) :: scalingtr,scalingtv,scalingrv
real(dp) :: Ktt,Krr,Kvv
real(dp) :: Ktr,Ktv,Krv

integer :: scalingIterations_max = 100
integer :: i

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities1(coords,momenta,clusters,&
        momentaT,momentaR,momentaV)

Ktt = 0.0d0
Krr = 0.0d0
Kvv = 0.0d0
Ktr = 0.0d0
Ktv = 0.0d0
Krv = 0.0d0
do i = 1, Natoms
    Ktt = Ktt + 0.5d0*sum(momentaT(:,i)**2)/masses(i)
    Krr = Krr + 0.5d0*sum(momentaR(:,i)**2)/masses(i)
    Kvv = Kvv + 0.5d0*sum(momentaV(:,i)**2)/masses(i)
    Ktr = Ktr + dot_product(momentaT(:,i),momentaR(:,i))/masses(i)
    Ktv = Ktv + dot_product(momentaT(:,i),momentaV(:,i))/masses(i)
    Krv = Krv + dot_product(momentaR(:,i),momentaV(:,i))/masses(i)
end do

scaling = 1.0d0
do i = 1, scalingIterations_max
    scalingtt = scaling**(2*alphat)
    scalingrr = scaling**(2*alphar)
    scalingvv = scaling**(2*alphav)
    scalingtr = scaling**(alphat+alphar)
    scalingtv = scaling**(alphat+alphav)
    scalingrv = scaling**(alphar+alphav)

    scalingT = scalingtt*Ktt + scalingrr*Krr + scalingvv*Kvv + &
               scalingtr*Ktr + scalingtv*Ktv + scalingrv*Krv
    scalingValue = (C1 * (H_baseline - V) / scalingT ) - 1.0d0

    print *, "iteration", i
    print *, "    scaling:", scaling
    print *, "      Value:", scalingValue

    if (abs(scalingValue) < scalingThreshold) exit

    scaling = scaling + scalingValue / &
                    ((C1 * (H_baseline - V) / (scalingT**2)) * &
                     (2*alphat*scalingtt*Ktt + &
                      2*alphar*scalingrr*Krr + &
                      2*alphav*scalingvv*Kvv + &
                      (alphat+alphar)*scalingtr*Ktr + &
                      (alphat+alphav)*scalingtv*Ktv + &
                      (alphar+alphav)*scalingrv*Krv) / scaling)

    if (scaling < 0) scaling = 0.0d0
end do

if (i < scalingIterations_max) then

    momenta = momentaT*scaling**(alphat) + &
              momentaR*scaling**(alphar) + &
              momentaV*scaling**(alphav)

    momentaIN = reshape(momenta,(/3*Natoms/))

end if

return

end subroutine scaleVelocities1

subroutine scaleVelocities1plus(coordsIN,&
                       momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaT,momentaR,momentaV

real(dp) :: scalingThreshold = 1.0d-6
real(dp) :: scalingT,scaling,f1,f2
real(dp) :: scalingr,scalingv
real(dp) :: scalingtt,scalingrr,scalingvv
real(dp) :: scalingtr,scalingtv,scalingrv
real(dp) :: f1tt,f1rr,f1vv
real(dp) :: f1tr,f1tv,f1rv
real(dp) :: f2tt,f2rr,f2vv
real(dp) :: f2tr,f2tv,f2rv
real(dp) :: totalMass
real(dp) :: DFa,DFb,DFc,DFd,detDF
real(dp) :: invDFa,invDFb,invDFc,invDFd,invdetDF
real(dp) :: DFthreshold = 1.0d-20
real(dp) :: minscalingT,minscaling,minf1,minf2

real(dp),dimension(3) :: ptt,prr,pvv

integer :: scalingIterations_max = 100
integer :: i

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities1(coords,momenta,clusters,&
        momentaT,momentaR,momentaV)

f1tt = 0.0d0
f1rr = 0.0d0
f1vv = 0.0d0
f1tr = 0.0d0
f1tv = 0.0d0
f1rv = 0.0d0
do i = 1, Natoms
    f1tt = f1tt + 0.5d0*sum(momentaT(:,i)**2)/masses(i)
    f1rr = f1rr + 0.5d0*sum(momentaR(:,i)**2)/masses(i)
    f1vv = f1vv + 0.5d0*sum(momentaV(:,i)**2)/masses(i)
    f1tr = f1tr + dot_product(momentaT(:,i),momentaR(:,i))/masses(i)
    f1tv = f1tv + dot_product(momentaT(:,i),momentaV(:,i))/masses(i)
    f1rv = f1rv + dot_product(momentaR(:,i),momentaV(:,i))/masses(i)
end do

f1tt = f1tt / C1
f1rr = f1rr / C1
f1vv = f1vv / C1
f1tr = f1tr / C1
f1tv = f1tv / C1
f1rv = f1rv / C1

ptt = 0.0d0
prr = 0.0d0
pvv = 0.0d0
totalMass = 0.0d0
do i = 1, Natoms
    ptt = ptt + momentaT(:,i)
    prr = prr + momentaR(:,i)
    pvv = pvv + momentaV(:,i)
    totalMass = totalMass + masses(i)
end do

f2tt = 0.5d0*sum(ptt**2)/(totalMass * C1)
f2rr = 0.5d0*sum(prr**2)/(totalMass * C1)
f2vv = 0.5d0*sum(pvv**2)/(totalMass * C1)
f2tr = dot_product(ptt,prr)/(totalMass * C1)
f2tv = dot_product(ptt,pvv)/(totalMass * C1)
f2rv = dot_product(prr,pvv)/(totalMass * C1)

scalingT = 1.0d0
scaling = 1.0d0
do i = 1, scalingIterations_max
    scalingr = scaling**alphar
    scalingv = scaling**alphav

    scalingtt = scalingT**2
    scalingrr = scalingr**2
    scalingvv = scalingv**2
    scalingtr = scalingT*scalingr
    scalingtv = scalingT*scalingv
    scalingrv = scalingr*scalingv

    f1 = scalingtt*f1tt + scalingrr*f1rr + scalingvv*f1vv + &
         scalingtr*f1tr + scalingtv*f1tv + scalingrv*f1rv - &
         (H_baseline - V)
    f2 = scalingtt*f2tt + scalingrr*f2rr + scalingvv*f2vv + &
         scalingtr*f2tr + scalingtv*f2tv + scalingrv*f2rv

    print *, "iteration", i
    print *, "    scaling:", scalingT, scaling
    print *, "      Value:", f1, f2

    if ((abs(f1) < scalingThreshold) .and.&
        (abs(f2) < scalingThreshold)) then
        minscalingT = scalingT
        minscaling = scaling
        minf1 = f1
        minf2 = f2
        exit
    end if

    if (abs(f2) < minf2) then
        if (abs(f1) < 0.1d0) then
            minscalingT = scalingT
            minscaling = scaling
            minf1 = f1
            minf2 = f2
        end if
    end if

    DFa = 2*scalingT*f1tt + scalingr*f1tr + scalingv*f1tv
    DFb = (2*alphaR*f1rr + 2*alphaV*f1vv + (alphaR+alphaV)*f1rv +&
           scalingT*(alphaR*f1tr + alphaV*f1tv)) / scaling
    DFc = 2*scalingT*f2tt + scalingr*f2tr + scalingv*f2tv
    DFd = (2*alphaR*f2rr + 2*alphaV*f2vv + (alphaR+alphaV)*f2rv +&
           scalingT*(alphaR*f2tr + alphaV*f2tv)) / scaling

    detDF = DFa*DFd - DFb*DFc
    if (abs(detDF) < DFthreshold) exit

    invdetDF = 1.0d0 / detDF
    invDFa = DFd*invdetDF
    invDFb = -DFb*invdetDF
    invDFc = -DFc*invdetDF
    invDFd = DFa*invdetDF

    scalingT = scalingT - (f1*invDFa + f2*invDFb)
    scaling = scaling - (f1*invDFc + f2*invDFd)

    if ((scalingT < 0).or.(scaling < 0)) exit
end do

if (abs(minf2) < 0.1d0) then

    momenta = momentaT*(minscalingT) + &
              momentaR*(minscaling**alphar) + &
              momentaV*(minscaling**alphav)

    momentaIN = reshape(momenta,(/3*Natoms/))

end if

return

end subroutine scaleVelocities1plus

subroutine scaleVelocities2(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaTCM,momentaTRe,momentaRV
real(dp),dimension(3) :: &
        velocityCMRe,velocityCMRV
real(dp) :: reducedmass

real(dp) :: scalingThreshold = 1.0d-12
real(dp) :: scalingTRe,scalingRV
real(dp) :: rere,rerv,rvrv

real(dp) :: f1,f2
real(dp) :: f1rere,f1rerv,f1rvrv
real(dp) :: f2rere,f2rerv,f2rvrv

real(dp) :: minTRe,minRV
real(dp) :: minf1,minf2

real(dp) :: DFa, DFb, DFc, DFd
real(dp) :: detDF,invdetDF
real(dp) :: DFthreshold = 1.0d-20
real(dp) :: invDFa, invDFb, invDFc, invDFd

integer :: scalingIterations_max = 100
integer :: i, j

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities2(coords,momenta,clusters,&
        momentaTCM,momentaTRe,momentaRV,&
        velocityCMRe,velocityCMRV,reducedmass)

f1rere = reducedmass*sum(velocityCMRe**2)/(2*C1)
f1rerv = reducedmass*dot_product(velocityCMRe,velocityCMRV)/C1
f1rvrv = reducedmass*sum(velocityCMRV**2)/(2*C1)

f2rere = 0.0d0
f2rerv = 0.0d0
f2rvrv = 0.0d0
do i = 1, Natoms
    f2rere = f2rere + 0.50d0*sum(momentaTRe(:,i)**2)/masses(i)
    f2rerv = f2rerv + dot_product(momentaTRe(:,i),momentaTCM(:,i)+momentaRV(:,i))/masses(i)
    f2rvrv = f2rvrv + 0.50d0*sum((momentaTCM(:,i)+momentaRV(:,i))**2)/masses(i)
end do

f2rere = f2rere / C1
f2rerv = f2rerv / C1
f2rvrv = f2rvrv / C1

minTRe = 1.0d0
minRV = 1.0d0
minf1 = 1.0d9
minf2 = 1.0d9

scalingTRe = 1.0d0
scalingRV = 1.0d0
do i = 1, scalingIterations_max
    rere = scalingTRe**2
    rerv = scalingTRe*scalingRV
    rvrv = scalingRV**2

    f1 = rere*f1rere + rerv*f1rerv + rvrv*f1rvrv - ERe_baseline
    f2 = rere*f2rere + rerv*f2rerv + rvrv*f2rvrv - (H_baseline - V)

    print *, "iteration", i
    print *, "    scaling:", scalingTRe, scalingRV
    print *, "      Value:", f1, f2
!   call getERe(reshape(coords,(/ 3*Natoms /)),&
!               reshape(momentaTRe*scalingTRe + &
!                       (momentaTCM+momentaRV)*scalingRV,(/ 3*Natoms /)),&
!               invdetDF)
!   print *, "        ERe:", invdetDF
!   print *, "   baseline:", ERe_baseline

    if ((abs(f1) < scalingThreshold) .and.&
        (abs(f2) < scalingThreshold)) then
        minTRe = scalingTRe
        minRV = scalingRV
        minf1 = f1
        minf2 = f2
        exit
    end if

    if (abs(f2) < minf2) then
        if (abs(f1) < 0.1d0) then
            minTRe = scalingTRe
            minRV = scalingRV
            minf1 = f1
            minf2 = f2
        end if
    end if

    DFa = 2*scalingTRe*f1rere + scalingRV*f1rerv
    DFb = scalingTRe*f1rerv + 2*scalingRV*f1rvrv
    DFc = 2*scalingTRe*f2rere + scalingRV*f2rerv
    DFd = scalingTRe*f2rerv + 2*scalingRV*f2rvrv

    detDF = DFa*DFd - DFb*DFc
    if (abs(detDF) < DFthreshold) exit

    invdetDF = 1.0d0 / detDF
    invDFa = DFd*invdetDF
    invDFb = -DFb*invdetDF
    invDFc = -DFc*invdetDF
    invDFd = DFa*invdetDF

    scalingTRe = scalingTRe - (f1*invDFa + f2*invDFb)
    scalingRV = scalingRV - (f1*invDFc + f2*invDFd)

!   if (scalingTre < 0) scalingTRe = 0.0d0
!   if (scalingRV < 0) scalingRV = 0.0d0

    if ((scalingTre < 0).or.(scalingRV < 0)) then
        exit
    end if
end do

!if (i < scalingIterations_max) then
if (abs(minf2) < 0.1d0) then

    momenta = momentaTRe*minTRe + &
              (momentaTCM+momentaRV)*minRV

    momentaIN = reshape(momenta,(/3*Natoms/))

end if

return

end subroutine scaleVelocities2

subroutine scaleVelocities2plus(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaTCM,momentaTRe,momentaRV
real(dp),dimension(3) :: &
        velocityTCM,velocityTRe,velocityCMRV
real(dp) :: reducedmass

real(dp) :: scalingThreshold = 1.0d-12
real(dp) :: scalingTCM,scalingTRe,scalingRV
real(dp) :: cmcm,cmre,cmrv,rere,rerv,rvrv

real(dp) :: f1,f2,f3
real(dp) :: f1cmcm,f1cmre,f1cmrv,f1rere,f1rerv,f1rvrv
real(dp) :: f2cmcm,f2cmre,f2cmrv,f2rere,f2rerv,f2rvrv
real(dp) :: f3cmcm,f3cmre,f3cmrv,f3rere,f3rerv,f3rvrv

real(dp) :: minTCM,minTRe,minRV
real(dp) :: minf1,minf2,minf3

real(dp) :: DFa, DFb, DFc
real(dp) :: DFd, DFe, DFf
real(dp) :: DFg, DFh, DFi
real(dp) :: detDF,invdetDF
real(dp) :: DFthreshold = 1.0d-20
real(dp) :: invDFa, invDFb, invDFc
real(dp) :: invDFd, invDFe, invDFf
real(dp) :: invDFg, invDFh, invDFi

integer :: scalingIterations_max = 100
integer :: i, j

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities2plus(coords,momenta,clusters,&
        momentaTCM,momentaTRe,momentaRV,&
        velocityTCM,velocityTRe,velocityCMRV,reducedmass)

f1cmcm = reducedmass*sum(velocityTCM**2)/(2*C1)
f1rere = reducedmass*sum(velocityTRe**2)/(2*C1)
f1rvrv = reducedmass*sum(velocityCMRV**2)/(2*C1)
f1cmre = reducedmass*dot_product(velocityTCM,velocityTRe)/C1
f1cmrv = reducedmass*dot_product(velocityTCM,velocityCMRV)/C1
f1rerv = reducedmass*dot_product(velocityTRe,velocityCMRV)/C1

f2cmcm = 0.0d0
f2cmre = 0.0d0
f2cmrv = 0.0d0
f2rere = 0.0d0
f2rerv = 0.0d0
f2rvrv = 0.0d0
do i = 1, Natoms
    f2cmcm = f2cmcm + 0.50d0*sum(momentaTCM(:,i)**2)/masses(i)
    f2rere = f2rere + 0.50d0*sum(momentaTRe(:,i)**2)/masses(i)
    f2rvrv = f2rvrv + 0.50d0*sum(momentaRV(:,i)**2)/masses(i)
    f2cmre = f2cmre + dot_product(momentaTCM(:,i),momentaTRe(:,i))/masses(i)
    f2cmrv = f2cmrv + dot_product(momentaTCM(:,i),momentaRV(:,i))/masses(i)
    f2rerv = f2rerv + dot_product(momentaTRe(:,i),momentaRV(:,i))/masses(i)
end do

f2cmcm = f2cmcm / C1
f2cmre = f2cmre / C1
f2cmrv = f2cmrv / C1
f2rere = f2rere / C1
f2rerv = f2rerv / C1
f2rvrv = f2rvrv / C1

minTCM = 1.0d0
minTRe = 1.0d0
minRV = 1.0d0
minf1 = 1.0d9
minf2 = 1.0d9
minf3 = 1.0d9

scalingTCM = 1.0d0
scalingTRe = 1.0d0
scalingRV = 1.0d0
do i = 1, scalingIterations_max
    cmcm = scalingTCM**2
    rere = scalingTRe**2
    rvrv = scalingRV**2
    cmre = scalingTCM*scalingTRe
    cmrv = scalingTCM*scalingRV
    rerv = scalingTRe*scalingRV

    f1 = cmcm*f1cmcm + rere*f1rere + rvrv*f1rvrv + &
         cmre*f1cmre + cmrv*f1cmrv + rerv*f1rerv - &
         ERe_baseline
    f2 = cmcm*f2cmcm + rere*f2rere + rvrv*f2rvrv + &
         cmre*f2cmre + cmrv*f2cmrv + rerv*f2rerv - &
         (H_baseline - V)
    f3 = cmcm*f3cmcm + rere*f3rere + rvrv*f3rvrv + &
         cmre*f3cmre + cmrv*f3cmrv + rerv*f3rerv

    print *, "iteration", i
    print *, "    scaling:", scalingTCM, scalingTRe, scalingRV
    print *, "      Value:", f1, f2, f3
!   call getERe(reshape(coords,(/ 3*Natoms /)),&
!               reshape(momentaTRe*scalingTRe + &
!                       (momentaTCM+momentaRV)*scalingRV,(/ 3*Natoms /)),&
!               invdetDF)
!   print *, "        ERe:", invdetDF
!   print *, "   baseline:", ERe_baseline

    if ((abs(f1) < scalingThreshold) .and.&
        (abs(f2) < scalingThreshold)) then
        minTCM = scalingTCM
        minTRe = scalingTRe
        minRV = scalingRV
        minf1 = f1
        minf2 = f2
        minf3 = f3
        exit
    end if

    if (abs(f2) < minf2) then
        if ((abs(f1) < 0.1d0).and.(abs(f3) < 0.1d0)) then
            minTCM = scalingTCM
            minTRe = scalingTRe
            minRV = scalingRV
            minf1 = f1
            minf2 = f2
            minf3 = f3
        end if
    end if

    ! Matrix is of this form:
    !
    !   A    B    C    - f1
    !
    !   D    E    F    - f2
    !
    !   G    H    I    - f3
    !
    !   |    |    |
    !  TCM  TRe   RV
    !

    DFa = 2*scalingTCM*f1cmcm +   scalingTRe*f1cmre +   scalingRV*f1cmrv
    DFb =   scalingTCM*f1cmre + 2*scalingTRe*f1rere +   scalingRV*f1rerv
    DFc =   scalingTCM*f1cmrv +   scalingTRe*f1rerv + 2*scalingRV*f1rvrv
    DFd = 2*scalingTCM*f2cmcm +   scalingTRe*f2cmre +   scalingRV*f2cmrv
    DFe =   scalingTCM*f2cmre + 2*scalingTRe*f2rere +   scalingRV*f2rerv
    DFf =   scalingTCM*f2cmrv +   scalingTRe*f2rerv + 2*scalingRV*f2rvrv
    DFg = 2*scalingTCM*f3cmcm +   scalingTRe*f3cmre +   scalingRV*f3cmrv
    DFh =   scalingTCM*f3cmre + 2*scalingTRe*f3rere +   scalingRV*f3rerv
    DFi =   scalingTCM*f3cmrv +   scalingTRe*f3rerv + 2*scalingRV*f3rvrv

    detDF = DFa*(DFe*DFi-DFf*DFh) - &
            DFb*(DFd*DFi-DFf*DFg) + &
            DFc*(DFd*DFh-DFe*DFg)
    if (abs(detDF) < DFthreshold) exit

    invdetDF = 1.0d0 / detDF

    invDFa = invdetDF*(DFe*DFi-DFf*DFh)
    invDFb = invdetDF*(DFf*DFg-DFd*DFi)
    invDFc = invdetDF*(DFd*DFh-DFe*DFg)
    invDFd = invdetDF*(DFc*DFh-DFb*DFi)
    invDFe = invdetDF*(DFa*DFi-DFc*DFg)
    invDFf = invdetDF*(DFb*DFg-DFa*DFh)
    invDFg = invdetDF*(DFb*DFf-DFc*DFe)
    invDFh = invdetDF*(DFc*DFd-DFa*DFf)
    invDFi = invdetDF*(DFa*DFe-DFb*DFd)

    scalingTCM = scalingTCM - (f1*invDFa + f2*invDFb + f3*invDFc)
    scalingTRe = scalingTRe - (f1*invDFd + f2*invDFe + f3*invDFf)
    scalingRV = scalingRV - (f1*invDFg + f2*invDFh + f3*invDFi)

!   if (scalingTre < 0) scalingTRe = 0.0d0
!   if (scalingRV < 0) scalingRV = 0.0d0

    if ((scalingTCM < 0).or.&
        (scalingTre < 0).or.&
        (scalingRV < 0)) then
        exit
    end if
end do

!if (i < scalingIterations_max) then
if ((abs(minf2) < 0.1d0).and.&
    (abs(minf3) < 0.1d0)) then

    momenta = momentaTRe*minTRe + &
              (momentaTCM+momentaRV)*minRV

    momentaIN = reshape(momenta,(/3*Natoms/))

end if

return

end subroutine scaleVelocities2plus

subroutine scaleVelocities3(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3,Natoms) :: &
        momentaTCM,momentaTRe,momentaRV
real(dp),dimension(3) :: &
        velocityCMRe,velocityCMRV
real(dp) :: reducedmass

real(dp) :: scalingThreshold = 1.0d-12
real(dp) :: scalingTRe,scalingRV

real(dp) :: f1,f2

real(dp) :: DFa, DFb, DFc, DFd
real(dp) :: invdetDF
real(dp) :: invDFa, invDFb, invDFc, invDFd

integer :: scalingIterations_max = 100
integer :: i, j

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities2(coords,momenta,clusters,&
        momentaTCM,momentaTRe,momentaRV,&
        velocityCMRe,velocityCMRV,reducedmass)

scalingTRe = 1.0d0
scalingRV = 1.0d0
do i = 1, scalingIterations_max
    momenta = momentaTRe*scalingTRe + &
              (momentaTCM+momentaRV)*scalingRV

    call getERe(reshape(coords,(/ 3*Natoms /)),&
                reshape(momenta,(/ 3*Natoms /)),&
                clusters,f1)
    f1 = f1 - ERe_baseline

    f2 = 0.0d0
    do j = 1, Natoms
        f2 = f2 + 0.5d0*sum((momenta(:,j))**2)/masses(j)
    end do
    f2 = (f2 / C1) - (H_baseline - V)

    print *, "iteration", i
    print *, "    scaling:", scalingTRe, scalingRV
    print *, "      Value:", f1, f2
    print *, "        ERe:", f1 + ERe_baseline
    print *, "   baseline:", ERe_baseline

    if ((abs(f1) < scalingThreshold) .and.&
        (abs(f2) < scalingThreshold)) exit

    call getDEReDTRe(momentaTRe,momentaTCM,momentaRV,&
                     clusters,scalingTRe,scalingRV,coords,DFa)
    call getDEReDRV(momentaTRe,momentaTCM,momentaRV,&
                    clusters,scalingTRe,scalingRV,coords,DFb)
    call getDHDTRe(momentaTRe,momentaTCM,momentaRV,&
                   scalingTRe,scalingRV,DFc)
    call getDHDRV(momentaTRe,momentaTCM,momentaRV,&
                  scalingTRe,scalingRV,DFd)

    print *, "DFa:", DFa
    print *, "DFb:", DFb
    print *, "DFc:", DFc
    print *, "DFd:", DFd

    invdetDF = 1.0d0 / (DFa*DFd - DFb*DFc)
    invDFa = DFd*invdetDF
    invDFb = -DFb*invdetDF
    invDFc = -DFc*invdetDF
    invDFd = DFa*invdetDF

    scalingTRe = scalingTRe - (f1*invDFa + f2*invDFb)
    scalingRV = scalingRV - (f1*invDFc + f2*invDFd)

    if (scalingTre < 0) scalingTRe = 0.0d0
    if (scalingRV < 0) scalingRV = 0.0d0
end do

if (i < scalingIterations_max) then
    momentaIN = reshape(momenta,(/3*Natoms/))
end if

return

end subroutine scaleVelocities3

subroutine getDEReDTRe(momentaTRe,momentaTCM,momentaRV,&
                       clusters,scalingTRe,scalingRV,&
                       coords,DEReDTRe)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3,Natoms),intent(in) :: &
        coords,momentaTRe,momentaTCM,momentaRV
integer,dimension(2,Natoms),intent(in) :: clusters
real(dp),intent(in) :: scalingTRe,scalingRV
real(dp),intent(out) :: DEReDTRe

real(dp),dimension(3,Natoms) :: momenta

real(dp) :: EReMINUS,ERePLUS

real(dp) :: DTRe = 1.0d-6

print *, "mom  TRe(1):",momentaTRe(:,1)

momenta = momentaTRe*(scalingTRe-DTRe) + &
          (momentaTCM+momentaRV)*scalingRV
print *, "momMINUS(1):",momenta(:,1)

call getERe(reshape(coords,(/ 3*Natoms /)),&
            reshape(momenta,(/ 3*Natoms /)),&
            clusters,EReMINUS)

momenta = momentaTRe*(scalingTRe+DTRe) + &
          (momentaTCM+momentaRV)*scalingRV
print *, " momPLUS:",momenta(:,1)

call getERe(reshape(coords,(/ 3*Natoms /)),&
            reshape(momenta,(/ 3*Natoms /)),&
            clusters,ERePLUS)

DEReDTRe = (ERePLUS - EReMINUS)/(2*DTRe)

return

end subroutine getDEReDTRe

subroutine getDEReDRV(momentaTRe,momentaTCM,momentaRV,&
                      clusters,scalingTRe,scalingRV,&
                      coords,DEReDRV)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3,Natoms),intent(in) :: &
        coords,momentaTRe,momentaTCM,momentaRV
integer,dimension(2,Natoms),intent(in) :: clusters
real(dp),intent(in) :: scalingTRe,scalingRV
real(dp),intent(out) :: DEReDRV

real(dp),dimension(3,Natoms) :: momenta

real(dp) :: EReMINUS,ERePLUS

real(dp) :: DRV = 1.0d-12

momenta = momentaTRe*scalingTRe + &
          (momentaTCM+momentaRV)*(scalingRV-DRV)

call getERe(reshape(coords,(/ 3*Natoms /)),&
            reshape(momenta,(/ 3*Natoms /)),&
            clusters,EReMINUS)

momenta = momentaTRe*scalingTRe + &
          (momentaTCM+momentaRV)*(scalingRV+DRV)

call getERe(reshape(coords,(/ 3*Natoms /)),&
            reshape(momenta,(/ 3*Natoms /)),&
            clusters,ERePLUS)

DEReDRV = (ERePLUS - EReMINUS)/(2*DRV)

return

end subroutine getDEReDRV

subroutine getDHDTRe(momentaTRe,momentaTCM,momentaRV,&
                     scalingTRe,scalingRV,&
                     DHDTRe)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: C1

real(dp),dimension(3,Natoms),intent(in) :: &
        momentaTRe,momentaTCM,momentaRV
real(dp),intent(in) :: scalingTRe,scalingRV
real(dp),intent(out) :: DHDTRe

real(dp),dimension(3,Natoms) :: momenta

real(dp) :: HMINUS,HPLUS

real(dp) :: DTRe = 1.0d-12

integer :: i

momenta = momentaTRe*(scalingTRe-DTRe) + &
          (momentaTCM+momentaRV)*scalingRV

HMINUS = 0.0d0
do i = 1, Natoms
    HMINUS = HMINUS + 0.5d0*sum(momenta(:,i)**2)/masses(i)
end do

momenta = momentaTRe*(scalingTRe+DTRe) + &
          (momentaTCM+momentaRV)*scalingRV

HPLUS = 0.0d0
do i = 1, Natoms
    HPLUS = HPLUS + 0.5d0*sum(momenta(:,i)**2)/masses(i)
end do

DHDTRe = (HPLUS - HMINUS)/(2*C1*DTRe)

return

end subroutine getDHDTRe

subroutine getDHDRV(momentaTRe,momentaTCM,momentaRV,&
                    scalingTRe,scalingRV,&
                    DHDRV)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: C1

real(dp),dimension(3,Natoms),intent(in) :: &
        momentaTRe,momentaTCM,momentaRV
real(dp),intent(in) :: scalingTRe,scalingRV
real(dp),intent(out) :: DHDRV

real(dp),dimension(3,Natoms) :: momenta

real(dp) :: HMINUS,HPLUS

real(dp) :: DRV = 1.0d-12

integer :: i

momenta = momentaTRe*scalingTRe + &
          (momentaTCM+momentaRV)*(scalingRV-DRV)

HMINUS = 0.0d0
do i = 1, Natoms
    HMINUS = HMINUS + 0.5d0*sum(momenta(:,i)**2)/masses(i)
end do

momenta = momentaTRe*scalingTRe + &
          (momentaTCM+momentaRV)*(scalingRV+DRV)

HPLUS = 0.0d0
do i = 1, Natoms
    HPLUS = HPLUS + 0.5d0*sum(momenta(:,i)**2)/masses(i)
end do

DHDRV = (HPLUS - HMINUS)/(2*C1*DRV)

return

end subroutine getDHDRV


subroutine scaleVelocities4(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta
real(dp),dimension(3) :: coordsCM

real(dp),dimension(3,Natoms) :: &
        momentaT,momentaR,momentaV
real(dp),dimension(3,3,Natoms) :: momentaTRV

real(dp),dimension(3,3) :: E1,E2,E3
real(dp) :: linear_kcal_ratio = 1.0d0
real(dp) :: angular_kcal_ratio = 1.0d0

real(dp),dimension(3,3) :: L11, L112
real(dp),dimension(3) :: L1
real(dp) :: L0,L2

real(dp),dimension(3,3) :: L11L11T,L112L112T

real(dp) :: gammaT,gammaR,gammaV,lambda
real(dp) :: DLgamma, DLlambda
real(dp) :: mingammaT,mingammaR,mingammaV,minlambda
real(dp) :: minDLgamma, minDLlambda
real(dp),dimension(3,1) :: gammavector

real(dp) :: gammaThreshold = 1.0d-5
real(dp) :: lambdaThreshold = 1.0d-5

real(dp),dimension(4) :: DL
real(dp),dimension(4,4) :: DDL, invDDL
real(dp) :: detDDL, invdetDDL
real(dp) :: DDLThreshold = 1.0d-7

integer :: scalingIterations_max = 100
integer :: i, j, m, n

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

!call decomposeVelocities2plus(coords,momenta,clusters,&
!        momentaTCM,momentaTRe,momentaRV,&
!        velocityTCM,velocityTRe,velocityCMRV,reducedmass)
call decomposeVelocities1(coords,&
             momenta,clusters,&
             momentaT,momentaR,momentaV)

do i = 1, Natoms
    momentaTRV(:,1,i) = momentaT(:,i)
    momentaTRV(:,2,i) = momentaR(:,i)
    momentaTRV(:,3,i) = momentaV(:,i)
end do

coordsCM = 0.0d0
do i = 1, Natoms
    coordsCM = coordsCM + masses(i)*coords(:,i)
end do
coordsCM = coordsCM / sum(masses(1:Natoms))

E1 = 0.0d0
E2 = 0.0d0
E3 = 0.0d0
do i = 1, Natoms
    coords(:,i) = coords(:,i) - coordsCM
    E1 = E1 + matmul(transpose(momentaTRV(:,:,i)),momentaTRV(:,:,i))/masses(i)
    E2 = E2 + momentaTRV(:,:,i)
    E3(1,:) = E3(1,:) + coords(2,i)*momentaTRV(3,:,i) &
                      - coords(3,i)*momentaTRV(2,:,i)
    E3(2,:) = E3(2,:) + coords(3,i)*momentaTRV(1,:,i) &
                      - coords(1,i)*momentaTRV(3,:,i)
    E3(3,:) = E3(3,:) + coords(1,i)*momentaTRV(2,:,i) &
                      - coords(2,i)*momentaTRV(1,:,i)
end do
E1 = E1 / C1
E2 = E2 * linear_kcal_ratio
E3 = E3 * angular_kcal_ratio

L11 = matmul(transpose(E2),E2) + matmul(transpose(E3),E3)
L1 = -2*reshape(matmul(reshape(L_baseline,(/ 1, 3 /)),E3),(/ 3 /))
L0 = sum(L_baseline**2)
L112 = E1
L2 = V - H_baseline

L11L11T = L11 + transpose(L11)
L112L112T = L112 + transpose(L112)

gammaT = 1.0d0
gammaR = 1.0d0
gammaV = 1.0d0
lambda = V - H_baseline

mingammaT = 1.0d0
mingammaR = 1.0d0
mingammaV = 1.0d0
minlambda = 1.0d9
minDLgamma = 1.0d9
minDLlambda = 1.0d9
do i = 1, scalingIterations_max

    gammavector(1,1) = gammaT
    gammavector(2,1) = gammaR
    gammavector(3,1) = gammaV

    DL(1:3) = reshape(matmul(lambda*L112L112T + L11L11T,&
                      gammavector),(/ 3 /)) + L1
!                    gammavector) + transpose(L1)
    DL(4) = dot_product(reshape(gammavector,(/ 3 /)),&
                   reshape(matmul(L112,gammavector),(/ 3 /))) + L2

    print *, "iteration", i
    print *, "    gamma:", gammaT, gammaR, gammaV
    print *, "       DL:", DL(1:4)

    DLgamma = sum(DL(1:3)**2)
    DLlambda = abs(DL(4))
    if ((DLgamma < gammaThreshold).and.&
        (DLlambda < lambdaThreshold)) then
        mingammaT = gammaT
        mingammaR = gammaR
        mingammaV = gammaV
        minlambda = lambda
        minDLgamma = DLgamma
        minDLlambda = DLlambda
        exit
    end if

    if (DLgamma < minDLgamma) then
        if (DLlambda < minDLlambda) then
            mingammaT = gammaT
            mingammaR = gammaR
            mingammaV = gammaV
            minlambda = lambda
            minDLgamma = DLgamma
            minDLlambda = DLlambda
        end if
    end if

    ! Matrix is of this form:
    !
    !   11   12   13   14   - DL1
    !
    !   21   22   23   24   - DL2
    !
    !   31   32   33   34   - DL3
    !
    !   41   42   43   44   - DL4
    !
    !   |    |    |    |
    !   gT   gR   gV   l
    !

    DDL(1:3,1:3) = lambda*L112L112T + L11L11T
    DDL(1:3,4) = reshape(matmul(L112L112T,gammavector),(/ 3 /))
    DDL(4,1:3) = DDL(1:3,4)
    DDL(4,4) = 0.0d0

    detDDL = 0.0d0
    do m = 1, 4
    do n = 1, 4
        invDDL(m,n) = cofactor3of4(DDL,m,n)
        if (m == 1) detDDL = detDDL + invDDL(m,n) * DDL(m,n)
    end do
    end do
    if (abs(detDDL) < DDLthreshold) exit
    invdetDDL = 1.0d0 / detDDL

    gammaT = gammaT - (&
             invDDL(1,1)*DL(1)+invDDL(1,2)*DL(2)+&
             invDDL(1,3)*DL(3)+invDDL(1,4)*DL(4)) * invdetDDL
    gammaR = gammaR - (&
             invDDL(2,1)*DL(1)+invDDL(2,2)*DL(2)+&
             invDDL(2,3)*DL(3)+invDDL(2,4)*DL(4)) * invdetDDL
    gammaV = gammaV - (&
             invDDL(3,1)*DL(1)+invDDL(3,2)*DL(2)+&
             invDDL(3,3)*DL(3)+invDDL(3,4)*DL(4)) * invdetDDL
    lambda = lambda - (&
             invDDL(4,1)*DL(1)+invDDL(4,2)*DL(2)+&
             invDDL(4,3)*DL(3)+invDDL(4,4)*DL(4)) * invdetDDL

    if ((gammaT < 0).or.&
        (gammaR < 0).or.&
        (gammaV < 0)) then
    print *, "    gammaT:", gammaT
    print *, "    gammaR:", gammaR
    print *, "    gammaV:", gammaV
    print *, "    lambda:", lambda
        exit
    end if
end do

!if (i < scalingIterations_max) then
if ((minDLgamma < gammaThreshold).and.&
    (minDLlambda < lambdaThreshold)) then

    momenta = momentaT*gammaT + &
              momentaR*gammaR + &
              momentaV*gammaV

    momentaIN = reshape(momenta,(/3*Natoms/))

end if

return

contains

    real(dp) function cofactor3of4(A44,i,j)
    use PARAMETERS
    real(dp),dimension(4,4),intent(in) :: A44
    real(dp),dimension(3,3) :: A33
    integer,intent(in) :: i, j
    integer :: m, n

    do m = 1, j-1
        do n = 1, i-1
            A33(m,n) = A44(m,n)
        end do

        do n = i+1, 4
            A33(m,n-1) = A44(m,n)
        end do
    end do

    do m = j+1, 4
        do n = 1, i-1
            A33(m-1,n) = A44(m,n)
        end do

        do n = i+1, 4
            A33(m-1,n-1) = A44(m,n)
        end do
    end do

    cofactor3of4 = &
       A33(1,1) * (A33(2,2)*A33(3,3)-A33(2,3)*A33(3,2)) - &
       A33(1,2) * (A33(2,1)*A33(3,3)-A33(2,3)*A33(3,1)) + &
       A33(1,3) * (A33(2,1)*A33(3,2)-A33(2,2)*A33(3,1))

    end function cofactor3of4

end subroutine scaleVelocities4

subroutine scaleVelocities5(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta
real(dp),dimension(3) :: coordsCM

real(dp),dimension(3,Natoms) :: &
        momentaT,momentaR,momentaV
real(dp),dimension(3,3,Natoms) :: momentaTRV

real(dp),dimension(3,3) :: E1,E2,E3
real(dp) :: linear_kcal_ratio = 1.0d0
real(dp) :: angular_kcal_ratio = 1.0d0

real(dp),dimension(3,3) :: L11, L112
real(dp),dimension(3) :: L1
real(dp) :: L0,L2

real(dp),dimension(3,3) :: L11L11T,L112L112T
real(dp),dimension(3,3) :: invL11L11T, Lintermediate, starL
real(dp) :: Lg
real(dp),dimension(1,1) :: LA, LB, LC
real(dp) :: LAA,LBB,LCC
real(dp) :: lambdadeterminant, lambdaplus, lambdaminus

real(dp),dimension(3,1) :: gammaplus,gammaminus
real(dp),dimension(3,1) :: DLgammaplus,DLgammaminus
real(dp),dimension(1,1) :: DLlambdaplus,DLlambdaminus

real(dp),dimension(3,1) :: Pdiff,Ldiff
real(dp),dimension(1,1) :: Tdiff

real(dp) :: gammaThreshold = 1.0d-5
real(dp) :: lambdaThreshold = 1.0d-5

integer :: scalingIterations_max = 100
integer :: i, j, m, n

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

call decomposeVelocities1(coords,&
             momenta,clusters,&
             momentaT,momentaR,momentaV)

do i = 1, Natoms
    momentaTRV(:,1,i) = momentaT(:,i)
    momentaTRV(:,2,i) = momentaR(:,i)
    momentaTRV(:,3,i) = momentaV(:,i)
end do

coordsCM = 0.0d0
do i = 1, Natoms
    coordsCM = coordsCM + masses(i)*coords(:,i)
end do
coordsCM = coordsCM / sum(masses(1:Natoms))

E1 = 0.0d0
E2 = 0.0d0
E3 = 0.0d0
do i = 1, Natoms
    coords(:,i) = coords(:,i) - coordsCM
    E1 = E1 + matmul(transpose(momentaTRV(:,:,i)),momentaTRV(:,:,i))/masses(i)
    E2 = E2 + momentaTRV(:,:,i)
    E3(1,:) = E3(1,:) + coords(2,i)*momentaTRV(3,:,i) &
                      - coords(3,i)*momentaTRV(2,:,i)
    E3(2,:) = E3(2,:) + coords(3,i)*momentaTRV(1,:,i) &
                      - coords(1,i)*momentaTRV(3,:,i)
    E3(3,:) = E3(3,:) + coords(1,i)*momentaTRV(2,:,i) &
                      - coords(2,i)*momentaTRV(1,:,i)
end do
E1 = E1 / C1
E2 = E2 * linear_kcal_ratio
E3 = E3 * angular_kcal_ratio

L11 = matmul(transpose(E2),E2) + matmul(transpose(E3),E3)
L1 = -2*reshape(matmul(reshape(L_baseline,(/ 1, 3 /)),E3),(/ 3 /))
L0 = sum(L_baseline**2)
L112 = E1
L2 = V - H_baseline

L11L11T = L11 + transpose(L11)
L112L112T = L112 + transpose(L112)

call get33inverse(L11L11T,invL11L11T)
Lintermediate = matmul(L112L112T,invL11L11T)
Lg = Lintermediate(1,1)+Lintermediate(2,2)+Lintermediate(3,3)

print *, ""
print *, "   L11L11T:"
print *, "        ", L11L11T(1,:)
print *, "        ", L11L11T(2,:)
print *, "        ", L11L11T(3,:)
print *, ""
print *, "invL11L11T:"
print *, "        ", invL11L11T(1,:)
print *, "        ", invL11L11T(2,:)
print *, "        ", invL11L11T(3,:)
print *, ""

starL = matmul(L11L11T,invL11L11T)

print *, "   product:"
print *, "        ", starL(1,:)
print *, "        ", starL(2,:)
print *, "        ", starL(3,:)
print *, ""

starL = matmul(invL11L11T,Lintermediate)
LA = matmul(matmul(matmul(matmul(&
            reshape(L1,(/ 1, 3 /)),&
            starL),&
            L112),&
            starL),&
            reshape(L1,(/ 3, 1 /)))
LB = matmul(matmul(matmul(matmul(&
            reshape(L1,(/ 1, 3 /)),&
            starL),&
            L112),&
            invL11L11T),&
            reshape(L1,(/ 3, 1 /))) * 2
!LB = matmul(matmul(matmul(matmul(&
!            reshape(L1,(/ 1, 3 /)),&
!            starL),&
!            L112),&
!            invL11L11T),&
!            reshape(L1,(/ 3, 1 /))) + &
!     matmul(matmul(matmul(matmul(&
!            reshape(L1,(/ 1, 3 /)),&
!            invL11L11T),&
!            L112),&
!            starL),&
!            reshape(L1,(/ 3, 1 /)))
LC = matmul(matmul(&
            reshape(L1,(/ 1, 3 /)),&
            starL),&
            reshape(L1,(/ 3, 1 /))) + L2

print *, "     starL:"
print *, "        ", starL(1,:)
print *, "        ", starL(2,:)
print *, "        ", starL(3,:)
print *, ""
print *, "  LA,LB,LC:"
print *, "        ", LA(1,:)
print *, "        ", LB(2,:)
print *, "        ", LC(3,:)
print *, ""
print *, "        Lg:"
print *, "        ", Lg
print *, ""

LAA = LA(1,1) + Lg*LB(1,1) + (Lg**2)*LC(1,1)
LBB = LB(1,1) + 2*Lg*LC(1,1)
LCC = LC(1,1)

lambdadeterminant = LBB**2 - 4 * LAA * LCC

if (lambdadeterminant < 0.0d0) then
    print *, ""
    print *, "Error in scaleVelocities5"
    print *, ""
    return
end if

lambdadeterminant = sqrt(lambdadeterminant)
lambdaplus = (-LBB + lambdadeterminant)/(2*LAA)
lambdaminus = (-LBB - lambdadeterminant)/(2*LAA)

gammaplus = matmul(invL11L11T + &
              (lambdaplus/(1.0d0+lambdaplus*Lg))*starL,&
              reshape(L1,(/ 3, 1 /)))
gammaminus = matmul(invL11L11T + &
              (lambdaminus/(1.0d0+lambdaminus*Lg))*starL,&
              reshape(L1,(/ 3, 1 /)))

DLgammaplus = matmul(L11L11T+lambdaplus*L112L112T,&
                     gammaplus) + reshape(L1,(/ 3, 1/))
DLlambdaplus = matmul(matmul(&
                      transpose(gammaplus),&
                      L112L112T),&
                      gammaplus) + L2
DLgammaminus = matmul(L11L11T+lambdaminus*L112L112T,&
                      gammaminus) + reshape(L1,(/ 3, 1/))
DLlambdaminus = matmul(matmul(&
                       transpose(gammaminus),&
                       L112L112T),&
                       gammaminus) + L2

print *, ""
print *, "Output of scaleVelocities5"
print *, ""

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "lambda plus:"
print *, "           ", lambdaplus
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "   DLlambda:"
print *, "           ", DLlambdaplus(1,1)
print *, "    DLgamma:"
print *, "           ", DLgammaplus(1,1)
print *, "           ", DLgammaplus(2,1)
print *, "           ", DLgammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

Pdiff = matmul(E2,gammaminus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaminus),&
                      E1),gammaminus) - (H_baseline - V)

print *, "lambda minus:"
print *, "           ", lambdaminus
print *, " gamma minus:"
print *, "           ", gammaminus(1,1)
print *, "           ", gammaminus(2,1)
print *, "           ", gammaminus(3,1)

print *, "   DLlambda:"
print *, "           ", DLlambdaminus(1,1)
print *, "    DLgamma:"
print *, "           ", DLgammaminus(1,1)
print *, "           ", DLgammaminus(2,1)
print *, "           ", DLgammaminus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

momenta = momentaT*gammaplus(1,1) + &
          momentaR*gammaplus(2,1) + &
          momentaV*gammaplus(3,1)

momentaIN = reshape(momenta,(/3*Natoms/))


starL = matmul(invL11L11T,Lintermediate)
print *, ""
print *, "   A(plus):"
print *, "        ", L11L11T(1,:) + lambdaplus*L112L112T(1,:)
print *, "        ", L11L11T(2,:) + lambdaplus*L112L112T(2,:)
print *, "        ", L11L11T(3,:) + lambdaplus*L112L112T(3,:)
print *, ""
print *, "invA(plus):"
print *, "        ", invL11L11T(1,:) + &
                     (lambdaplus/(1.0d0+lambdaplus*Lg))*starL(1,:)
print *, "        ", invL11L11T(2,:) + &
                     (lambdaplus/(1.0d0+lambdaplus*Lg))*starL(2,:)
print *, "        ", invL11L11T(3,:) + &
                     (lambdaplus/(1.0d0+lambdaplus*Lg))*starL(3,:)
print *, ""

starL = matmul(L11L11T+lambdaplus*L112L112T,&
               invL11L11T+(lambdaplus/(1.0d0+lambdaplus*Lg))*starL)

print *, "   product:"
print *, "        ", starL(1,:)
print *, "        ", starL(2,:)
print *, "        ", starL(3,:)
print *, ""


return

contains

    subroutine get33inverse(A33,invA33)
    use PARAMETERS
    real(dp),dimension(3,3),intent(in) :: A33
    real(dp),dimension(3,3),intent(out) :: invA33
    real(dp) :: detA33, invdetA33
    integer :: i, j, m1, m2, n1, n2

    do i = 1, 3
    do j = 1, 3

    m1 = i + 1
    if (m1 > 3) m1 = 1
    n1 = j + 1
    if (n1 > 3) n1 = 1

    m2 = m1 + 1
    if (m2 > 3) m2 = 1
    n2 = n1 + 1
    if (n2 > 3) n2 = 1

    invA33(i,j) = (&
        A33(m1,n1)*A33(m2,n2) - &
        A33(m1,n2)*A33(m2,n1)) !&
!           * ((-1)**(i+j))

    end do
    end do

    detA33 = &
        A33(1,1)*invA33(1,1) + &
        A33(1,2)*invA33(1,2) + &
        A33(1,3)*invA33(1,3)

    if (abs(detA33) < 1.0d-30) then
        invdetA33 = 0.0d0
        invA33 = 0.0d0
    else
        invdetA33 = 1.0d0 / detA33
        invA33 = transpose(invA33) * invdetA33

        print *, ""
        print *, "A(invA):"
        print *, matmul(A33,invA33)
    end if

    end subroutine get33inverse


end subroutine scaleVelocities5

subroutine scaleVelocities6(coordsIN,&
                   momentaIN,clusters,&
                   successfulScaling)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters
logical,optional,intent(out) :: successfulScaling

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta
real(dp),dimension(3) :: coordsCM

real(dp),dimension(3,Natoms) :: &
        momentaT,momentaR,momentaV
real(dp),dimension(3,3,Natoms) :: momentaTRV

real(dp),dimension(3,3) :: E1,E2,E3
real(dp) :: linear_kcal_ratio = 1.0d0
real(dp) :: angular_kcal_ratio = 1.0d0

real(dp),dimension(3,3) :: L11, L112, identity33
real(dp),dimension(3) :: L1
real(dp) :: L0,L2

real(dp),dimension(3,3) :: L11L11T,L112L112T
real(dp),dimension(3,3) :: invL11L11T, Lintermediate, starL
real(dp) :: Lg
real(dp),dimension(1,1) :: LA, LB, LC
real(dp) :: LAA,LBB,LCC
real(dp) :: lambdadeterminant, lambdaplus, lambdaminus, lambdaold

real(dp),dimension(3,1) :: gammaplus,gammaminus
real(dp),dimension(3,1) :: DLgammaplus,DLgammaminus
real(dp),dimension(1,1) :: DLlambdaplus,DLlambdaminus
real(dp),dimension(3,1) :: dummy1,dummy2

real(dp) :: fthreshold = 1.0d-6
real(dp),dimension(1,1) :: fplus,dfplus,fold
real(dp),dimension(1,1) :: fminus,dfminus
logical :: bisection_flag

real(dp),dimension(3,1) :: Pdiff,Ldiff
real(dp),dimension(1,1) :: Tdiff

real(dp) :: gammaThreshold = 1.0d-5
real(dp) :: lambdaThreshold = 1.0d-5

integer :: scalingIterations_max = 100
integer :: i, j, m, n

if (present(successfulScaling)) successfulScaling = .false.

! If there is too much potential energy then
! there is no way to scale the momenta to
! match the kinetic energy
if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

! Decompose the momenta into its
! translational,
! rotational, and
! vibrational components
call decomposeVelocities1(coords,&
             momenta,clusters,&
             momentaT,momentaR,momentaV)

! Load these momenta into a block
! matrix for ease of use later
do i = 1, Natoms
    momentaTRV(:,1,i) = momentaT(:,i)
    momentaTRV(:,2,i) = momentaR(:,i)
    momentaTRV(:,3,i) = momentaV(:,i)
end do

! Calculate the position of the
! center of mass for angular momentum
! calculations later
coordsCM = 0.0d0
do i = 1, Natoms
    coordsCM = coordsCM + masses(i)*coords(:,i)
end do
coordsCM = coordsCM / sum(masses(1:Natoms))

! These are as labelled in the paper
E1 = 0.0d0
E2 = 0.0d0
E3 = 0.0d0
do i = 1, Natoms
    coords(:,i) = coords(:,i) - coordsCM

    ! E1 is the kinetic energy operator
    E1 = E1 + matmul(transpose(momentaTRV(:,:,i)),momentaTRV(:,:,i))/masses(i)

    ! E2 is the center of mass momentum operator
    E2 = E2 + momentaTRV(:,:,i)

    ! E3 is the relative angular momentum operator
    E3(1,:) = E3(1,:) + coords(2,i)*momentaTRV(3,:,i) &
                      - coords(3,i)*momentaTRV(2,:,i)
    E3(2,:) = E3(2,:) + coords(3,i)*momentaTRV(1,:,i) &
                      - coords(1,i)*momentaTRV(3,:,i)
    E3(3,:) = E3(3,:) + coords(1,i)*momentaTRV(2,:,i) &
                      - coords(2,i)*momentaTRV(1,:,i)
end do
E1 = E1 / (2*C1)
E2 = E2 * linear_kcal_ratio
E3 = E3 * angular_kcal_ratio

identity33(1,:) = (/ 1.0d0, 0.0d0, 0.0d0 /)
identity33(2,:) = (/ 0.0d0, 1.0d0, 0.0d0 /)
identity33(3,:) = (/ 0.0d0, 0.0d0, 1.0d0 /)

! These are as labelled in the paper

! Regular
if (.false.) then
L11 = matmul(transpose(E2),E2) + matmul(transpose(E3),E3)
L1 = -2*reshape(matmul(reshape(L_baseline,(/ 1, 3 /)),E3),(/ 3 /))
L0 = sum(L_baseline**2)

! With regularization
else
L11 = matmul(transpose(E2),E2) + matmul(transpose(E3),E3) + &
      identity33
L1 = -2*(reshape(matmul(reshape(L_baseline,(/ 1, 3 /)),E3),(/ 3 /)) +&
         reshape((/ 1.0d0,1.0d0,1.0d0 /),(/ 3 /)))
L0 = sum(L_baseline**2) +&
          3.0d0
end if
L112 = E1
L2 = V - H_baseline

! These are the "tilde" operators in the
! paper (tilde L11 and tilde L112, respectively)
L11L11T = L11 + transpose(L11)
L112L112T = L112 + transpose(L112)

print *, ""
print *, "Output of scaleVelocities6"
print *, ""

! The iterations progress either with the
! Newton-Raphson method or with the bisection
! method;
! The bisection method is only used later in
! the iterative cycle after finding two
! values of opposite sign
bisection_flag = .false.

! The initial value of lambda is based on
! off the value of lambda associated with
! gamma = 1 (no scaling)
gammaplus = 1.0d0
lambdaplus = sqrt(sum((reshape(L1, (/ 3, 1/)) + &
                  matmul(L11L11T,gammaplus))**2) /&
                  sum(matmul(L112L112T,gammaplus)**2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "i=", 0
print *, ""
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

! This is for a test
dummy1 = matmul(L11L11T,gammaplus) + reshape(L1, (/ 3, 1 /))
dummy2 = matmul(L112L112T,gammaplus)

print *, ""
print *, "Absolute gamma values:"
print *, dummy1(1,1), dummy2(1,1)
print *, dummy1(2,1), dummy2(2,1)
print *, dummy1(3,1), dummy2(3,1)
print *, "Relative gamma values:"
print *, dummy1(1,1) / dummy2(1,1)
print *, dummy1(2,1) / dummy2(2,1)
print *, dummy1(3,1) / dummy2(3,1)
print *, ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The initial value is NOT unique and
! may have either a positive or negative
! value
call get33inverse(L11L11T + lambdaplus*L112L112T, starL)
if (sum((matmul(starL,reshape(L1,(/ 3, 1 /)))**2)) > 1.5) then
    lambdaplus = -lambdaplus
end if

! For now, only a maximum number of 100
! iterations are allowed
do i = 1, 100
print *, "i=",i
print *, ""

call get33inverse(L11L11T + lambdaplus*L112L112T, starL)
gammaplus = -matmul(starL,&
              reshape(L1,(/ 3, 1 /)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MISTAKE ?

!Lintermediate = matmul(matmul(&
!                starL,L112L112T),starL)
Lintermediate = matmul(matmul(&
                starL,L112),starL)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fplus = matmul(matmul(&
        reshape(L1,(/ 1, 3 /)),Lintermediate),&
              reshape(L1,(/ 3, 1 /))) + L2

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "lambda plus:"
print *, "           ", lambdaplus
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

if (abs(fplus(1,1)) < fthreshold) exit

!dfplus = matmul(matmul(matmul(matmul(&
!         reshape(L1,(/ 1, 3 /)),Lintermediate),&
!         L112),starL),&
!              reshape(L1,(/ 3, 1 /)))*2

if (bisection_flag) then

    if (fplus(1,1)*fminus(1,1) > 0) then
        fminus = fplus
        lambdaminus = lambdaplus
    else
        fold = fplus
        lambdaold = lambdaplus
    end if

    lambdaplus = (lambdaminus+lambdaold)/2

else

    if ((i > 1).and.(fplus(1,1)*fminus(1,1) < 0)) then
        bisection_flag= .true.

        lambdaold = lambdaminus
        fold = fminus
        lambdaminus = lambdaplus
        fminus = fplus
       
    else
        dfplus = &
            -matmul(matmul(matmul(matmul(&
             reshape(L1,(/ 1, 3 /)),Lintermediate),&
             L112),starL),&
                  reshape(L1,(/ 3, 1 /))) - &
             matmul(matmul(matmul(matmul(&
             reshape(L1,(/ 1, 3 /)),starL),&
             L112),Lintermediate),&
                  reshape(L1,(/ 3, 1 /)))
    
        print *, "     dfplus:"
        print *, "           ", dfplus(1,1)
    
    end if

    lambdaminus = lambdaplus
    fminus = fplus
    lambdaplus = lambdaplus - &
          fplus(1,1) / dfplus(1,1)
end if

print *, "      fplus:"
print *, "           ", fplus(1,1)

DLgammaplus = matmul(L11L11T+lambdaplus*L112L112T,&
                     gammaplus) + reshape(L1,(/ 3, 1/))
DLlambdaplus = matmul(matmul(&
                      transpose(gammaplus),&
                      L112L112T),&
                      gammaplus) + L2

print *, "   DLlambda:"
print *, "           ", DLlambdaplus(1,1)
print *, "    DLgamma:"
print *, "           ", DLgammaplus(1,1)
print *, "           ", DLgammaplus(2,1)
print *, "           ", DLgammaplus(3,1)

end do

!if (abs(fplus(1,1)) < fthreshold) then
if ((abs(Tdiff(1,1)) < fthreshold).and.&
    (all(gammaplus(1:3,1)>0))     ) then
    print *, "SUCCESSFUL scaling"
    momenta = momentaT*gammaplus(1,1) + &
              momentaR*gammaplus(2,1) + &
              momentaV*gammaplus(3,1)
    momentaIN = reshape(momenta,(/3*Natoms/))

    if (present(successfulScaling)) successfulScaling = .true.
else
    print *, " CANCELLED scaling"

call get33inverse(L11L11T + lambdaplus*L112L112T, starL)
gammaplus = 1.0d0
Lintermediate = matmul(matmul(&
                starL,L112L112T),starL)
fplus = matmul(matmul(&
        reshape(L1,(/ 1, 3 /)),Lintermediate),&
              reshape(L1,(/ 3, 1 /))) + L2

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "lambda plus:"
print *, "           ", lambdaplus
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

print *, "      fplus:"
print *, "           ", fplus(1,1)

DLgammaplus = matmul(L11L11T+lambdaplus*L112L112T,&
                     gammaplus) + reshape(L1,(/ 3, 1/))
DLlambdaplus = matmul(matmul(&
                      transpose(gammaplus),&
                      L112L112T),&
                      gammaplus) + L2

print *, "   DLlambda:"
print *, "           ", DLlambdaplus(1,1)
print *, "    DLgamma:"
print *, "           ", DLgammaplus(1,1)
print *, "           ", DLgammaplus(2,1)
print *, "           ", DLgammaplus(3,1)

end if

return

contains

    subroutine get33inverse(A33,invA33)
    use PARAMETERS
    real(dp),dimension(3,3),intent(in) :: A33
    real(dp),dimension(3,3),intent(out) :: invA33
    real(dp) :: detA33, invdetA33
    integer :: i, j, m1, m2, n1, n2

    do i = 1, 3
    do j = 1, 3

    m1 = i + 1
    if (m1 > 3) m1 = 1
    n1 = j + 1
    if (n1 > 3) n1 = 1

    m2 = m1 + 1
    if (m2 > 3) m2 = 1
    n2 = n1 + 1
    if (n2 > 3) n2 = 1

    invA33(i,j) = (&
        A33(m1,n1)*A33(m2,n2) - &
        A33(m1,n2)*A33(m2,n1)) !&
!           * ((-1)**(i+j))

    end do
    end do

    detA33 = &
        A33(1,1)*invA33(1,1) + &
        A33(1,2)*invA33(1,2) + &
        A33(1,3)*invA33(1,3)

    if (abs(detA33) < 1.0d-30) then
        invdetA33 = 0.0d0
        invA33 = 0.0d0
    else
        invdetA33 = 1.0d0 / detA33
        invA33 = transpose(invA33) * invdetA33

        print *, ""
        print *, "A(invA):"
        print *, matmul(A33,invA33)
    end if

    end subroutine get33inverse

end subroutine scaleVelocities6

subroutine scaleVelocities7(coordsIN,&
                   momentaIN,clusters)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta
real(dp),dimension(3,Natoms) :: momentaRV

real(dp),dimension(3) :: momentumTCM
real(dp) :: totalMass, Kprior, scalingT

integer :: i

if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

totalMass = 0.0d0
momentumTCM = 0.0d0
do i = 1, Natoms
    totalMass = totalMass + masses(i)
    momentumTCM = momentumTCM + momenta(:,i)
end do
momentumTCM = momentumTCM / totalMass

Kprior = 0.0d0
do i = 1, Natoms
    momentaRV(:,i) = momenta(:,i) - momentumTCM*masses(i)
    Kprior = Kprior + sum((momentaRV(:,i))**2)/masses(i)
end do
Kprior = Kprior / 2

scalingT = sqrt(C1*(H_baseline - V)/Kprior)

print *, "scaleVelocities14 Ouput:"
print *, "  scaling:", scalingT

momenta = momentaRV*scalingT
momentaIN = reshape(momenta,(/3*Natoms/))

return

end subroutine scaleVelocities7

subroutine scaleVelocities8(coordsIN,&
                   momentaIN,momentaError,&
                   successfulScaling)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
real(dp),dimension(3*Natoms),intent(in) :: momentaError
logical,optional,intent(out) :: successfulScaling

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta

real(dp),dimension(3) :: coordsCM
real(dp),dimension(3,Natoms) :: &
        momentaK,momentaG

real(dp),dimension(3,10,Natoms) :: momentaMatrix
real(dp),dimension(3,3,Natoms) :: coordsCrossMatrix

real(dp),dimension(1,10) :: PP, LP
real(dp),dimension(3,10) :: BP, RP
real(dp),dimension(10,10) :: BP2, RP2, MP, NP

real(dp) :: sing,cosg,sinb,cosb,sina,cosa
real(dp) :: sinasinb,sinacosb,cosasinb,cosacosb
real(dp) :: sinasing,sinacosg,cosasing,cosacosg
real(dp) :: sinbsing,sinbcosg,cosbsing,cosbcosg
real(dp) :: sinasinbsing,sinacosbsing,cosasinbsing,cosacosbsing
real(dp) :: sinasinbcosg,sinacosbcosg,cosasinbcosg,cosacosbcosg

real(dp),dimension(1,10) :: L1, L10
real(dp),dimension(10,10) :: L2, L20

real(dp),dimension(10,1) :: fd, fg, fb, fa
real(dp),dimension(10,1) :: fdg,fdb,fda,fgb,fga,fba
real(dp),dimension(10,1) :: fdd,fgg,fbb,faa
real(dp),dimension(10,1) :: fx

real(dp),dimension(5,5) :: ddL
real(dp),dimension(5,1) :: deltax, dL
real(dp),dimension(5) :: x

real(dp),dimension(3,1) :: Pdiff, Ldiff
real(dp),dimension(1,1) :: Tdiff

real(dp) :: Rfactor = 1.0d2

real(dp) :: linear_kcal_ratio = 1.0d0
real(dp) :: angular_kcal_ratio = 1.0d0

real(dp) :: TdiffThreshold = 1.0d-6
real(dp) :: dLThreshold = (1.0d-5)**2

integer :: scalingIterations_max = 100
integer :: i, j, m, n

if (present(successfulScaling)) successfulScaling = .false.

! If there is too much potential energy then
! there is no way to scale the momenta to
! match the kinetic energy
if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

! Decompose the momenta into the
! "known" (K) and "guess" (G)
! components
momentaG = reshape(momentaError,(/3,Natoms/))
momentaK = momenta - momentaG

! Load these momenta into a block
! matrix for ease of use later
momentaMatrix = 0.0d0
do i = 1, Natoms
    momentaMatrix(:,1,i) = momentaK(:,i)
    momentaMatrix(1,2:4,i) = momentaG(:,i)
    momentaMatrix(2,5:7,i) = momentaG(:,i)
    momentaMatrix(3,8:10,i) = momentaG(:,i)
end do

! Calculate the position of the
! center of mass for angular momentum
! calculations later
coordsCM = 0.0d0
do i = 1, Natoms
    coordsCM = coordsCM + masses(i)*coords(:,i)
end do
coordsCM = coordsCM / sum(masses(1:Natoms))

! Load these coordinates into a
! block matrix for ease of use later
coordsCrossMatrix = 0.0d0
do i = 1, Natoms
    coords(:,i) = coords(:,i) - coordsCM

    coordsCrossMatrix(1,2,i) = -coords(3,i)
    coordsCrossMatrix(1,3,i) =  coords(2,i)

    coordsCrossMatrix(2,1,i) =  coords(3,i)
    coordsCrossMatrix(2,3,i) = -coords(1,i)

    coordsCrossMatrix(3,1,i) = -coords(2,i)
    coordsCrossMatrix(3,2,i) =  coords(1,i)
end do

BP = 0.0d0
RP = 0.0d0
MP = 0.0d0
NP = 0.0d0
do i = 1, Natoms

    BP = BP + momentaMatrix(:,:,i)
    RP = RP + matmul(coordsCrossMatrix(:,:,i),&
                     momentaMatrix(:,:,i))
    MP = MP + matmul(transpose(momentaMatrix(:,:,i)),&
                     momentaMatrix(:,:,i)) / masses(i)
    NP = NP + matmul(transpose(momentaMatrix(:,:,i)),&
                     momentaMatrix(:,:,i))
end do

PP = 0.0d0
LP =-2*matmul(reshape(L_baseline,(/ 1, 3 /)),RP) ! BP)

MP = MP / (2*C1)
PP = PP * linear_kcal_ratio
LP = LP * angular_kcal_ratio

BP2 = matmul(transpose(BP),BP)
RP2 = matmul(transpose(RP),RP)



print *, ""
print *, "Output of scaleVelocities8"

! The normal method
if (.false.) then
    L10 = PP + LP
    L20 = 2*(BP2 + RP2)
    print *, "  regularization OFF"

! With regularization
else
    L10 = PP + LP - 2 * matmul(reshape(&
        (/ 1.0d0, 1.0d0, 0.0d0, 0.0d0, &
                  0.0d0, 1.0d0, 0.0d0, &
                  0.0d0, 0.0d0, 1.0d0 /),(/1,10/)),NP) * Rfactor
    L20 = 2*(BP2 + RP2 + NP*Rfactor)
    print *, "  regularization  ON ... Rfactor:", Rfactor

end if

print *, ""

write(6,FMT="(A36,3x,A36)") "      momentaK (known)      ",&
                            "      momentaG (guess)      "
do i = 1, Natoms
    write(6,FMT="(3(F11.7,1x),3x,3(F11.7,1x))") momentaK(:,i), momentaG(:,i)
end do

print *, ""

! The Newton-Raphson method is used

x(1) = 1.0d0
x(2) = 0.0d0
x(3) = 0.0d0
x(4) = 0.0d0
x(5) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! For now, only a maximum number of 100
! iterations are allowed
do i = 1, 100

print *, ""
print *, "i=",i
print *, ""

sing = sin(x(2))
cosg = cos(x(2))

sinb = sin(x(3))
cosb = cos(x(3))

sina = sin(x(4))
cosa = cos(x(4))

sinasinb = sina*sinb
sinacosb = sina*cosb
cosasinb = cosa*sinb
cosacosb = cosa*cosb

sinasing = sina*sing
sinacosg = sina*cosg
cosasing = cosa*sing
cosacosg = cosa*cosg

sinbsing = sinb*sing
sinbcosg = sinb*cosg
cosbsing = cosb*sing
cosbcosg = cosb*cosg

sinasinbsing = sinasinb*sing
sinacosbsing = sinacosb*sing
cosasinbsing = cosasinb*sing
cosacosbsing = cosacosb*sing

sinasinbcosg = sinasinb*cosg
sinacosbcosg = sinacosb*cosg
cosasinbcosg = cosasinb*cosg
cosacosbcosg = cosacosb*cosg

! Fd (derivative w.r.t capital gamma)

fd(1,1) = 0.0d0

fd(2,1) = cosacosb
fd(3,1) = cosasinbsing - sinacosg
fd(4,1) = cosasinbcosg + sinasing

fd(5,1) = sinacosb
fd(6,1) = sinasinbsing + cosacosg
fd(7,1) = sinasinbcosg - cosasing

fd(8,1) = -sinb
fd(9,1) = cosbsing
fd(10,1) = cosbcosg

! Fg (derivative w.r.t. lowercase gamma)

fg(1,1) = 0.0d0

fg(2,1) = 0.0d0
fg(3,1) = fd(4,1) * x(1) ! (cosasinbcosg + sinasing) * x(1)
fg(4,1) = -fd(3,1) * x(1) ! (-cosasinbsing + sinacosg) * x(1)

fg(5,1) = 0.0d0
fg(6,1) = fd(7,1) * x(1) ! (sinasinbcosg - cosasing) * x(1)
fg(7,1) = -fd(6,1) * x(1) ! (-sinasinbsing - cosacosg) * x(1)

fg(8,1) = 0.0d0
fg(9,1) = (cosbcosg) * x(1)
fg(10,1) = (-cosbsing) * x(1)

! Fb (derivative w.r.t. beta)

fb(1,1) = 0.0d0

fb(2,1) = (-cosasinb) * x(1)
fb(3,1) = (cosacosbsing) * x(1)
fb(4,1) = (cosacosbcosg) * x(1)

fb(5,1) = (-sinasinb) * x(1)
fb(6,1) = (sinacosbsing) * x(1)
fb(7,1) = (sinacosbcosg) * x(1)

fb(8,1) = (-cosb) * x(1)
fb(9,1) = (-sinbsing) * x(1)
fb(10,1) = (-sinbcosg) * x(1)

! Fa (derivative w.r.t. alpha)

fa(1,1) = 0.0d0

fa(2,1) = (-sinacosb) * x(1)
fa(3,1) = fg(7,1) ! (-sinasinbsing - cosacosg) * x(1)
fa(4,1) = -fg(6,1) ! (-sinasinbcosg + cosasing) * x(1)

fa(5,1) = (cosacosb) * x(1)
fa(6,1) = -fg(4,1) ! (cosasinbsing - sinacosg) * x(1)
fa(7,1) = fg(3,1) ! (cosasinbcosg + sinasing) * x(1)

fa(8,1) = 0.0d0
fa(9,1) = 0.0d0
fa(10,1) = 0.0d0

! Fdg (double derivative w.r.t. capital gamma and lowercase gamma)

fdg(1,1) = 0.0d0

fdg(2,1) = 0.0d0
fdg(3,1) = fd(4,1) ! (cosasinbcosg + sinasing)
fdg(4,1) = -fd(3,1) ! (-cosasinbsing + sinacosg)

fdg(5,1) = 0.0d0
fdg(6,1) = fd(7,1) ! (sinasinbcosg - cosasing)
fdg(7,1) = -fd(6,1)  ! (-sinasinbsing - cosacosg)

fdg(8,1) = 0.0d0
fdg(9,1) = (cosbcosg)
fdg(10,1) = (-cosbsing)

! Fdb (double derivative w.r.t. capital gamma and beta)

fdb(1,1) = 0.0d0

fdb(2,1) = (-cosasinb)
fdb(3,1) = (cosacosbsing)
fdb(4,1) = (cosacosbcosg)

fdb(5,1) = (-sinasinb)
fdb(6,1) = (sinacosbsing)
fdb(7,1) = (sinacosbcosg)

fdb(8,1) = (-cosb)
fdb(9,1) = (-sinbsing)
fdb(10,1) = (-sinbcosg)

! Fda (double derivative w.r.t. capital gamma and alpha)

fda(1,1) = 0.0d0

fda(2,1) = (-sinacosb)
fda(3,1) = fdg(7,1) ! (-sinasinbsing - cosacosg)
fda(4,1) = -fdg(6,1) ! (-sinasinbcosg + cosasing)

fda(5,1) = (cosacosb)
fda(6,1) = -fdg(4,1) ! (cosasinbsing - sinacosg)
fda(7,1) = fdg(3,1) ! (cosasinbcosg + sinasing)

fda(8,1) = 0.0d0
fda(9,1) = 0.0d0
fda(10,1) = 0.0d0

! Fgb (double derivative w.r.t. lowercase gamma and beta)

fgb(1,1) = 0.0d0

fgb(2,1) = 0.0d0
fgb(3,1) = fb(4,1) ! (cosacosbcosg) * x(1)
fgb(4,1) = -fb(3,1) ! (-cosacosbsing) * x(1)

fgb(5,1) = 0.0d0
fgb(6,1) = fb(7,1) ! (sinacosbcosg) * x(1)
fgb(7,1) = -fb(6,1) ! (-sinacosbsing) * x(1)

fgb(8,1) = 0.0d0
fgb(9,1) = fb(10,1) ! (-sinbcosg) * x(1)
fgb(10,1) = -fb(9,1) ! (sinbsing) * x(1)

! Fga (double derivative w.r.t. lowercase gamma and alpha)

fga(1,1) = 0.0d0

fga(2,1) = 0.0d0
fga(3,1) = fa(4,1) ! (-sinasinbcosg + cosasing) * x(1)
fga(4,1) = -fa(3,1) ! (sinasinbsing + cosacosg) * x(1)

fga(5,1) = 0.0d0
fga(6,1) = fa(7,1) ! (cosasinbcosg + sinasing) * x(1)
fga(7,1) = -fa(6,1) ! (-cosasinbsing + sinacosg) * x(1)

fga(8,1) = 0.0d0
fga(9,1) = 0.0d0
fga(10,1) = 0.0d0

! Fba (derivative w.r.t. beta and alpha)

fba(1,1) = 0.0d0

fba(2,1) = -fb(5,1) ! (sinasinb) * x(1)
fba(3,1) = -fb(6,1) ! (-sinacosbsing) * x(1)
fba(4,1) = -fb(7,1) ! (-sinacosbcosg) * x(1)

fba(5,1) = fb(2,1) ! (-cosasinb) * x(1)
fba(6,1) = fb(3,1) ! (cosacosbsing) * x(1)
fba(7,1) = fb(4,1) ! (cosacosbcosg) * x(1)

fba(8,1) = 0.0d0
fba(9,1) = 0.0d0
fba(10,1) = 0.0d0

! fx

fx(1,1) = 1.0d0

fx(2,1) = fa(5,1) ! (cosacosb) * x(1)
fx(3,1) = fa(6,1) ! (cosasinbsing - sinacosg) * x(1)
fx(4,1) = fa(7,1) ! (cosasinbcosg + sinasing)* x(1)

fx(5,1) = -fa(2,1) ! (sinacosb) * x(1)
fx(6,1) = -fa(3,1) ! (sinasinbsing + coscosg) * x(1)
fx(7,1) = -fa(4,1) ! (sinasinbcosg - cosasing) * x(1)

fx(8,1) = (-sinb) * x(1)
fx(9,1) = -fg(10,1) ! (cosbsing) * x(1)
fx(10,1) = fg(9,1) ! (cosbcosg) * x(1)

! Fdd (double derivative w.r.t. capital gamma)

fdd = 0.0d0

! Fgg (double derivative w.r.t. lowercase gamma)

fgg = -fx

fgg(1,1) = 0.0d0

fgg(2,1) = 0.0d0

fgg(5,1) = 0.0d0

fgg(8,1) = 0.0d0

! Fbb (double derivative w.r.t. beta)

fbb = -fx

fbb(1,1) = 0.0d0

fbb(3,1) = (-cosasinbsing) * x(1)
fbb(4,1) = (-cosasinbcosg) * x(1)

fbb(6,1) = (-sinasinbsing) * x(1)
fbb(7,1) = (-sinasinbcosg) * x(1)

! Faa (double derivative w.r.t. lowercase gamma)

faa = -fx

faa(1,1) = 0.0d0

faa(8,1) = 0.0d0

faa(9,1) = 0.0d0

faa(10,1) = 0.0d0




L1 = L10 ! PP + LP
L2 = L20 - 2 * x(5) * MP ! 2*(BP2 + RP2 - x(5) * MP)


dL(1:1,1:1) = matmul(matmul(transpose(fx),L2),fd) + matmul(L1,fd)
dL(2:2,1:1) = matmul(matmul(transpose(fx),L2),fg) + matmul(L1,fg)
dL(3:3,1:1) = matmul(matmul(transpose(fx),L2),fb) + matmul(L1,fb)
dL(4:4,1:1) = matmul(matmul(transpose(fx),L2),fa) + matmul(L1,fa)
dL(5:5,1:1) = H_baseline - V - matmul(matmul(transpose(fx),MP),fx)

ddL(1:1,1:1) = matmul(matmul(transpose(fx),L2),fdd) + &
               matmul(matmul(transpose(fd),L2),fd) + matmul(L1,fdd)
ddL(1:1,2:2) = matmul(matmul(transpose(fx),L2),fdg) + &
               matmul(matmul(transpose(fd),L2),fg) + matmul(L1,fdg)
ddL(1:1,3:3) = matmul(matmul(transpose(fx),L2),fdb) + &
               matmul(matmul(transpose(fd),L2),fb) + matmul(L1,fdb)
ddL(1:1,4:4) = matmul(matmul(transpose(fx),L2),fda) + &
               matmul(matmul(transpose(fd),L2),fa) + matmul(L1,fda)
ddL(1:1,5:5) = -2*matmul(matmul(transpose(fx),MP),fd)

ddL(2:2,1:1) = ddL(1:1,2:2)
ddL(2:2,2:2) = matmul(matmul(transpose(fx),L2),fgg) + &
               matmul(matmul(transpose(fg),L2),fg) + matmul(L1,fgg)
ddL(2:2,3:3) = matmul(matmul(transpose(fx),L2),fgb) + &
               matmul(matmul(transpose(fg),L2),fb) + matmul(L1,fgb)
ddL(2:2,4:4) = matmul(matmul(transpose(fx),L2),fga) + &
               matmul(matmul(transpose(fg),L2),fa) + matmul(L1,fga)
ddL(2:2,5:5) = -2*matmul(matmul(transpose(fx),MP),fg)

ddL(3:3,1:1) = ddL(1:1,3:3)
ddL(3:3,2:2) = ddL(2:2,3:3)
ddL(3:3,3:3) = matmul(matmul(transpose(fx),L2),fbb) + &
               matmul(matmul(transpose(fb),L2),fb) + matmul(L1,fbb)
ddL(3:3,4:4) = matmul(matmul(transpose(fx),L2),fba) + &
               matmul(matmul(transpose(fb),L2),fa) + matmul(L1,fba)
ddL(3:3,5:5) = -2*matmul(matmul(transpose(fx),MP),fb)

ddL(4:4,1:1) = ddL(1:1,4:4)
ddL(4:4,2:2) = ddL(2:2,4:4)
ddL(4:4,3:3) = ddL(3:3,4:4)
ddL(4:4,4:4) = matmul(matmul(transpose(fx),L2),faa) + &
               matmul(matmul(transpose(fa),L2),fa) + matmul(L1,faa)
ddL(4:4,5:5) = -2*matmul(matmul(transpose(fx),MP),fa)

ddL(5:5,1:1) = ddL(1:1,5:5)
ddL(5:5,2:2) = ddL(2:2,5:5)
ddL(5:5,3:3) = ddL(3:3,5:5)
ddL(5:5,4:4) = ddL(4:4,5:5)
ddL(5:5,5:5) = 0.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *, "      Gamma:"
print *, x(1)
print *, "      gamma:"
print *, x(2)
print *, "       beta:"
print *, x(3)
print *, "      alpha:"
print *, x(4)
print *, "     lambda:"
print *, x(5)

Pdiff = matmul(BP,fx)
Ldiff = matmul(RP,fx) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(transpose(fx),MP),fx) - (H_baseline - V)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

print *, "       |dL|:"
print *, sum(dL(1:5,1)**2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

print *, "        ddL:"
do j = 1, 5
!   write(6,FMT="(5x,5(E11.8,1x))") ddL(j,:)
    write(6,FMT=*) ddL(j,:)
end do

print *, "         dL:"
do j = 1, 5
!   write(6,FMT="(5x,E11.8,1x)") ddL(j,1)
    write(6,FMT=*) dL(j,1)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (sum(dL(1:5,1)**2) < dLthreshold) exit

call squareAxb(ddL,5,dL,deltax)

x = x - reshape(deltax,(/5/))

print *, "         dx:"
do j = 1, 5
!   write(6,FMT="(5x,E11.8,1x)") deltax(j,1)
    write(6,FMT=*) deltax(j,1)
end do

dL = matmul(ddL,deltax)

print *, "     ddL*dx:"
do j = 1, 5
!   write(6,FMT="(5x,E11.8,1x)") ddL(j,1)
    write(6,FMT=*) dL(j,1)
end do


end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!if ((sum(dL(1:5,1)**2) < dLthreshold).and.&
!if ((abs(Tdiff(1,1)) < TdiffThreshold).and.&
!    (x(1)>0)                     ) then
if (abs(Tdiff(1,1)) < TdiffThreshold) then
    print *, "SUCCESSFUL scaling"

    do i = 1, Natoms
        momenta(:,i) = reshape(matmul(momentaMatrix(:,:,i),fx),(/3/))
    end do
    momentaIN = reshape(momenta,(/3*Natoms/))

    if (present(successfulScaling)) successfulScaling = .true.
else
    print *, " CANCELLED scaling"

end if

return

end subroutine scaleVelocities8

subroutine scaleVelocities9(coordsIN,&
                   momentaIN,clusters,&
                   successfulScaling)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
COMMON/CONSTN/C1

integer :: NATOMS
double precision :: T, V, H
double precision :: C1

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(inout) :: momentaIN
integer,dimension(2,Natoms),intent(in) :: clusters
logical,optional,intent(out) :: successfulScaling

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta
real(dp),dimension(3) :: coordsCM

real(dp),dimension(3,Natoms) :: &
        momentaT,momentaR,momentaV
real(dp),dimension(3,3,Natoms) :: momentaTRV

real(dp),dimension(3,3) :: E1,E2,E3
real(dp) :: linear_kcal_ratio = 1.0d0
real(dp) :: angular_kcal_ratio = 1.0d0

real(dp),dimension(3,3) :: L11, L112, identity33
real(dp),dimension(3) :: L1
real(dp) :: L0,L2

real(dp),dimension(3,3) :: L11L11T,L112L112T
real(dp),dimension(3,3) :: invL11L11T, Lintermediate, starL
real(dp) :: Lg
real(dp),dimension(1,1) :: LA, LB, LC
real(dp) :: LAA,LBB,LCC
real(dp) :: lambdadeterminant, lambdaplus, lambdaminus, lambdaold

real(dp),dimension(3,1) :: gammaplus,gammaminus
real(dp),dimension(3,1) :: DLgammaplus,DLgammaminus
real(dp),dimension(1,1) :: DLlambdaplus,DLlambdaminus
real(dp),dimension(3,1) :: dummy1,dummy2

real(dp) :: fthreshold = 1.0d-6
real(dp),dimension(1,1) :: fplus,dfplus,fold
real(dp),dimension(1,1) :: fminus,dfminus
logical :: bisection_flag

real(dp),dimension(3,1) :: Pdiff,Ldiff
real(dp),dimension(1,1) :: Tdiff

real(dp) :: gammaThreshold = 1.0d-5
real(dp) :: lambdaThreshold = 1.0d-5

integer :: scalingIterations_max = 100
integer :: i, j, m, n

if (present(successfulScaling)) successfulScaling = .false.

! If there is too much potential energy then
! there is no way to scale the momenta to
! match the kinetic energy
if (H_baseline < V) return

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

! Decompose the momenta into its
! translational,
! rotational, and
! vibrational components
call decomposeVelocities3(coords,&
             momenta,clusters,&
             momentaT,momentaR,momentaV)

! Load these momenta into a block
! matrix for ease of use later
do i = 1, Natoms
    momentaTRV(:,1,i) = momentaT(:,i)
    momentaTRV(:,2,i) = momentaR(:,i)
    momentaTRV(:,3,i) = momentaV(:,i)
end do

! Calculate the position of the
! center of mass for angular momentum
! calculations later
coordsCM = 0.0d0
do i = 1, Natoms
    coordsCM = coordsCM + masses(i)*coords(:,i)
end do
coordsCM = coordsCM / sum(masses(1:Natoms))

! These are as labelled in the paper
E1 = 0.0d0
E2 = 0.0d0
E3 = 0.0d0
do i = 1, Natoms
    coords(:,i) = coords(:,i) - coordsCM

    ! E1 is the kinetic energy operator
    E1 = E1 + matmul(transpose(momentaTRV(:,:,i)),momentaTRV(:,:,i))/masses(i)

    ! E2 is the center of mass momentum operator
    E2 = E2 + momentaTRV(:,:,i)

    ! E3 is the relative angular momentum operator
    E3(1,:) = E3(1,:) + coords(2,i)*momentaTRV(3,:,i) &
                      - coords(3,i)*momentaTRV(2,:,i)
    E3(2,:) = E3(2,:) + coords(3,i)*momentaTRV(1,:,i) &
                      - coords(1,i)*momentaTRV(3,:,i)
    E3(3,:) = E3(3,:) + coords(1,i)*momentaTRV(2,:,i) &
                      - coords(2,i)*momentaTRV(1,:,i)
end do
E1 = E1 / (2*C1)
E2 = E2 * linear_kcal_ratio
E3 = E3 * angular_kcal_ratio

identity33(1,:) = (/ 1.0d0, 0.0d0, 0.0d0 /)
identity33(2,:) = (/ 0.0d0, 1.0d0, 0.0d0 /)
identity33(3,:) = (/ 0.0d0, 0.0d0, 1.0d0 /)

! These are as labelled in the paper

! Regular
if (.false.) then
L11 = matmul(transpose(E2),E2) + matmul(transpose(E3),E3)
L1 = -2*reshape(matmul(reshape(L_baseline,(/ 1, 3 /)),E3),(/ 3 /))
L0 = sum(L_baseline**2)

! With regularization
else
L11 = matmul(transpose(E2),E2) + matmul(transpose(E3),E3) + &
      identity33
L1 = -2*(reshape(matmul(reshape(L_baseline,(/ 1, 3 /)),E3),(/ 3 /)) +&
         reshape((/ 1.0d0,1.0d0,1.0d0 /),(/ 3 /)))
L0 = sum(L_baseline**2) +&
          3.0d0
end if
L112 = E1
L2 = V - H_baseline

! These are the "tilde" operators in the
! paper (tilde L11 and tilde L112, respectively)
L11L11T = L11 + transpose(L11)
L112L112T = L112 + transpose(L112)

print *, ""
print *, "Output of scaleVelocities9"
print *, ""

! The iterations progress either with the
! Newton-Raphson method or with the bisection
! method;
! The bisection method is only used later in
! the iterative cycle after finding two
! values of opposite sign
bisection_flag = .false.

! The initial value of lambda is based on
! off the value of lambda associated with
! gamma = 1 (no scaling)
gammaplus = 1.0d0
lambdaplus = sqrt(sum((reshape(L1, (/ 3, 1/)) + &
                  matmul(L11L11T,gammaplus))**2) /&
                  sum(matmul(L112L112T,gammaplus)**2))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "i=", 0
print *, ""
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

! This is for a test
dummy1 = matmul(L11L11T,gammaplus) + reshape(L1, (/ 3, 1 /))
dummy2 = matmul(L112L112T,gammaplus)

print *, ""
print *, "Absolute gamma values:"
print *, dummy1(1,1), dummy2(1,1)
print *, dummy1(2,1), dummy2(2,1)
print *, dummy1(3,1), dummy2(3,1)
print *, "Relative gamma values:"
print *, dummy1(1,1) / dummy2(1,1)
print *, dummy1(2,1) / dummy2(2,1)
print *, dummy1(3,1) / dummy2(3,1)
print *, ""

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The initial value is NOT unique and
! may have either a positive or negative
! value
call get33inverse(L11L11T + lambdaplus*L112L112T, starL)
if (sum((matmul(starL,reshape(L1,(/ 3, 1 /)))**2)) > 1.5) then
    lambdaplus = -lambdaplus
end if

! For now, only a maximum number of 100
! iterations are allowed
do i = 1, 100
print *, "i=",i
print *, ""

call get33inverse(L11L11T + lambdaplus*L112L112T, starL)
gammaplus = -matmul(starL,&
              reshape(L1,(/ 3, 1 /)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MISTAKE ?

!Lintermediate = matmul(matmul(&
!                starL,L112L112T),starL)
Lintermediate = matmul(matmul(&
                starL,L112),starL)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fplus = matmul(matmul(&
        reshape(L1,(/ 1, 3 /)),Lintermediate),&
              reshape(L1,(/ 3, 1 /))) + L2

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "lambda plus:"
print *, "           ", lambdaplus
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

if (abs(fplus(1,1)) < fthreshold) exit

!dfplus = matmul(matmul(matmul(matmul(&
!         reshape(L1,(/ 1, 3 /)),Lintermediate),&
!         L112),starL),&
!              reshape(L1,(/ 3, 1 /)))*2

if (bisection_flag) then

    if (fplus(1,1)*fminus(1,1) > 0) then
        fminus = fplus
        lambdaminus = lambdaplus
    else
        fold = fplus
        lambdaold = lambdaplus
    end if

    lambdaplus = (lambdaminus+lambdaold)/2

else

    if ((i > 1).and.(fplus(1,1)*fminus(1,1) < 0)) then
        bisection_flag= .true.

        lambdaold = lambdaminus
        fold = fminus
        lambdaminus = lambdaplus
        fminus = fplus
       
    else
        dfplus = &
            -matmul(matmul(matmul(matmul(&
             reshape(L1,(/ 1, 3 /)),Lintermediate),&
             L112),starL),&
                  reshape(L1,(/ 3, 1 /))) - &
             matmul(matmul(matmul(matmul(&
             reshape(L1,(/ 1, 3 /)),starL),&
             L112),Lintermediate),&
                  reshape(L1,(/ 3, 1 /)))
    
        print *, "     dfplus:"
        print *, "           ", dfplus(1,1)
    
    end if

    lambdaminus = lambdaplus
    fminus = fplus
    lambdaplus = lambdaplus - &
          fplus(1,1) / dfplus(1,1)
end if

print *, "      fplus:"
print *, "           ", fplus(1,1)

DLgammaplus = matmul(L11L11T+lambdaplus*L112L112T,&
                     gammaplus) + reshape(L1,(/ 3, 1/))
DLlambdaplus = matmul(matmul(&
                      transpose(gammaplus),&
                      L112L112T),&
                      gammaplus) + L2

print *, "   DLlambda:"
print *, "           ", DLlambdaplus(1,1)
print *, "    DLgamma:"
print *, "           ", DLgammaplus(1,1)
print *, "           ", DLgammaplus(2,1)
print *, "           ", DLgammaplus(3,1)

end do

!if (abs(fplus(1,1)) < fthreshold) then
if ((abs(Tdiff(1,1)) < fthreshold).and.&
    (all(gammaplus(1:3,1)>0))     ) then
    print *, "SUCCESSFUL scaling"
    momenta = momentaT*gammaplus(1,1) + &
              momentaR*gammaplus(2,1) + &
              momentaV*gammaplus(3,1)
    momentaIN = reshape(momenta,(/3*Natoms/))

    if (present(successfulScaling)) successfulScaling = .true.
else
    print *, " CANCELLED scaling"

call get33inverse(L11L11T + lambdaplus*L112L112T, starL)
gammaplus = 1.0d0
Lintermediate = matmul(matmul(&
                starL,L112L112T),starL)
fplus = matmul(matmul(&
        reshape(L1,(/ 1, 3 /)),Lintermediate),&
              reshape(L1,(/ 3, 1 /))) + L2

Pdiff = matmul(E2,gammaplus)
Ldiff = matmul(E3,gammaplus) - reshape(L_baseline, (/ 3, 1 /))
Tdiff = matmul(matmul(&
                      transpose(gammaplus),&
                      E1),gammaplus) - (H_baseline - V)

print *, "lambda plus:"
print *, "           ", lambdaplus
print *, " gamma plus:"
print *, "           ", gammaplus(1,1)
print *, "           ", gammaplus(2,1)
print *, "           ", gammaplus(3,1)

print *, "      Pdiff:"
print *, Pdiff
print *, "      Ldiff:"
print *, Ldiff
print *, "      Tdiff:"
print *, Tdiff

print *, "      fplus:"
print *, "           ", fplus(1,1)

DLgammaplus = matmul(L11L11T+lambdaplus*L112L112T,&
                     gammaplus) + reshape(L1,(/ 3, 1/))
DLlambdaplus = matmul(matmul(&
                      transpose(gammaplus),&
                      L112L112T),&
                      gammaplus) + L2

print *, "   DLlambda:"
print *, "           ", DLlambdaplus(1,1)
print *, "    DLgamma:"
print *, "           ", DLgammaplus(1,1)
print *, "           ", DLgammaplus(2,1)
print *, "           ", DLgammaplus(3,1)

end if

return

contains

    subroutine get33inverse(A33,invA33)
    use PARAMETERS
    real(dp),dimension(3,3),intent(in) :: A33
    real(dp),dimension(3,3),intent(out) :: invA33
    real(dp) :: detA33, invdetA33
    integer :: i, j, m1, m2, n1, n2

    do i = 1, 3
    do j = 1, 3

    m1 = i + 1
    if (m1 > 3) m1 = 1
    n1 = j + 1
    if (n1 > 3) n1 = 1

    m2 = m1 + 1
    if (m2 > 3) m2 = 1
    n2 = n1 + 1
    if (n2 > 3) n2 = 1

    invA33(i,j) = (&
        A33(m1,n1)*A33(m2,n2) - &
        A33(m1,n2)*A33(m2,n1)) !&
!           * ((-1)**(i+j))

    end do
    end do

    detA33 = &
        A33(1,1)*invA33(1,1) + &
        A33(1,2)*invA33(1,2) + &
        A33(1,3)*invA33(1,3)

    if (abs(detA33) < 1.0d-30) then
        invdetA33 = 0.0d0
        invA33 = 0.0d0
    else
        invdetA33 = 1.0d0 / detA33
        invA33 = transpose(invA33) * invdetA33

        print *, ""
        print *, "A(invA):"
        print *, matmul(A33,invA33)
    end if

    end subroutine get33inverse

end subroutine scaleVelocities9





subroutine getLWOclusters(coordsIN,momentaIN,Lout)
use PARAMETERS
use ANALYSIS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3*Natoms),intent(in) :: coordsIN
real(dp),dimension(3*Natoms),intent(in) :: momentaIN
real(dp),dimension(3),intent(out) :: Lout

real(dp),dimension(3,Natoms) :: coords
real(dp),dimension(3,Natoms) :: momenta
real(dp),dimension(3) :: coordsCM
integer :: i

coords = reshape(coordsIN,(/3,Natoms/))
momenta = reshape(momentaIN,(/3,Natoms/))

coordsCM = 0.0d0
do i = 1, Natoms
    coordsCM = coordsCM + masses(i)*coords(:,i)
end do
coordsCM = coordsCM / sum(masses(1:Natoms))

Lout = 0.0d0
do i = 1, Natoms
    coords(:,i) = coords(:,i) - coordsCM
    Lout(1) = Lout(1) + coords(2,i)*momenta(3,i) &
                      - coords(3,i)*momenta(2,i)
    Lout(2) = Lout(2) + coords(3,i)*momenta(1,i) &
                      - coords(1,i)*momenta(3,i)
    Lout(3) = Lout(3) + coords(1,i)*momenta(2,i) &
                      - coords(2,i)*momenta(1,i)
end do

return

end subroutine getLWOclusters







subroutine decomposeVelocities1(coords,&
                      momenta,clusters,&
             momentaT,momentaR,momentaV)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3,Natoms),intent(in) :: &
        coords, momenta
real(dp),dimension(3,Natoms),intent(out) :: &
        momentaT,momentaR,momentaV

integer,dimension(2,Natoms),intent(in) :: clusters

real(dp),dimension(3) :: deltacoords,&
        velocityCM1, CM1,&
        velocityCM2, CM2
real(dp) :: totalM1, totalM2

integer :: i, j

velocityCM1 = 0.0d0
CM1 = 0.0d0
totalM1 = 0.0d0
velocityCM2 = 0.0d0
CM2 = 0.0d0
totalM2 = 0.0d0
do j = 1, Natoms
    velocityCM1 = velocityCM1 + momenta(:,j)*clusters(1,j)
    CM1 = CM1 + masses(j)*coords(:,j)*clusters(1,j)
    totalM1 = totalM1 + masses(j)*clusters(1,j)
    velocityCM2 = velocityCM2 + momenta(:,j)*clusters(2,j)
    CM2 = CM2 + masses(j)*coords(:,j)*clusters(2,j)
    totalM2 = totalM2 + masses(j)*clusters(2,j)
end do

if (any(clusters(1,:) > 0)) then
    velocityCM1 = velocityCM1 / totalM1
    CM1 = CM1 / totalM1
end if
if (any(clusters(2,:) > 0)) then
    velocityCM2 = velocityCM2 / totalM2
    CM2 = CM2 / totalM2
end if

do j = 1, Natoms
    momentaT(:,j) = masses(j)*(velocityCM1*clusters(1,j)+&
                               velocityCM2*clusters(2,j))
end do

momentaR = momenta - momentaT

do j = 1, Natoms
    deltacoords = coords(:,j) - (CM1*clusters(1,j)+&
                                 CM2*clusters(2,j))
    momentaV(:,j) = dot_product(momentaR(:,j),deltacoords)*&
                        deltacoords / sum(deltacoords**2)
end do

momentaR = momentaR - momentaV

return

end subroutine decomposeVelocities1

subroutine decomposeVelocities2(coords,&
                      momenta,clusters,&
             momentaTCM,momentaTRe,momentaRV,&
             velocityCMRe,velocityCMRV,reducedmass)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3,Natoms),intent(in) :: &
        coords, momenta
integer,dimension(2,Natoms),intent(in) :: clusters
real(dp),dimension(3,Natoms),intent(out) :: &
        momentaTRe,momentaTCM,momentaRV
real(dp),dimension(3),intent(out) :: &
        velocityCMRe,velocityCMRV
real(dp),intent(out) :: reducedmass

real(dp),dimension(3) ::&
        momentumCM1,&
        momentumCM2
real(dp) :: totalM1, totalM2

integer :: i, j

momentumCM1 = 0.0d0
totalM1 = 0.0d0
momentumCM2 = 0.0d0
totalM2 = 0.0d0
do j = 1, Natoms
    momentumCM1 = momentumCM1 + momenta(:,j)*clusters(1,j)
    totalM1 = totalM1 + masses(j)*clusters(1,j)
    momentumCM2 = momentumCM2 + momenta(:,j)*clusters(2,j)
    totalM2 = totalM2 + masses(j)*clusters(2,j)
end do

do j = 1, Natoms
    momentaTCM(:,j) = (momentumCM1 + momentumCM2)/2
    momentaTRe(:,j) = momentumCM1*clusters(1,j) +&
                      momentumCM2*clusters(2,j) -&
                      momentaTCM(:,j)
end do

momentaRV = momenta - momentaTCM - momentaTRe

velocityCMRe = 0.0d0
velocityCMRV = 0.0d0
do j = 1, Natoms
    velocityCMRe = velocityCMRe + momentaTRe(:,j) /&
        (totalM1*clusters(1,j)-totalM2*clusters(2,j))
    velocityCMRV = velocityCMRV + &
        (momentaTCM(:,j) + momentaRV(:,j)) /&
        (totalM1*clusters(1,j)-totalM2*clusters(2,j))
end do

reducedmass = totalM1*totalM2/(totalM1+totalM2)

return

end subroutine decomposeVelocities2

subroutine decomposeVelocities2plus(coords,&
                      momenta,clusters,&
             momentaTCM,momentaTRe,momentaRV,&
             velocityTCM,velocityTRe,velocityCMRV,reducedmass)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3,Natoms),intent(in) :: &
        coords, momenta
integer,dimension(2,Natoms),intent(in) :: clusters
real(dp),dimension(3,Natoms),intent(out) :: &
        momentaTRe,momentaTCM,momentaRV
real(dp),dimension(3),intent(out) :: &
        velocityTCM,velocityTRe,velocityCMRV
real(dp),intent(out) :: reducedmass

real(dp),dimension(3) ::&
        momentumCM1,&
        momentumCM2
real(dp) :: totalM1, totalM2

integer :: i, j

momentumCM1 = 0.0d0
totalM1 = 0.0d0
momentumCM2 = 0.0d0
totalM2 = 0.0d0
do j = 1, Natoms
    momentumCM1 = momentumCM1 + momenta(:,j)*clusters(1,j)
    totalM1 = totalM1 + masses(j)*clusters(1,j)
    momentumCM2 = momentumCM2 + momenta(:,j)*clusters(2,j)
    totalM2 = totalM2 + masses(j)*clusters(2,j)
end do

do j = 1, Natoms
    momentaTCM(:,j) = (momentumCM1 + momentumCM2)/2
    momentaTRe(:,j) = momentumCM1*clusters(1,j) +&
                      momentumCM2*clusters(2,j) -&
                      momentaTCM(:,j)
end do

momentaRV = momenta - momentaTCM - momentaTRe

velocityTCM = 0.0d0
velocityTRe = 0.0d0
velocityCMRV = 0.0d0
do j = 1, Natoms
    velocityTCM = velocityTCM + momentaTCM(:,j) /&
        (totalM1*clusters(1,j)-totalM2*clusters(2,j))
    velocityTRe = velocityTRe + momentaTRe(:,j) /&
        (totalM1*clusters(1,j)-totalM2*clusters(2,j))
    velocityCMRV = velocityCMRV + momentaRV(:,j) /&
        (totalM1*clusters(1,j)-totalM2*clusters(2,j))
end do

reducedmass = totalM1*totalM2/(totalM1+totalM2)

return

end subroutine decomposeVelocities2plus

subroutine decomposeVelocities3(coords,&
                      momenta,clusters,&
             momentaT,momentaR,momentaV)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS

integer :: NATOMS

real(dp),dimension(3,Natoms),intent(in) :: &
        coords, momenta
real(dp),dimension(3,Natoms),intent(out) :: &
        momentaT,momentaR,momentaV

integer,dimension(2,Natoms),intent(in) :: clusters
integer :: clusterSum

real(dp),dimension(3,Natoms) :: deltacoords
real(dp),dimension(3) :: &
        velocityCM1, CM1,&
        velocityCM2, CM2
real(dp) :: totalM1, totalM2

real(dp) :: dum1
real(dp),dimension(3,1) :: angularMomentum
real(dp),dimension(3) :: &
        angularVelocity1, angularVelocity2
real(dp),dimension(3,3) :: &
        momentOfInertia,invMomentOfInertia

integer :: i, j

velocityCM1 = 0.0d0
CM1 = 0.0d0
totalM1 = 0.0d0
velocityCM2 = 0.0d0
CM2 = 0.0d0
totalM2 = 0.0d0
do j = 1, Natoms
    velocityCM1 = velocityCM1 + momenta(:,j)*clusters(1,j)
    CM1 = CM1 + masses(j)*coords(:,j)*clusters(1,j)
    totalM1 = totalM1 + masses(j)*clusters(1,j)
    velocityCM2 = velocityCM2 + momenta(:,j)*clusters(2,j)
    CM2 = CM2 + masses(j)*coords(:,j)*clusters(2,j)
    totalM2 = totalM2 + masses(j)*clusters(2,j)
end do

do j = 1, Natoms
    if (clusters(1,j) == 1) then
        momentaT(:,j) = masses(j)*velocityCM1
    else
        momentaT(:,j) = masses(j)*velocityCM2
    end if
    deltacoords(:,j) = coords(:,j) - (CM1*clusters(1,j)+&
                                      CM2*clusters(2,j))
end do

momentaV = momenta - momentaT

clusterSum = sum(clusters(1,:))
if (clusterSum > 0) then
    velocityCM1 = velocityCM1 / totalM1
    CM1 = CM1 / totalM1

    if (clusterSum > 1) then
    momentOfInertia = 0.0d0
    angularMomentum = 0.0d0
    do j = 1, Natoms
        if (clusters(1,j) == 0) cycle
        dum1 = sum(deltacoords(:,j)**2)
        momentOfInertia = momentOfInertia -&
             masses(j)*(matmul(reshape(deltacoords(:,j),(/3,1/)),&
                               reshape(deltacoords(:,j),(/1,3/))))
        do i = 1, 3
            momentOfInertia(i,i) = momentOfInertia(i,i) + dum1
        end do
        angularMomentum(1,1) = angularMomentum(1,1) +&
            deltacoords(2,j)*momentaV(3,j) -&
            deltacoords(3,j)*momentaV(2,j)
        angularMomentum(2,1) = angularMomentum(2,1) +&
            deltacoords(3,j)*momentaV(1,j) -&
            deltacoords(1,j)*momentaV(3,j)
        angularMomentum(3,1) = angularMomentum(3,1) +&
            deltacoords(1,j)*momentaV(2,j) -&
            deltacoords(2,j)*momentaV(1,j)
    end do

    call get33inverse(momentOfInertia,invMomentOfInertia)
    angularVelocity1 = reshape(matmul(&
             invMomentOfInertia,angularMomentum),(/3/))
    else
    angularVelocity1 = 0.0d0
    end if
end if

clusterSum = sum(clusters(2,:))
if (clusterSum > 0) then
    velocityCM2 = velocityCM2 / totalM2
    CM2 = CM2 / totalM2

    if (clusterSum > 1) then
    momentOfInertia = 0.0d0
    angularMomentum = 0.0d0
    do j = 1, Natoms
        if (clusters(2,j) == 0) cycle
        dum1 = sum(deltacoords(:,j)**2)
        momentOfInertia = momentOfInertia -&
             masses(j)*(matmul(reshape(deltacoords(:,j),(/3,1/)),&
                               reshape(deltacoords(:,j),(/1,3/))))
        do i = 1, 3
            momentOfInertia(i,i) = momentOfInertia(i,i) + dum1
        end do
        angularMomentum(1,1) = angularMomentum(1,1) +&
            deltacoords(2,j)*momentaV(3,j) -&
            deltacoords(3,j)*momentaV(2,j)
        angularMomentum(2,1) = angularMomentum(2,1) +&
            deltacoords(3,j)*momentaV(1,j) -&
            deltacoords(1,j)*momentaV(3,j)
        angularMomentum(3,1) = angularMomentum(3,1) +&
            deltacoords(1,j)*momentaV(2,j) -&
            deltacoords(2,j)*momentaV(1,j)
    end do

    call get33inverse(momentOfInertia,invMomentOfInertia)
    angularVelocity2 = reshape(matmul(&
             invMomentOfInertia,angularMomentum),(/3/))
    else
    angularVelocity2 = 0.0d0
    end if
end if


do j = 1, Natoms
    if (clusters(1,j) == 1) then
        momentaR(1,j) = masses(j)*&
                        angularVelocity1(2)*deltacoords(3,j)-&
                        angularVelocity1(3)*deltacoords(2,j)
        momentaR(2,j) = masses(j)*&
                        angularVelocity1(3)*deltacoords(1,j)-&
                        angularVelocity1(1)*deltacoords(3,j)
        momentaR(3,j) = masses(j)*&
                        angularVelocity1(1)*deltacoords(2,j)-&
                        angularVelocity1(2)*deltacoords(1,j)
    else
        momentaR(1,j) = masses(j)*&
                        angularVelocity2(2)*deltacoords(3,j)-&
                        angularVelocity2(3)*deltacoords(2,j)
        momentaR(2,j) = masses(j)*&
                        angularVelocity2(3)*deltacoords(1,j)-&
                        angularVelocity2(1)*deltacoords(3,j)
        momentaR(3,j) = masses(j)*&
                        angularVelocity2(1)*deltacoords(2,j)-&
                        angularVelocity2(2)*deltacoords(1,j)
    end if
end do

momentaV = momentaV - momentaR

return

contains

    subroutine get33inverse(A33,invA33)
    use PARAMETERS
    real(dp),dimension(3,3),intent(in) :: A33
    real(dp),dimension(3,3),intent(out) :: invA33
    real(dp) :: detA33, invdetA33
    integer :: i, j, m1, m2, n1, n2

    do i = 1, 3
    do j = 1, 3

    m1 = i + 1
    if (m1 > 3) m1 = 1
    n1 = j + 1
    if (n1 > 3) n1 = 1

    m2 = m1 + 1
    if (m2 > 3) m2 = 1
    n2 = n1 + 1
    if (n2 > 3) n2 = 1

    invA33(i,j) = (&
        A33(m1,n1)*A33(m2,n2) - &
        A33(m1,n2)*A33(m2,n1)) !&
!           * ((-1)**(i+j))

    end do
    end do

    detA33 = &
        A33(1,1)*invA33(1,1) + &
        A33(1,2)*invA33(1,2) + &
        A33(1,3)*invA33(1,3)

    if (abs(detA33) < 1.0d-30) then
        invdetA33 = 0.0d0
        invA33 = 0.0d0
    else
        invdetA33 = 1.0d0 / detA33
        invA33 = transpose(invA33) * invdetA33

        print *, ""
        print *, "A(invA):"
        print *, matmul(A33,invA33)
    end if

    end subroutine get33inverse



end subroutine decomposeVelocities3





!subroutine clusterAtoms(coords,velocities,&
subroutine clusterAtoms(coords,&
                        clusters1,clusters2)
use PARAMETERS
use FUNCTIONS
implicit none

COMMON/FORCES/NATOMS
COMMON/PRLIST/T,V,H
!COMMON/PQDOT/P(NDA3),QDOT(NDA3),W(NDA)

integer :: NATOMS
double precision :: T, V, H
!double precision :: W

real(dp),dimension(3,Natoms),intent(in) :: &
        coords
!       coords, velocities
integer,dimension(1,Natoms),intent(out) :: &
        clusters1
integer,dimension(2,Natoms),intent(out) :: &
        clusters2

integer,parameter :: K_max = 3
integer :: NK, K
integer,parameter :: CONV_max = 20
integer :: CONV

real(dp),dimension(3) :: mincoords, maxcoords

real(dp),dimension(K_max,3) :: centroid
real(dp),dimension(K_max,Natoms) :: dis
integer,dimension(K_max,Natoms) :: map
integer,dimension(K_max) :: dum

integer :: i,j

do i = 1, 3
    mincoords(i) = minval(coords(i,:),1)
    maxcoords(i) = maxval(coords(i,:),1)
end do

do NK = 1, K_max

    do K = 1, NK
        do i = 1, 3
            centroid(K,i) = mincoords(i) + &
                (K/(1.0*(NK+2)))*(maxcoords(i)-mincoords(i))
!               rand()*(maxcoords(i)-mincoords(i))
        end do
    end do

    dis = 0.0d0
    do K = 1, NK
        do j = 1, Natoms
           dis(K,j) = sum((coords(:,j)-centroid(K,:))**2)
        end do
    end do

    if (NK > 1) then
        map = 0
        do j = 1, Natoms
           map(minloc(dis(1:NK,j),1),j) = 1
        end do
    else
        map = 1
    end if

    do K = 1, NK
        dum(K) = sum(map(K,:))
    end do

    if (minval(dum(1:NK),1) < 0.5) cycle

    do CONV = 1, CONV_max

        centroid = 0.0d0
        do K = 1, NK
            do j = 1, Natoms
                do i = 1, 3
                    centroid(K,i) = centroid(K,i)+coords(i,j)*map(K,j)
                end do
            end do

            centroid(K,:) = centroid(K,:) / sum(map(K,:))
        end do

        do K = 1, NK
            do j = 1, Natoms
                dis(K,j) = sum((coords(:,j)-centroid(K,:))**2)
            end do
        end do

        map = 0
        do j = 1, Natoms
           map(minloc(dis(1:NK,j),1),j) = 1
        end do

    end do

    if (NK == 1) then
        clusters1 = map(1:NK,:)
    else if (NK == 2) then
        clusters2 = map(1:NK,:)
    else
    end if

end do

return

end subroutine clusterAtoms

end module velocityScaling
