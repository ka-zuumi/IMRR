program main
  use PARAMETERS
  use ANALYSIS
  use VARIABLES
  use SIMILARITY
  use interactMultipleGrids

  integer :: NC

  real(dp) :: NWChemConversion=0.04184D0*627.5095D0/0.529177249D0
  real(dp), dimension(3*Natoms) :: Q,P,PDOT
  real(dp),dimension(3,Natoms) :: approx_gradient
  real(dp) :: T, V, H

  integer :: NI = 3 * Natoms

  integer :: i, j, k

  call initializeGG()

  call initializeTrajectory(1)

  NC = 0

  do 

    !This reads VENUS/NWChem output
    call readNextQPandPDOT(Q,P,V,PDOT)
    PDOT = -PDOT * NWChemConversion

    !This reads an (extended) XYZ file
    !call readNextQandPDOTfromXYZ(Q,PDOT)
    !PDOT = -PDOT

    call getVarsHBrCO2(Q(1:3*NATOMS),5,vals,&
                 Nvar,BOND_LABELLING_DATA)

    if ((NC > 500).and.&
        (any(vals(1:Nvar) > var_maxvar))) STOP

    if (checkNoState_OnlyHistory_flag) then
      call checkNoState_OnlyHistory(&
            reshape(Q(1:NI),(/3,NI/3/)),&
            approx_gradient)
    else
      call checkState_PCM(vals,&
            reshape(Q(1:NI),(/3,NI/3/)),&
            approx_gradient,I,J,K)
    end if

    candidate_gradient = -candidate_gradient * NWChemConversion
    approx_gradient = -approx_gradient * NWChemConversion


    !We may add to the grid HERE at the half-step
    !(notice we are NOT adding the approx_gradient
    !but the PDOT from the previous iteration)

    if (grid_addition > 0) then
        call addState_new(vals,&
                Q(1:NI),&
                PDOT(1:NI),V)
    end if

    !For history testing, we make sure to skip the
    !first few steps

    if (NC < Nhistory_max) temp_reject_flag = .true.


    do i = 1, Nvar
      previous_var_index(i) =&
         int(vals(i)*divisor(i,Norder_order(1)+1))
    end do
    Vprev= V

    write(filechannel3,FMT="(2(F7.3,1x),I2,1x,8(ES10.3,1x))")&
        vals(1), vals(2),&
        Ninterpolation,&
        largest_weighted_SIs2(1),&
        largest_weighted_SIs(1),&
        SIbuffer1(1,1), interpolated_SIs(1),&
        SIbuffer1(2,1), interpolated_SIs(2),&
        sqrt(3*sum(&
                 (PDOT(1:NI)-reshape(candidate_gradient,(/ NI /))&
                 )**2)/NI) / NWChemConversion,&
        sqrt(3*sum(&
                 (PDOT(1:NI)-reshape(approx_gradient,(/ NI /))&
                 )**2)/NI) / NWChemConversion

    if (Nhistory_max > 0) then
      if (Nhistory < Nhistory_max) Nhistory = Nhistory + 1
      do k=1,ni
        if (Nhistory > 1) then
          qhistory(k,2:Nhistory)=qhistory(k,1:Nhistory-1)
          phistory(k,2:Nhistory)=phistory(k,1:Nhistory-1)
          pdothistory(k,2:Nhistory)=pdothistory(k,1:Nhistory-1)
        end if

        qhistory(k,1)=q(k)
        phistory(k,1)=p(k)
        pdothistory(k,1)=pdot(k)
      enddo
    end if

    NC = NC + 1
  end do

  call unsetAllocations()

  call finalizeGridCheck()


end program main

subroutine readNextQP(Q,P,V)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

double precision :: T, H

double precision,dimension(3*Natoms),intent(out) :: Q, P
double precision,intent(out) :: V
character(200) :: readtrajectory_aline
integer :: iostate
integer :: i,j

do
    read(filechannel6,iostat=iostate,&
         FMT="(A)") readtrajectory_aline
    if (iostate /= 0) exit
    if (readtrajectory_aline(1:5) == " XXXX") then
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline
        if (readtrajectory_aline(3:7) /= "THE C") exit
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline
        if (iostate /= 0) exit

        read(filechannel6,FMT=&
        "(18x,ES17.9,22x,ES17.9,A)")&
        T, V, readtrajectory_aline
        read(filechannel6,FMT=&
        "(18x,ES17.9,A)")&
        H, readtrajectory_aline
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline

        do j = 1, Natoms
            read(filechannel6,FMT=&
            "(3(F11.7,1x),3x,3(F11.7,1x),A)")&
            Q(j*3-2:j*3), P(j*3-2:j*3)
        end do

        exit
    end if
end do

return
end subroutine readNextQP

subroutine readNextQPandPDOT(Q,P,V,PDOT)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

double precision :: T, H

double precision,dimension(3*Natoms),intent(out) :: Q, P, PDOT
double precision,intent(out) :: V
double precision,dimension(3*Natoms) :: Q2

real(dp),parameter :: bohr2angstrom = 0.529177d0
character(200) :: readtrajectory_aline
integer :: i,j

logical :: QP_have_been_read_in
integer :: iostate

QP_have_been_read_in = .false.
do
    read(filechannel6,iostat=iostate,&
         FMT="(A)") readtrajectory_aline
    if (iostate /= 0) exit
    if (readtrajectory_aline(59:66) == "gradient") then
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline
        if (iostate /= 0) exit

        do j = 1, Natoms
            read(filechannel6,iostat=iostate,FMT=&
            "(10x,3(1x,F10.6),1x,3(1x,F10.6))")&
            Q2(j*3-2:j*3), PDOT(j*3-2:j*3)
            if (iostate /= 0) cycle
        end do

    elseif (readtrajectory_aline(1:5) == " XXXX") then
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline
        if (readtrajectory_aline(3:7) /= "THE C") exit
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline
        if (iostate /= 0) exit

        read(filechannel6,FMT=&
        "(18x,ES17.9,22x,ES17.9,A)")&
        T, V, readtrajectory_aline
        read(filechannel6,FMT=&
        "(18x,ES17.9,A)")&
        H, readtrajectory_aline
        read(filechannel6,iostat=iostate,&
             FMT="(A)") readtrajectory_aline

        do j = 1, Natoms
            read(filechannel6,FMT=&
            "(3(F11.7,1x),3x,3(F11.7,1x),A)")&
            Q(j*3-2:j*3), P(j*3-2:j*3)
!   print *, Q(j*3-2) / Q2(j*3-2)
!   print *, Q(j*3-1) / Q2(j*3-1)
!   print *, Q(j*3-0) / Q2(j*3-0)
        end do

        QP_have_been_read_in = .true.

        exit
    end if
end do

return
end subroutine readNextQPandPDOT

subroutine readNextQandPDOTfromXYZ(Q,PDOT)
use PARAMETERS
use FUNCTIONS
use ANALYSIS
use SIMILARITY
implicit none

double precision,dimension(3*Natoms),intent(out) :: Q, PDOT
character(200) :: readtrajectory_aline
integer :: i,j

integer :: iostate

read(filechannel6,FMT=*,iostat=iostate) i
if (iostate /= 0) then
  close(filechannel6)
  stop
end if

read(filechannel6,FMT=*)
do j = 1, Natoms
  read(filechannel6,FMT="(2x,A)") readtrajectory_aline
  read(readtrajectory_aline,FMT=*) &
            Q(j*3-2:j*3), PDOT(j*3-2:j*3)
end do

return
end subroutine readNextQandPDOTfromXYZ

