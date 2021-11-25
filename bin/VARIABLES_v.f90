!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               VARIABLES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This modules defines the collective variables used to
!               grid the library
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!               getVarsMaxMin                   coords                  intent(in),real(dp),dim(3,Natoms)
!                                               Natoms                  intent(in),integer
!                                               vals                    intent(out),real(dp),dim(Nvar)
!                                               Nvar                    intent(in),integer
!                                               labelling               intent(out),integer,dim(Natoms)
!
!               getDistanceSquared              coord1                  intent(in),real(dp),dim(3)
!                                               coord2                  intent(in),real(dp),dim(3)
!                                               d                       intent(out),real(dp)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               getDistanceSquared            VARIABLES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module VARIABLES
use DOUBLE
use PARAMETERS
implicit none
integer,dimension(4,2),parameter :: BOND_LABELS_TENTATIVE = reshape((/ 1, 2, 1, 2, &
                                                                       3, 4, 4, 3 /),    (/ 4, 2 /))
real(dp),dimension(Nvar) :: vals

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getVarsMaxMin
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine gets the collective variables of a frame
!
!               This subroutine assumes this is a TWO variable griding
!
!               Currently, there are two main segments, one for if this is a H-H2
!               system and one for if this is a H2-H2 system; to switch back and
!               forth, you must manually comment and uncomment the respective
!               segments
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               coords                          REAL(DP),DIM(3,Natoms)          The coordinates of the frame
!               Natoms                          INTEGER                         The number of atoms in the system
!               Nvar                            INTEGER                         The number of variables in the grid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               vals                            REAL(DP),DIM(Nvals)             The variables defining the frame
!               labelling                       INTEGER,DIM(Natoms)             The labeling scheme used to account for atom
!                                                                               indistinguishability
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!               length1                         REAL(DP)                        USED IN THE  H - H2 SYSTEM
!               length2                         REAL(DP)                        USED IN THE  H - H2 SYSTEM
!
!               lengths                         REAL(DP),DIM(4)                 USED IN THE H2 - H2 SYSTEM
!               min1_length                     INTEGER                         USED IN THE H2 - H2 SYSTEM
!               min2_length                     INTEGER                         USED IN THE H2 - H2 SYSTEM
!               max_length                      INTEGER                         USED IN THE H2 - H2 SYSTEM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine getVarsMaxMin(coords,Natoms,vals,Nvar,labelling)
implicit none
!
! for GNU fortran:
!integer,intent(in) :: Natoms, Nvar
! for intel fortran:
integer*8,intent(in) :: Natoms, Nvar
!
!

real(dp),dimension(3,Natoms),intent(in) :: coords
real(dp),dimension(Nvar),intent(out) :: vals
integer,dimension(Natoms),intent(inout) :: labelling
real(dp) :: length1,length2
!real(dp),dimension(4) :: lengths
!integer :: min1_length, min2_length, max_length
!integer :: i, start_label


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! H - H2  (test)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!call getDistanceSquared(coords(:,1),coords(:,2),length1)
!call getDistanceSquared(coords(:,1),coords(:,3),length2)
!
!if (length1 < length2) then
!	labelling(2) = 2
!	labelling(3) = 3
!	vals(1) = sqrt(length1)
!else
!	labelling(2) = 3
!	labelling(3) = 2
!	vals(1) = sqrt(length2)
!end if
!
!vals(2) = 1.0d0 + dot_product(coords(:,1)-coords(:,2),coords(:,1)-coords(:,3)) /&
!                  sqrt(length1*length2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! H - H2  (normal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!In this labeling scheme, we implicitly assume the nonbonded hydrogen
!is the 1st atom (as it is when I define it in PHYSICS)
!Otherwise, make sure you label it as the first atom

!Get the length between atoms 1 and 2
call getDistanceSquared(coords(:,3),coords(:,1),length1)

!Get the length between atoms 1 and 3
call getDistanceSquared(coords(:,3),coords(:,2),length2)

!If length1 < length2, then the 2nd atom must be closest to the
!nonbonded hydrogen, so it is labeled as 2nd
if (length1 < length2) then
	labelling(2) = 2
	labelling(3) = 3
	vals(1) = sqrt(length1)
	vals(2) = sqrt(length2)

!Otherwise, the 2nd atom is the farthest from the
!nonbonded hydrogen, so it is labeled as 3rd
else
	labelling(2) = 3
	labelling(3) = 2
	vals(1) = sqrt(length2)
	vals(2) = sqrt(length1)
end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! H2 - H2 (normal)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!First, get all lengths between nonbonded hydrogen
!!In PHYSICS, we define nonbonded hydrogen pairs by rows in the variable
!!BOND_LABELS_TENTATIVE; consequently, this variable must be 4x2 in dimension
!
!do i = 1, 4
!	call getDistanceSquared(coords(:,BOND_LABELS_TENTATIVE(i,1)),&
!				coords(:,BOND_LABELS_TENTATIVE(i,2)),lengths(i))
!end do
!
!!Next, locate the MIN, 2nd MIN, and MAX distance between nonbonded hydrogens
!!(Bear with me, this actually works)
!!Basically, because there are only 4 lengths, we can reduce this to
!!a psuedo-case statement that only checks whether the last two lengths
!!(between 1,4 and 3,2) are 1) lower, 2) in-between, or 3) greater
!!than the established MIN and MAX between the first two lengths
!
!min1_length = minloc(lengths(1:2),1)
!min2_length = 3 - min1_length
!max_length = min2_length
!do i = 3, 4
!	if (lengths(i) < lengths(min2_length)) then
!		if (lengths(i) < lengths(min1_length)) then
!			min2_length = min1_length
!			min1_length = i
!		else
!			min2_length = i
!		end if
!	else if (lengths(i) > lengths(max_length)) then
!		max_length = i
!	else
!	end if
!end do
!
!!Unfortunately, because of how we define the variable and grid, we
!!must take two square roots here
!
!vals(1) = sqrt(lengths(min1_length))
!vals(2) = sqrt(lengths(max_length))
!
!!Now, on to labeling
!!If the MIN AND MAX lengths both involve the same hydrogen then
!!we label that hydrogen as the first atom
!
!start_label = 0
!do i = 1, 4
!	if (any(BOND_LABELS_TENTATIVE(min1_length,:) == i, 1) .and. &
!	  (any(BOND_LABELS_TENTATIVE(max_length,:) == i, 1))) then
!		start_label = i
!	end if
!end do
!
!!If not, then we label the first hydrogen as the hydrogen involved
!!in both the MIN AND 2nd MIN lengths
!
!if (start_label == 0) then
!	do i = 1, 4
!		if (any(BOND_LABELS_TENTATIVE(min1_length,:) == i, 1) .and. &
!		  (any(BOND_LABELS_TENTATIVE(min2_length,:) == i, 1))) then
!			start_label = i
!		end if
!	end do
!end if
!
!!After deciding the first atom, we label the second atom
!!as the hydrogen not bonded to the first atom (using minimal length)
!
!labelling(1) = start_label
!if (BOND_LABELS_TENTATIVE(min1_length,1) == start_label) then
!	labelling(2) = BOND_LABELS_TENTATIVE(min1_length,2)
!else
!	labelling(2) = BOND_LABELS_TENTATIVE(min1_length,1)
!end if
!
!!Then we label the third atom to be the hydrogen not bonded
!!to the first atom (using maximal length)
!
!if (BOND_LABELS_TENTATIVE(max_length,1) == start_label) then
!	labelling(3) = BOND_LABELS_TENTATIVE(max_length,2)
!else
!	labelling(3) = BOND_LABELS_TENTATIVE(max_length,1)
!end if
!
!!And the fourth atom must be whichever index
!!we have not used yet, which is 4+3+2+1 - (x1+x2+x3)
!
!labelling(4) = 10 - sum(labelling(1:3))

end subroutine getVarsMaxMin


subroutine getVarsHBrCO2(coords,Natoms,vals,Nvar,labelling)
implicit none
integer,intent(in) :: Natoms, Nvar
real(dp),dimension(3,Natoms),intent(in) :: coords
real(dp),dimension(Nvar),intent(out) :: vals
integer,dimension(Natoms),intent(inout) :: labelling
real(dp) :: length1,length2
integer :: i

!Var1 will be the distance between the
!hydrogen (1) and the carbon (4)

call getDistanceSquared(&
    coords(:,1),coords(:,4),vals(1))

vals(1) = sqrt(vals(1))


!Var2 will be the minimum of the distances
!between the bromine (2) and the two
!oxygens (3,5)

call getDistanceSquared(&
    coords(:,2),coords(:,3),length1)
call getDistanceSquared(&
    coords(:,2),coords(:,5),length2)

vals(2) = sqrt(min(length1,length2))

!Who cares about labelling ... am I right?

do i = 1, Natoms
    labelling(i) = i
end do

end subroutine getVarsHBrCO2



!Variable one is the distance between the midpoints of two bonds (ONLY FOR H2 - H2)
subroutine getVar1(coords,Natoms,var1)
implicit none
integer,intent(in) :: Natoms
real(dp), dimension(3,Natoms),intent(in) :: coords
real(dp), intent(out) :: var1
real(dp) :: length1, length2, length3, length4

!call getDistanceSquared(coords(:,1)+coords(:,2),coords(:,3)+coords(:,4),var1)
!var1 = 0.5 * sqrt(var1)

!call getDistanceSquared(coords(:,1),coords(:,3),length1)
!call getDistanceSquared(coords(:,2),coords(:,4),length2)
!call getDistanceSquared(coords(:,1),coords(:,4),length3)
!call getDistanceSquared(coords(:,2),coords(:,3),length4)
!var1 = sqrt(min(length1,length2,length3,length4))

call getDistanceSquared(coords(:,1),coords(:,2),length1)
call getDistanceSquared(coords(:,1),coords(:,3),length2)
var1 = sqrt(min(length1,length2))

end subroutine getVar1

!Variable two is the cosine of the angle between the first H2 bond vector and
!						 the second H2 bond vector
!but shifted over by one, so that all values are positive
subroutine getVar2(coords,Natoms,var2)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
!real(dp), dimension(3) :: bond_vector1, bond_vector2, bond_CM_vector
!real(dp) :: bond1_length, bond2_length, bond_CM_length
!real(dp), dimension(4) :: bond_lengths       !NOTE: This is not generic yet
!integer :: bond_pair
real(dp) :: length1, length2, length3, length4
real(dp), intent(out) :: var2

call getDistanceSquared(coords(:,1),coords(:,2),length1)
call getDistanceSquared(coords(:,1),coords(:,3),length2)
var2 = sqrt(max(length1,length2))

!call getDistanceSquared(coords(:,1),coords(:,3),length1)
!call getDistanceSquared(coords(:,2),coords(:,4),length2)
!call getDistanceSquared(coords(:,1),coords(:,4),length3)
!call getDistanceSquared(coords(:,2),coords(:,3),length4)
!var2 = sqrt(max(length1,length2,length3,length4))

!call getDistanceSquared(coords(:,1),coords(:,3),bond_lengths(1))
!call getDistanceSquared(coords(:,2),coords(:,4),bond_lengths(2))
!call getDistanceSquared(coords(:,1),coords(:,4),bond_lengths(3))
!call getDistanceSquared(coords(:,2),coords(:,3),bond_lengths(4))
!bond_pair = minloc(bond_lengths,dim=1)
!
!if (bond_pair == 1) then
!	bond_vector1 = coords(:,2) - coords(:,1)
!	bond_vector2 = coords(:,4) - coords(:,3)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!else if (bond_pair == 2) then
!	bond_vector1 = coords(:,1) - coords(:,2)
!	bond_vector2 = coords(:,3) - coords(:,4)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!else if (bond_pair == 3) then
!	bond_vector1 = coords(:,2) - coords(:,1)
!	bond_vector2 = coords(:,3) - coords(:,4)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!else
!	bond_vector1 = coords(:,1) - coords(:,2)
!	bond_vector2 = coords(:,4) - coords(:,3)
!	bond1_length = sum(bond_vector1**2)
!	bond2_length = sum(bond_vector2**2)
!
!end if
!
!var2 = 1.0d0 + dot_product(bond_vector1,bond_vector2)/sqrt(bond1_length*bond2_length)

end subroutine getVar2

!Variables three and four are the distance between atoms 1 and 2 and 1 and 3, respectively
subroutine getVar3(coords,Natoms,var3)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var3

call getDistanceSquared(coords(:,1),coords(:,2),var3)
var3 = sqrt(var3)
end subroutine getVar3

subroutine getVar4(coords,Natoms,var4)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var4

call getDistanceSquared(coords(:,1),coords(:,3),var4)
var4 = sqrt(var4)
end subroutine getVar4

!Variables five and six are the distance between atoms 1 and 3 and 2 and 4, respectively
subroutine getVar5(coords,Natoms,var5)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var5

call getDistanceSquared(coords(:,1),coords(:,3),var5)
var5 = sqrt(var5)
end subroutine getVar5

subroutine getVar6(coords,Natoms,var6)
implicit none
integer, intent(in) :: Natoms
real(dp), dimension(3,Natoms), intent(in) :: coords
real(dp), intent(out) :: var6

call getDistanceSquared(coords(:,2),coords(:,4),var6)
var6 = sqrt(var6)
end subroutine getVar6



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINE
!               getDistanceSquared
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This subroutine gets the distance (squared) between two atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       INPUT                           KIND                            DESCRIPTION
!
!               coord1                          REAL(DP),DIM(3)                 The coordinates of the first atom
!               coord2                          REAL(DP),DIM(3)                 The coordinates of the second atom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       OUTPUT                          KIND                            DESCRIPTION
!
!               d                               REAL(DP)                        The distance between the two atoms
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       IMPORTANT VARIABLES             KIND                            DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILES                             FILETYPE                      DESCRIPTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getDistanceSquared(coord1,coord2,d)
implicit none
real(dp), dimension(3), intent(in) :: coord1, coord2
real(dp), intent(out) :: d

d = (coord1(1)-coord2(1))**2 + (coord1(2)-coord2(2))**2 + &
        (coord1(3)-coord2(3))**2

end subroutine getDistanceSquared

end module VARIABLES
