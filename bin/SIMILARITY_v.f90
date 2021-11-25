!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       MODULE
!               SIMILARITY
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       PURPOSE
!               This module contains variables dictating the physics of the simulations
!               as well as what ensemble we are sampling initial conditions from and
!               some subroutines used to smooth this along
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       FILECHANNELS                    ACTION
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       SUBROUTINES                     ARGUMENTS               KIND
!
!               InitialSetup3                   coords                  intent(out),real(dp),dim(3,Natoms)
!                                               velocities              intent(out),real(dp),dim(3,Natoms)
!
!               BondedForce                     coords1                 intent(in),real(dp),dim(3)
!                                               coords2                 intent(in),real(dp),dim(3)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               r                       intent(in),real(dp),optional
!
!               NonBondedForce                  coords1                 intent(in),real(dp),dim(3)
!                                               coords2                 intent(in),real(dp),dim(3)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               gradient1               intent(out),real(dp),dim(3,Natoms)
!                                               r                       intent(in),real(dp),optional
!
!               cross                           A                       intent(in),real(dp),dim(3)
!                                               B                       intent(in),real(dp),dim(3)
!                                               AcrossB                 intent(out),real(dp),dim(3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!       CALLS                           MODULE
!
!               cross                           PHYSICS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module SIMILARITY
use DOUBLE
use PARAMETERS

integer,parameter :: NSIs = 3
integer,parameter :: Ntransform = 1

!real(dp),dimension(3,Natoms) :: coords_target
!real(dp),dimension(Natoms,Natoms) :: CM_target
real(dp),allocatable :: coords_target(:,:)
real(dp),allocatable :: CM_target(:,:)

!real(dp),dimension(3,Natoms,NSIs) :: coords_transformed
real(dp),allocatable :: coords_transformed(:,:,:)
real(dp),dimension(3,3,NSIs) :: U_transformed
real(dp),dimension(NSIs) :: SIs

real(dp),dimension(NSIs) :: default_SIs
real(dp),allocatable :: SIbuffer1(:,:)

real(dp),dimension(NSIs) :: min_SIs
real(dp),dimension(NSIs) :: largest_SIs
real(dp),dimension(NSIs) :: largest_weighted_SIs
real(dp),dimension(NSIs) :: largest_weighted_SIs2
real(dp),dimension(NSIs) :: interpolated_SIs


contains


subroutine getSIs(coords_in,coords_out,U_out,SIs_out)
    use PARAMETERS
    implicit none

    real,dimension(Natoms) :: ones

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(3,Natoms),intent(out) :: coords_out
    real(dp),dimension(3,3),intent(out) :: U_out
    real(dp),dimension(NSIs),intent(out) :: SIs_out

    ones = 1.0

    call getRMSD(coords_in,1)
    call getCMD(charges,coords_in,2)
    call getCMD(ones,coords_in,3)
!   call getCMD(ones,coords_in,2)
!   call getCMD(masses,coords_in,3)

    coords_out = coords_transformed(:,:,Ntransform)
    U_out = U_transformed(:,:,Ntransform)
    SIs_out = SIs

    return

end subroutine getSIs

subroutine setTarget(coords_in)
    use PARAMETERS
    use ANALYSIS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    integer :: i,j

    default_SIs(1) = 01.000100d0
    default_SIs(2) = 20.000100d0
    default_SIs(3) = 40.000100d0

    if (accept_worst) then
        do i = 1, NSIs
            SIbuffer1(i,:) = 0.0d0
        end do
    else
        do i = 1, NSIs
            SIbuffer1(i,:) = default_SIs(i)
        end do
    end if

    coords_target = coords_in
    call getCoulombMatrix(coords_in,CM_target)

    return

end subroutine setTarget


subroutine getCoulombMatrixJacobianProduct(coords_in,CMJP)
    use PARAMETERS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(3,Natoms),intent(out) :: CMJP
    real(dp),dimension(Natoms,Natoms) :: CM
    integer :: i, j, k, l
    CMJP = 0.0d0
    CM = 0.0d0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i < j) then
            cycle
        else if (i == j) then
            CM(i,j) = 1.0d0
        else
            CM(i,j) = 1.0d0 / sum((coords_in(:,i)-coords_in(:,j))**2)

            ! Now get the Jacobian product
            do k = 1, 3
                CMJP(k,i) = CMJP(k,i) + &
                    (coords_in(k,i)-coords_in(k,j)) * CM(i,j)
                CMJP(k,j) = CMJP(k,j) + &
                    (coords_in(k,j)-coords_in(k,i)) * CM(i,j)
            end do
        end if
    end do
    end do
    return
end subroutine getCoulombMatrixJacobianProduct

subroutine getCoulombMatrixJacobianProductForInputCLS(coords1,coords2,inputCLScol)
    use PARAMETERS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords1, coords2
    real(dp),dimension(3*Natoms),intent(out) :: inputCLScol
    real(dp),dimension(3,Natoms) :: CMJP
    call getCoulombMatrixJacobianProduct(coords2,CMJP)
    inputCLScol(1:3*Natoms) = reshape(CMJP,(/3*Natoms/))
    call getCoulombMatrixJacobianProduct(coords1,CMJP)
    inputCLScol(1:3*Natoms) = inputCLScol(1:3*Natoms) - reshape(CMJP,(/3*Natoms/))
    return
end subroutine getCoulombMatrixJacobianProductForInputCLS


subroutine getCoulombMatrix(coords_in,CM)
    use PARAMETERS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(Natoms,Natoms) :: CM

    integer :: i, j

    CM = 0.0d0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i < j) then
            cycle
        else if (i == j) then
            CM(i,j) = 1.0d0
        else
            CM(i,j) = 1.0d0 / &
                sqrt(sum((coords_in(:,i)-coords_in(:,j))**2))
        end if
    end do
    end do

    return

end subroutine getCoulombMatrix

subroutine getCoulombMatrixVector_test(coords_in,CMvec)
    use PARAMETERS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension((Natoms*(Natoms-1))/2),intent(out) :: CMvec

    integer :: index1
    integer :: i, j, k

    index1 = 0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i <= j) then
            cycle
        else
            index1 = index1 + 1
            CMvec(index1) = charges(i)*charges(j) / &
                sqrt(sum((coords_in(:,i)-coords_in(:,j))**2))
        end if
    end do
    end do

    return

end subroutine getCoulombMatrixVector_test

subroutine getCoulombMatrixVector(coords_in,CMvec)
    use PARAMETERS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension((Natoms*(Natoms-1))/2),intent(out) :: CMvec

    integer :: index1
    integer :: i, j, k

    index1 = 0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i <= j) then
            cycle
        else
            index1 = index1 + 1
            CMvec(index1) = 1.0d0 / &
                sqrt(sum((coords_in(:,i)-coords_in(:,j))**2))
        end if
    end do
    end do

    return

end subroutine getCoulombMatrixVector

subroutine getCoulombMatrixGradient(coords_in,CMgrad)
    use PARAMETERS
    implicit none

    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(Natoms,Natoms) :: CM
    real(dp),dimension((Natoms*(Natoms-1))/2,3*Natoms),intent(out) :: CMgrad

    integer :: index1
    integer :: index2, index3
    integer :: i, j, k

    CMgrad = 0.0d0

    index1 = 0
    CM = 0.0d0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i <= j) then
            cycle
        else
            index1 = index1 + 1
            CM(i,j) = 1.0d0 /&
                (sqrt(sum((coords_in(:,i)-coords_in(:,j))**2))**3)

            index2 = 3*(i-1)
            index3 = 3*(j-1)
            do k = 1, 3
                index2 = index2 + 1
                index3 = index3 + 1
                CMgrad(index1,index2) = (coords_in(k,i)-coords_in(k,j)) * CM(i,j)
                CMgrad(index1,index3) = - CMgrad(index1,index2)
            end do
        end if
    end do
    end do

    return

end subroutine getCoulombMatrixGradient

subroutine getCMD(coefficients,coords_in,NSI)
    use PARAMETERS
    implicit none

    integer,intent(in) :: NSI
    real,dimension(Natoms),intent(in) :: coefficients
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(Natoms,Natoms) :: CM
    real(dp) :: CMdiff
    integer :: i,j

    call getCoulombMatrix(coords_in,CM)

    CMdiff = 0.0d0
    do i = 1, Natoms
    do j = 1, Natoms
        if (i > j) then
            CMdiff = CMdiff + &
                     coefficients(i) * coefficients(j) * &
                     (CM(i,j) - CM_target(i,j))**2
        end if

    end do
    end do

    SIs(NSI) = sqrt(CMdiff/Natoms)
    coords_transformed(:,:,NSI) = coords_in
    U_transformed(:,:,NSI) = reshape((/ 1, 0, 0, &
                                        0, 1, 0, &
                                        0, 0, 1 /),(/ 3, 3 /))

    return

end subroutine getCMD

subroutine getRMSD(coords_in,NSI)
    use ls_rmsd_original
    use PARAMETERS
    implicit none

    integer,intent(in) :: NSI
    real(dp),dimension(3,Natoms),intent(in) :: coords_in
    real(dp),dimension(3) :: x_center,y_center
    integer :: i

    call rmsd_dp(Natoms,coords_in,coords_target,&
                 1,U_transformed(:,:,NSI),&
                 x_center,y_center,SIs(NSI))

    do i = 1, 3
        coords_transformed(i,:,NSI) =&
            coords_in(i,:) - x_center(i)
    end do

    coords_transformed(:,:,NSI) = matmul(&
        U_transformed(:,:,NSI),&
        coords_transformed(:,:,NSI))

    do i = 1, 3
        coords_transformed(i,:,NSI) =&
            coords_transformed(i,:,NSI) + y_center(i)
    end do

    return

end subroutine getRMSD





subroutine getFC(coords1,gradient1,&
                 coords2,gradient2)
use PARAMETERS
use ANALYSIS
implicit none


real(dp),dimension(3*Natoms),intent(in) :: coords1,gradient1
real(dp),dimension(3*Natoms),intent(in) :: coords2,gradient2

real(dp),dimension(9*Natoms*Natoms,3*Natoms) :: WP,WI
real(dp),dimension(3*Natoms) :: delta_coords
real(dp),dimension(3*Natoms) :: delta_gradient
real(dp) :: rmsd1, rmsd2
real(dp) :: FC_min = 1.0d-5

real(dp),dimension(NSIs) :: SIs_new
real(dp),dimension(3,3) :: U_new
real(dp),dimension(3,Natoms) :: coords_new,gradient_new

real(dp),dimension(9*Natoms*Natoms) :: FCvector
!real(dp),intent(out) :: FC
real(dp) :: FC1, FC2

integer :: i, j

print *, "FC coords in"
do i = 1, Natoms*3
    print *, coords1(i), coords2(i)
end do

delta_coords = coords1-coords2
rmsd1 = sqrt(sum(delta_coords**2)/Natoms)

delta_coords = delta_coords + sign(FC_min,delta_coords)
delta_coords = 1.0d0/delta_coords

delta_gradient = gradient1-gradient2

WP = 0.0d0
WI = 0.0d0
do i = 1, 3*Natoms
    WP((i-1)*3*Natoms+1:i*3*Natoms,i) = &
        delta_coords
    do j = 1, 3*Natoms
        WI((i-1)*3*Natoms+j,j) = &
            delta_coords(j)
    end do
end do

FCvector = matmul((WP-WI),delta_gradient)
FC1 = sqrt(sum(FCvector**2)/(Natoms*Natoms))

call setTarget(reshape(coords1,(/3,Natoms/)))
call getSIs(reshape(coords2,(/3,Natoms/)),&
            coords_new,U_new,SIs_new)

delta_coords = reshape(coords_new,(/3*Natoms/))-coords1
rmsd2 = sqrt(sum(delta_coords**2)/Natoms)

delta_coords = delta_coords + sign(FC_min,delta_coords)
delta_coords = 1.0d0/delta_coords

gradient_new = reshape(gradient2,(/3,Natoms/))
gradient_new = matmul(U_new,gradient_new)
delta_gradient = reshape(gradient_new,(/3*Natoms/))-gradient2

WP = 0.0d0
WI = 0.0d0
do i = 1, 3*Natoms
    WP((i-1)*3*Natoms+1:i*3*Natoms,i) = &
        delta_coords
    do j = 1, 3*Natoms
        WI((i-1)*3*Natoms+j,j) = &
            delta_coords(j)
    end do
end do

FCvector = matmul((WP-WI),delta_gradient)
FC2 = sqrt(sum(FCvector**2)/(Natoms*Natoms))



open(6666,file=gridpath5//"FC.dat",position="append")
write(6666,FMT="(4ES10.2)") rmsd1, FC1, rmsd2, FC2
close(6666)

return
end subroutine getFC

end module SIMILARITY
