module FUNCTIONS
use DOUBLE
implicit none

!real(dp), dimension(10,1000) :: RKHS_A
!real(dp), dimension(10,12000) :: RKHS_dA
!real(dp), dimension(12000) :: RKHS_F
!real(dp),dimension(12000,12000) :: RKHS_Qt, RKHS_L
real(dp), allocatable :: RKHS_A(:,:)
real(dp), allocatable :: RKHS_dA(:,:)
real(dp), allocatable :: RKHS_F(:)
real(dp), allocatable :: RKHS_Qt(:,:), RKHS_L(:,:)
logical :: RKHSarrays_are_not_allocated = .true.

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      PARTITION FUNCTION (QuickSort Alogrithm by Tony Hoare)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   matrix A, dim (rows,cols)       "values to be sorted"
!               matrix B, dim (rows,1)          "pointers to orignal indexes"
!               integer start_index             "in case of sorting only a part of the martix A, 
!                                                the starting index (before sorting) of the to-be-sorted part"
!               integer end_index               "in case of sorting only a part of the martix A, 
!                                                the final index (before sorting) of the to-be-sorted part"
!               integer var                     "the index of the column of the martix A that is used as sorting criteria"
!      OUTPUT:  integer pivot_index             "the index of the row that partitions matrix A according to var"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Column "var" of matrix A is partitioned
!       but only for the sub-matrix A (start_index:end_index)
!       Any position changes in A are mimicked in B.
!       The value of end_index is moved to pivot_index so that
!       all values of A(start:pivot) are <= A(pivot) and
!       all values of A(pivot:end) are >= A(pivot).
!       Pivot_index is then returned
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine partition(A,B,rows,cols,start_index,end_index,var,pivot_index)
implicit none
integer, intent(in) :: start_index,end_index,rows,cols,var
integer, intent(out) :: pivot_index
real(dp), dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B
integer :: i, j
integer :: test_value,pivot_value

pivot_value = A(end_index,var)
j = start_index-1
do i = start_index, end_index-1
        test_value = A(i,var)
        if (test_value.le.pivot_value) then
                j = j + 1
                call swapD(A,rows,cols,i,j)
                call swapI(B,rows,1,i,j)
        end if
end do
pivot_index = j+1
call swapD(A,rows,cols,pivot_index,end_index)
call swapI(B,rows,1,pivot_index,end_index)

end subroutine partition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      SWAP FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IN/OUT:   matrix A, dim (rows,cols)       "values"
!      INPUT:    integer i1                      "spot1"
!                integer i2                      "spot2"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Spots 1 and 2 are swapped in matrix A
!       
!       swapD is for matrices of type real(dp) (vals)
!       swapR is for matrices of type real (vals)
!       swapI is for matrices of type integer (indexer)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine swapD(A,rows,cols,i1, i2)
implicit none
integer, intent(in) :: rows, cols, i1, i2
real(dp), dimension(rows,cols), intent(out) :: A
integer :: i
real(dp), dimension(cols) :: val1, val2

do i = 1, cols
        val1(i) = A(i1,i)
        val2(i) = A(i2,i)
end do
do i = 1, cols
        A(i2,i) = val1(i)
        A(i1,i) = val2(i)
end do

end subroutine swapD

subroutine swapR(A,rows,cols,i1, i2)
implicit none
integer, intent(in) :: rows, cols, i1, i2
real, dimension(rows,cols), intent(out) :: A
integer :: i
real, dimension(cols) :: val1, val2

do i = 1, cols
        val1(i) = A(i1,i)
        val2(i) = A(i2,i)
end do
do i = 1, cols
        A(i2,i) = val1(i)
        A(i1,i) = val2(i)
end do

end subroutine swapR

subroutine swapI(A,rows,cols,i1, i2)
implicit none
integer, intent(in) :: rows, cols, i1, i2
integer, dimension(rows,cols), intent(out) :: A
integer :: i
integer, dimension(cols) :: val1, val2

do i = 1, cols
        val1(i) = A(i1,i)
        val2(i) = A(i2,i)
end do
do i = 1, cols
        A(i2,i) = val1(i)
        A(i1,i) = val2(i)
end do

end subroutine swapI


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      QSORT2 FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      IN/OUT:  matrix A, dim (rows,cols)       "values to be sorted"
!               matrix B, dim (rows,1)          "pointers to orignal indexes"
!  The following inputs control the sorting
!      INPUT:   integer start_index             "in case of sorting only a part of the martix A, 
!                                                the starting index (before sorting) of the to-be-sorted part"
!               integer end_index               "in case of sorting only a part of the martix A, 
!                                                the final index (before sorting) of the to-be-sorted part"
!               integer var                     "the index of the column of the martix A that is used as sorting criteria"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Column "var" of matrix A is sorted but only for the sub-matrix A (start_index:end_index).
!       Any position changes in A are mimicked in B.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  ex.  A = [ 1.45  2.39 /       B = [ 1 /
!             0.98  7.45 /             2 /
!             2.11  3.00   ]           3   ]
!
!       QSORT2(A,B,3,2,   1,3,1) produces:
!       A = [ 0.98  7.45 /       B = [ 2 /
!             1.45  2.39 /             1 /
!             2.11  3.00   ]           3   ]
!
!       QSORT(A,B,3,2,    1,3,2) produces:
!       A = [ 1.45  2.39 /       B = [ 1 /
!             2.11  3.00 /             3 /
!             0.98  7.45   ]           2   ]
!
!       QSORT2(A,B,3,2,   1,2,1) produces:
!       A = [ 0.98  7.45 /       B = [ 2 /
!             1.45  2.39 /             1 /
!             2.11  3.00   ]           3   ]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine qsort2(A,B,rows,cols,start_index,end_index,var)
implicit none
integer, intent(in) :: rows, cols, start_index, end_index, var
integer :: pivot_index
real(dp), dimension(rows,cols), intent(out) :: A
integer, dimension(rows,1), intent(out) :: B

call partition(A,B,rows,cols,start_index,end_index,var,pivot_index)
if (start_index /= pivot_index) &
                call qsort2(A,B,rows,cols,start_index,pivot_index-1,var)
if (end_index /= pivot_index) &
                call qsort2(A,B,rows,cols,pivot_index+1,end_index,var)

end subroutine qsort2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      GRIDER FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      INPUT:   matrix A, dim (rows,cols)       "already-sorted values to be grided"
!               real gridline_spacing           "length of grid spacing"
!               real gridline_start             "value of first gridline"
!               integer max_gridlines           "the number of gridlines"
!               integer start_index             "start"
!               integer end_index               "end"
!               integer var                     "which column"
!      OUTPUT:  array grid, dim(Ngridmax)       "indexes of A for gridlines"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Column var of matrix A is assumed to be sorted
!       This functions grids submatrix A(start:end) into grid
!       There is a plain test description in addCells4 as well
!       
!       As a visual example, let's grid
!               A = [1.1, 2.7, 4.1, 4.2, 8.9, 9.6]
!
!       grider(grid1,A, 1.0, 0.0, 10, 6,1, 1, 6, 1)
!       grider(grid2,A, 2.0, 0.0, 10, 6,1, 1, 6, 1)
!       grider(grid3,A, 1.0, 2.0, 10, 6,1, 1, 6, 1)
!       grider(grid4,A, 1.0, 0.0,  5, 6,1, 1, 6, 1)
!       grider(grid5,A, 1.0, 0.0, 10, 6,1, 4, 6, 1)
!
!           A = | | 1.1 | 2.7 | | 4.1 4.2 | | | | 8.9 | 9.6 |
!       grid1 = 1,1,    2,    3,3,        5,5,5,5     6    
!
!           A = | 1.1 | 2.7 | 4.1 4.2 | | 8.9 9.6 | | | | | |
!       grid2 = 1,    2,    3,        5,5,        7,7,7,7,7
!
!           A = 1.1 | 2.7 | | 4.1 4.2 | | | | 8.9 | 9.6 | | |
!       grid3 =     2,    3,3,        5,5,5,5,    6,    7,7
!
!           A = | | 1.1 | 2.6 | | 4.1 4.2 |
!       grid4 = 1,1,    2,    3,3,        
!
!           A = | | 1.1 | 2.7 | | 4.1 4.2 | | | | 8.9 | 9.6 |
!       grid1 = 4,4,    4,    4,4,        5,5,5,5     6    
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine grider(grid,A,gridline_spacing,gridline_start,max_gridlines,&
                  rows,cols,start_index,end_index,var)
implicit none
integer, intent(in) :: rows, cols, var, max_gridlines
integer, intent(in) :: start_index, end_index
integer :: A_index, grid_index, i
real, intent(in) :: gridline_spacing, gridline_start
real, dimension(rows,cols), intent(in) :: A
real :: gridline
integer, dimension(max_gridlines), intent(out) :: grid

A_index = start_index
grid_index = 1
gridline = gridline_start + gridline_spacing
do
! RS: Question
! RS: gridline is the lower boundary of the grid 
! RS: b.c. grid_index-1 = 0 at the beginning
! RS: Is this correct? Shouldn't you use the upper boundary?
! RS: the first round until (A_index == end_index) is pointless
! RS: Since A is already sorted, I think we can use a much algorithm to grid
! RS: Like the "finding the root" homework
! RS: Let's talk about this tomorrow

! RS: Do you plan to change this later?
!                       KF: yes, it will be the next push, most likely
!       gridline = gridline_start + (grid_index-1)*gridline_spacing
        if (gridline < A(A_index,var)) then
                grid(grid_index) = A_index
                grid_index = grid_index + 1
                gridline = gridline + gridline_spacing
        else if (A_index == end_index) then
!               do
!                       if (grid_index > max_gridlines) exit
!                       grid(grid_index) = A_index + 1
!                       grid_index = grid_index + 1
!               end do
!               exit
                A_index = A_index + 1
                do i = grid_index, max_gridlines
                        grid(i) = A_index
                end do
                exit
        else
                A_index = A_index + 1
        end if
end do

end subroutine grider



subroutine chooseINT(Nintegers,lowerBound,upperBound,randomIntegers)
implicit none
integer,intent(in) :: Nintegers, lowerBound, upperBound
integer,dimension(Nintegers),intent(out) :: randomIntegers
real :: rand_num
integer :: rand_int
integer :: Nrange
logical :: stop_flag
integer :: i,j

Nrange = upperBound - lowerBound + 1

do i = 1, Nintegers
        do
                stop_flag = .false.

!
! for GNU fortran:
!               rand_num = rand()
! for intel fortran:
                call random_number(rand_num)
!
!
                rand_int = floor(rand_num*Nrange)
                if ((rand_int < 0).or.(rand_int == Nrange)) cycle
                do j = 1, i-1
                        if (rand_int == randomIntegers(j)) stop_flag = .true.
                end do
                if (stop_flag) cycle
                randomIntegers(i) = rand_int
                exit
        end do
end do

randomIntegers = randomIntegers + lowerBound

end subroutine chooseINT

subroutine QRDECOMP(A,m,n,Q,R)
    implicit none
    integer,intent(in) :: m,n
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(m,m), intent(out) :: Q
    real(dp), dimension(m,n), intent(out) :: R
    real(dp), dimension(m,m) :: Qn
    real(dp), dimension(m,n) :: Rn

    real(dp), dimension(m) :: x
    real(dp) :: u, xnorm, tau
    integer :: i, j, t

    do i = 1, m
        do j = 1, m
            if (i == j) then
                Q(i,j) = 1
            else
                Q(i,j) = 0
            end if
        end do
    end do

    Rn = A

!print *, "QRstart"
!print *, "need to do ", min(m-1,n), "rows"
    do t = 1, min(m-1,n)
!print *, "row", t

        do i = 1, m
            do j = 1, m
                if (i == j) then
                    Qn(i,j) = 1
                else
                    Qn(i,j) = 0
                end if
            end do
        end do

        x = Rn(:,t)
        xnorm = sqrt(sum(x(t:m)**2,dim=1))
        x(t) = x(t) - sign(xnorm,Rn(t,t))
        xnorm = sqrt(sum(x(t:m)**2,dim=1))
        x = x / xnorm

        Qn(t:m,t:m) = Qn(t:m,t:m) - 2*matmul(&
                reshape(x(t:m),(/m-t+1,1/)),&
                reshape(x(t:m),(/1,m-t+1/)))

        Rn = matmul(transpose(Qn),Rn)

        Q = matmul(Q,transpose(Qn))

    end do
!print *, "QRend"

    R = matmul(transpose(Q),A)

end subroutine QRDECOMP

subroutine testQRDECOMP()
    implicit none
    integer,parameter :: m = 3
    integer,parameter :: n = 3
    character(1) :: ntext
    character(1) :: mtext
    real(dp), dimension(m,n) :: A, B, R
    real(dp), dimension(m,m) :: Qtotal, Q

    integer i,j,k

    write(mtext,FMT="(I1)") m
    write(ntext,FMT="(I1)") n

    do i = 1, m
        do j = 1, n
!
! for GNU fortran:
!           A(i,j) = rand()
! for intel fortran:
            call random_number(A(i,j))
!
!
            R(i,j) = A(i,j)

            if (i == j) then
                Qtotal(i,j) = 1
            else
                Qtotal(i,j) = 0
            end if
        end do
    end do

    A(1,1) = 1.0d0
    A(1,2) = 2.0d0
    A(1,3) = 3.0d0

    A(2,1) = 2.0d0
    A(2,2) = 4.0d0
    A(2,3) = 6.0d0

    print *, ""
    print *, "A"
    print *, ""

    do i = 1, n
        write(6,FMT="("//mtext//"(F14.8))") A(:,i)
    end do

    call QRDECOMP(A,m,n,Q,R)

    A = matmul(transpose(Q),A)
    B = matmul(Q,R)

    print *, ""
    print *, "Q"
    print *, ""

    do j = 1, m
        write(6,FMT="("//mtext//"(F14.8))") Q(:,j)
    end do

    print *, ""
    print *, "R"
    print *, ""

    do j = 1, n
        write(6,FMT="("//mtext//"(F14.8))") R(:,j)
    end do

    print *, ""
    print *, "QR"
    print *, ""

    do j = 1, n
        write(6,FMT="("//mtext//"(F14.8))") B(:,j)
    end do

    print *, ""
    print *, "QtA"
    print *, ""

    do j = 1, n
        write(6,FMT="("//mtext//"(F14.8))") A(:,j)
    end do

end subroutine testQRDECOMP

subroutine CLS(A,m,n,C,p,d,b,x)
    implicit none
    integer,intent(in) :: m, n, p
    character(1) :: ntext
    character(1) :: ptext
    character(1) :: mptext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(p,n), intent(in) :: C
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(p), intent(in) :: d
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(m+p,n) :: S,R
    real(dp), dimension(m,m+p) :: Q1,Q1R
    real(dp), dimension(p,m+p) :: Q2,Q2R
    real(dp), dimension(m+p,m+p) :: Q2Tq
    real(dp), dimension(m+p,p) :: Q2T,Q2Tr
    real(dp), dimension(m+p,m+p) :: Q,Qt,QQ
    real(dp), dimension(m+p,n) :: QR
    real(dp), dimension(m+p) :: u,ce
    real(dp), dimension(p) :: w
    real(dp), dimension(m+p) :: y
    real(dp), dimension(m) :: bout
    real(dp), dimension(p,1) :: dout
    real(dp), dimension(m+p,1) :: ceout

    logical :: error_flag

    integer :: i, j

    write(ntext,FMT="(I1)") n
    write(ptext,FMT="(I1)") p
    write(mptext,FMT="(I1)") m+p

    do i = 1, m
        do j = 1, n
            S(i,j) = A(i,j)
        end do
        y(i) = b(i)
    end do

    do i = 1, p
        do j = 1, n
            S(m+i,j) = C(i,j)
        end do
        y(m+i) = d(i)
    end do

    call QRDECOMP(S,m+p,n,Q,R)

    print *, ""
    print *, "S ( 'stacked' [ A^T C^T ]^T )"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ntext//"(F10.6))") S(i,:)
    end do

    do i = 1, n
        if (all(abs(R(:,n)) < 1.0d-9)) then
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""
            print *, "ERROR: S is not left-invertible"
            print *, " (linearly independent columns)"
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""
        end if
    end do

    Qt = transpose(Q)

    print *, ""
    print *, "Q"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//mptext//"(F10.6))") Q(i,:)
    end do

    print *, ""
    print *, "R"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ntext//"(F10.6))") R(i,:)
    end do

    QQ = matmul(Qt,Q)

    print *, ""
    print *, "Q^T Q"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//mptext//"(F10.6))") QQ(i,:)
    end do

    QR = matmul(Q,R)

    print *, ""
    print *, "QR"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ntext//"(F10.6))") QR(i,:)
    end do

    Q1 = Q(1:m,:)

    print *, ""
    print *, "Q1"
    print *, ""

    do i = 1, m
        write(6,FMT="("//mptext//"(F10.6))") Q1(i,:)
    end do

    Q1R = matmul(Q1,R)

    print *, ""
    print *, "Q1R"
    print *, ""

    do i = 1, m
        write(6,FMT="("//ntext//"(F10.6))") Q1R(i,:)
    end do

    Q2 = Q(m+1:m+p,:)

    print *, ""
    print *, "Q2"
    print *, ""

    do i = 1, p
        write(6,FMT="("//mptext//"(F10.6))") Q2(i,:)
    end do

    Q2R = matmul(Q2,R)

    print *, ""
    print *, "Q2R"
    print *, ""

    do i = 1, p
        write(6,FMT="("//ntext//"(F10.6))") Q2R(i,:)
    end do

    call QRDECOMP(transpose(Q2),m+p,p,Q2Tq,Q2Tr)

    print *, ""
    print *, "Q2Tq"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//mptext//"(F10.6))") Q2Tq(i,:)
    end do

    print *, ""
    print *, "Q2Tr"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ptext//"(F10.6))") Q2Tr(i,:)
    end do

    Q2T = matmul(Q2Tq,Q2Tr)

    print *, ""
    print *, "Q2Tq Q2Tr"
    print *, ""

    do i = 1, m+p
        write(6,FMT="("//ptext//"(F10.6))") Q2T(i,:)
    end do

    if (abs(R(m+p,m+p-1)) < 1.0d-10) then
        print *, ""
        print *, "!!!!!!!!!!!!!!!!!!!!!!!"
        print *, ""
        print *, "Exact solution found!"
        print *, ""
        print *, "!!!!!!!!!!!!!!!!!!!!!!!"
        print *, ""

        print *, ""
        print *, "d"
        print *, ""

        do i = 1, p
            write(6,FMT="(F10.6)") d(i)
        end do

        call BackSubstitute(&
            transpose(Q2Tr),p,m+p,u,d,&
            error_flag)

        dout = matmul(transpose(Q2Tr),&
                  reshape(u,(/m+p,1/)))

        print *, ""
        print *, "dout"
        print *, ""

        do i = 1, p
            write(6,FMT="(F10.6)") dout(i,1)
        end do

!       ce = reshape(&
!            matmul(transpose(Q2Tq),&
!            matmul(transpose(Q1),&
!            reshape(b,(/m,1/)))),&
!            (/m+p/)) - u

        ce  = -u

        print *, ""
        print *, "ce"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") ce(i)
        end do

        call ForwardSubstitute(Q2Tr,m+p,p,&
        w,reshape(ce,(/m+p/)))

!       w = -1.0d0

        ceout = matmul(Q2Tr,&
                   reshape(w,(/p,1/)))

        print *, ""
        print *, "ceout"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") ceout(i,1)
        end do

        y = reshape(&
            matmul(transpose(Q1),&
                reshape(b,(/m,1/))) - &
            matmul(transpose(Q2),&
                reshape(w,(/p,1/))),&
            (/m+p/))

        print *, ""
        print *, "y"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") y(i)
        end do

        call ForwardSubstitute(R,m+p,n,x,y)

        print *, ""
        print *, "x"
        print *, ""

        do i = 1, n
            write(6,FMT="(F10.6)") x(i)
        end do

        y = reshape(matmul(R,&
                reshape(x,(/n,1/))),&
                (/m+p/))

        print *, ""
        print *, "R x"
        print *, ""

        do i = 1, m+p
            write(6,FMT="(F10.6)") y(i)
        end do



        bout = matmul(A,x) - b

        print *, ""
        print *, "A x - b"
        print *, ""

        do i = 1, m
            write(6,FMT="(F10.6)") bout(i)
        end do

        x = 1.0d0 / n

        bout = matmul(A,x) - b

        print *, ""
        print *, "suppose x = 1/n"
        print *, ""
        print *, "A x - b"
        print *, ""

        do i = 1, m
            write(6,FMT="(F10.6)") bout(i)
        end do

    end if

end subroutine CLS

subroutine CLS2(A,m,n,C,p,d,b,x)
    implicit none
    integer,intent(in) :: m, n, p
    character(2) :: ntext
    character(2) :: ptext
    character(2) :: mptext
    character(2) :: nptext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(p,n), intent(in) :: C
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(p), intent(in) :: d
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(n+p,n+p) :: S,R
    real(dp), dimension(n+p,n+p) :: Q,Qt,QQ
    real(dp), dimension(n+p,n+p) :: QR
    real(dp), dimension(n,n) :: H
    real(dp), dimension(n+p,n+p) :: E
    real(dp), dimension(n,1) :: z
    real(dp), dimension(n+p,1) :: w
    real(dp), dimension(n+p) :: y, u
    real(dp), dimension(m,1) :: bout

    integer :: i, j

!   write(ntext,FMT="(I2)") n
!   write(ptext,FMT="(I2)") p
!   write(mptext,FMT="(I2)") m+p
!   write(nptext,FMT="(I2)") n+p

!   print *, ""
!   print *, "A"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="("//ntext//"(F14.10))") A(i,:)
!   end do


!   print *, ""
!   print *, "b"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F14.10)") b(i)
!   end do

    H = matmul(transpose(A),A)

    E(1:n,1:n) = H
    E(1:n,n+1:n+p) = transpose(C)
    E(n+1:n+p,1:n) = C
    E(n+1:n+p,n+1:n+p) = 0.0d0

    call QRDECOMP(E,n+p,n+p,Q,R)

!   print *, ""
!   print *, "E"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") E(i,:)
!   end do

!    do i = 1, n
!        if (all(abs(R(:,n)) < 1.0d-9)) then
!            print *, ""
!            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!            print *, ""
!            print *, "ERROR: S is not left-invertible"
!            print *, " (linearly independent columns)"
!            print *, ""
!            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!            print *, ""
!        end if
!    end do

    Qt = transpose(Q)

!   print *, ""
!   print *, "Q"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") Q(i,:)
!   end do

!   print *, ""
!   print *, "R"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") R(i,:)
!   end do

!   QQ = matmul(Qt,Q)

!   print *, ""
!   print *, "Q^T Q"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") QQ(i,:)
!   end do

!   QR = matmul(Q,R)

!   print *, ""
!   print *, "QR"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="("//nptext//"(F10.6))") QR(i,:)
!   end do

    z = matmul(transpose(A),&
           reshape(b,(/m,1/)))

    w(1:n,1) = z(1:n,1)
    w(n+1:n+p,1) = d(:)

!   print *, ""
!   print *, "w"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") w(i,1)
!   end do

    w = matmul(Qt,w)
    y = reshape(w,(/n+p/))

!   print *, ""
!   print *, "y = Qt w"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") y(i)
!   end do

    call ForwardSubstitute(R,n+p,n+p,u,y)

!   print *, ""
!   print *, "u"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") u(i)
!   end do

!   w = matmul(R,reshape(u,(/n+p,1/)))
!   y = reshape(w,(/n+p/))

!   print *, ""
!   print *, "y = R u"
!   print *, ""

!   do i = 1, n+p
!       write(6,FMT="(F10.6)") y(i)
!   end do

    x = u(1:n)

!   bout = matmul(A,reshape(x,(/n,1/)))

!   print *, ""
!   print *, "bout = A x"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F14.10)") bout(i,1)
!   end do

!   print *, ""
!   print *, "||bout-b|| = ", sqrt(sum((bout(:,1)-b)**2))
!   print *, "||bout-b|| = ", sqrt(sum((bout(1:12,1)-b(1:12))**2)/4)
!   print *, ""

!   x = 1.0d0 / n

!   bout = matmul(A,reshape(x,(/n,1/)))

!   print *, ""
!   print *, "suppose x = 1/n"
!   print *, ""
!   print *, "bout = A x"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F10.6)") bout(i,1)
!   end do

!   print *, ""
!   print *, "||bout-b|| = ", sqrt(sum((bout(:,1)-b)**2))
!   print *, "||bout-b|| = ", sqrt(sum((bout(1:12,1)-b(1:12))**2)/4)
!   print *, ""

!   x = u(1:n)
!   x(1) = x(1) + 0.01d0 * u(1)
!   x = x / (sum(x))

!   bout = matmul(A,reshape(x,(/n,1/)))

!   print *, ""
!   print *, "suppose x = normalized(x + 0.01 * x1 * e1)"
!   print *, ""
!   print *, "bout = A x"
!   print *, ""

!   do i = 1, m
!       write(6,FMT="(F10.6)") bout(i,1)
!   end do

!   print *, ""
!   print *, "||bout-b|| = ", sqrt(sum((bout(:,1)-b)**2))
!   print *, "||bout-b|| = ", sqrt(sum((bout(1:12,1)-b(1:12))**2)/4)
!   print *, ""

!   x = u(1:n)

end subroutine CLS2

subroutine CLSwK(A,m,n,C,p,d,b,x)
    implicit none
    integer,intent(in) :: m, n, p
    character(2) :: ntext
    character(2) :: ptext
    character(2) :: mptext
    character(2) :: nptext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(p,n), intent(in) :: C
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(p), intent(in) :: d
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(n+p,n+p) :: S,R
    real(dp), dimension(n+p,n+p) :: Q,Qt,QQ
    real(dp), dimension(n+p,n+p) :: QR
    real(dp), dimension(n,n) :: H
    real(dp), dimension(n+p,n+p) :: E
    real(dp), dimension(n,1) :: z
    real(dp), dimension(n+p,1) :: w
    real(dp), dimension(n+p) :: y, u
    real(dp), dimension(m,1) :: bout

    integer :: kID = 1
    logical :: use_constraint = .false.

    integer :: i, j

    select case(kID)

    ! The original CLS2 algorithm:
    case(0)

        H = matmul(transpose(A),A)
        z = matmul(transpose(A),&
               reshape(b,(/m,1/)))

    ! Just take the squared distance bteween
    ! the two inputs
    case(1)

        do i = 1, n
            H(i,i) = 0.0d0
            do j = i+1, n
                H(i,j) = sum((A(:,i)-A(:,j))**2)
                H(j,i) = H(i,j)
            end do
        end do

        do i = 1, n
            z(i,1) = sum((A(:,i)-b)**2)
        end do

    ! The Matern kernel (as in the sGDML paper)
   !case (2)

    end select

    if (use_constraint) then
        E(1:n,1:n) = H
        E(1:n,n+1:n+p) = transpose(C)
        E(n+1:n+p,1:n) = C
        E(n+1:n+p,n+1:n+p) = 0.0d0

        call QRDECOMP(E,n+p,n+p,Q,R)

        Qt = transpose(Q)

        w(1:n,1) = z(1:n,1)
        w(n+1:n+p,1) = d(:)

        w = matmul(Qt,w)
        y = reshape(w,(/n+p/))

        call ForwardSubstitute(R,n+p,n+p,u,y)

        x = u(1:n)

    else
        call QRDECOMP(H,n,n,Q(1:n,1:n),R(1:n,1:n))

        Qt(1:n,1:n) = transpose(Q(1:n,1:n))

        w(1:n,1:1) = matmul(Qt(1:n,1:n),z)

        call ForwardSubstitute(R(1:n,1:n),n,n,x,w(1:n,1))

    end if

end subroutine CLSwK


subroutine LS(A,m,n,b,x)
    implicit none
    integer,intent(in) :: m, n
    character(2) :: ntext
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(m), intent(in) :: b
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(n,n) :: H, Q, R, Qt
    real(dp), dimension(n,1) :: z
    real(dp), dimension(n,1) :: w
    real(dp), dimension(n) :: y, u

    integer :: i, j

    H = matmul(transpose(A),A)

    call QRDECOMP(H,n,n,Q,R)

    Qt = transpose(Q)

    z = matmul(transpose(A),&
           reshape(b,(/m,1/)))

    w(1:n,1) = z(1:n,1)

    w = matmul(Qt,w)
    y = reshape(w,(/n/))

    call ForwardSubstitute(R,n,n,u,y)

    x = u(1:n)

end subroutine LS

subroutine squareAxb(A,n,b,x)
    implicit none
    integer,intent(in) :: n
    real(dp), dimension(n,n), intent(in) :: A
    real(dp), dimension(n), intent(in) :: b
    real(dp), dimension(n), intent(out) :: x

    real(dp), dimension(n,n) :: Q, R, Qt
    real(dp), dimension(n,1) :: w
    real(dp), dimension(n) :: y, u

    integer :: i, j

    call QRDECOMP(A,n,n,Q,R)

    Qt = transpose(Q)
    w = reshape(b,(/n,1/))
    w = matmul(Qt,w)
    y = reshape(w,(/n/))

    call ForwardSubstitute(R,n,n,u,y)

    x = u(1:n)

end subroutine squareAxb

! Minimize the Frobenius norm
subroutine mFN(Omega,m,n,A,k,B)
    implicit none
    integer,intent(in) :: m, n, k
    real(dp), dimension(n,k), intent(in) :: A
    real(dp), dimension(m,k), intent(in) :: B
    real(dp), dimension(m,n), intent(out) :: Omega

    real(dp), dimension(n,n) :: H
    real(dp), dimension(m,n) :: G, OmegaP
    real(dp), dimension(n,n) :: L, R, Q, Qt

    logical :: error_flag

    integer :: i, j
    
    H = matmul(A, transpose(A))
    G = matmul(B, transpose(A))

    call QRDECOMP(H,n,n,Q,R)
    Qt = transpose(Q)
    L = transpose(R)

!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   print *, "A:"
!   do i = 1, n
!       print *, i, A(i,1:k)
!   end do

!   print *, "Q:"
!   do i = 1, n
!       print *, i, Q(i,1:n)
!   end do

!   print *, "R:"
!   do i = 1, n
!       print *, i, R(i,1:k)
!   end do

!   H = abs(matmul(Q,R) - A)
!   print *, "|QR-A|:"
!   do i = 1, n
!       print *, i, H(i,1:k)
!   end do
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do i = 1, m
        call BackSubstitute(L,n,n,&
                OmegaP(i,1:n),G(i,1:n),&
                error_flag)
        if (error_flag) exit
    end do

    if (error_flag) then
        print *, "ERROR in MNF"
        Omega = 0.0d0
        do i = 1, m
            Omega(i,i) = 1.0d0
        end do
    else
        Omega = matmul(OmegaP,Qt)
    end if

end subroutine mFN

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solve for the vector-valued function as
! an optimization problem where:
!
! (1) A will be of size mxn containing the
!   n descriptors of size
!   m
!
! (2) Omega will be of size mxmn with an
!   1xn block matrix structure with
!   mxm matrix-valued kernels for each
!
! (3) The kernel can be separated into the
!     sum of two components, the:
!   (dx)(dx^T) component and the
!   I_mn component
!     where dx = x1 - x2 for k(x1,x2)
subroutine sVVF(A,m,n,a0,Omega)
    implicit none
    integer,intent(in) :: m, n
!   real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(m*n), intent(in) :: A
    real(dp), dimension(m) :: a0
    real(dp), dimension(m,m*n), intent(out) :: Omega

    real(dp), dimension(m,m*n) :: Wtmp
    real(dp),dimension(m,m) :: Im
    real(dp),dimension(m*n,m*n) :: K, Q, R, Qt, L
    real(dp),dimension(m) :: dx
    real(dp) :: d, OPc, IPc

    integer :: startindex1,endindex1
    integer :: startindex2,endindex2

    logical :: error_flag

    integer :: i, j

    Im = 0.0d0
    do i = 1, m
        Im(i,i) = 1.0d0
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        dx = A(startindex1:endindex1) - a0
        d = sqrt(sum(dx**2))
        call sOPIP(d,OPc,IPc)
        Omega(1:m,startindex1:endindex1) = &
            matmul(reshape(dx,(/ m, 1 /)),&
                   reshape(dx,(/ 1, m /))) * OPc + &
            Im * IPc
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        startindex2 = startindex1
        endindex2 = endindex1
        do j = i, n
!           dx = A(1:m,i)-A(1:m,j)
            dx = A(startindex1:endindex1) - &
                 A(startindex2:endindex2)
            d = sqrt(sum(dx**2))
            call sOPIP(d,OPc,IPc)
            K(startindex1:endindex1,&
              startindex2:endindex2) = &
                matmul(reshape(dx,(/ m, 1 /)),&
                       reshape(dx,(/ 1, m /))) * OPc + &
                Im * IPc

            K(startindex2:endindex2,&
              startindex1:endindex1) = &
                K(startindex1:endindex1,&
                  startindex2:endindex2)

            startindex2 = startindex2 + m
            endindex2 = endindex2 + m
        end do
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do

    call QRDECOMP(K,m*n,m*n,Q,R)
    Qt = transpose(Q)
    L = transpose(R)

    do i = 1, m
        call BackSubstitute(L,m*n,m*n,&
                Wtmp(i,1:m*n),Omega(i,1:m*n),&
                error_flag)
        if (error_flag) exit
    end do

    if (error_flag) then
        print *, "ERROR in sVVF"
        Omega = 0.0d0
        do i = 1, m
            Omega(i,i) = 1.0d0
        end do
    else
        Omega = matmul(Wtmp,Qt)
    end if

end subroutine sVVF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Solve for the vector-valued function as
! an optimization problem where:
!
! (1) A will be of size mxn containing the
!   n descriptors of size
!   m
!
! (2) Omega will be of size mxmn with an
!   1xn block matrix structure with
!   mxm matrix-valued kernels for each
!
! (3) The kernel can be separated into the
!     sum of two components, the:
!   (dx)(dx^T) component and the
!   I_mn component
!     where dx = x1 - x2 for k(x1,x2)
subroutine sVVFwJ(A,dA,m,d,n,a0,da0,Omega)
    implicit none
    integer,intent(in) :: m, d, n
    real(dp), dimension(d,n), intent(in) :: A
    real(dp), dimension(d,m*n), intent(in) :: dA
    real(dp), dimension(d), intent(in) :: a0
    real(dp), dimension(d,m), intent(in) :: da0
    real(dp), dimension(m,m*n), intent(out) :: Omega

    real(dp), dimension(m,m*n) :: Wtmp
    real(dp),dimension(d,d) :: Id
    real(dp),dimension(m*n,m*n) :: K, Q, R, Qt, L
    real(dp),dimension(d) :: dx
    real(dp) :: dist, OPc, IPc

    real(dp),parameter :: lam = 1.0d-15

    integer :: startindex1,endindex1
    integer :: startindex2,endindex2

    logical :: error_flag

    integer,save :: sVVFwJ_counter = 0
    character(15) :: sVVFwJ_filename

    integer :: i, j

    Id = 0.0d0
    do i = 1, d
        Id(i,i) = 1.0d0
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        dx = A(1:d,i) - a0
        dist = sqrt(sum(dx**2))
        call sOPIP(dist,OPc,IPc)
        Omega(1:m,startindex1:endindex1) = &
            matmul(matmul(transpose(da0),&
           (&
            matmul(reshape(dx,(/ d, 1 /)),&
                   reshape(dx,(/ 1, d /))) * OPc + &
            Id * IPc &
           )&
                         ),dA(1:d,startindex1:endindex1))
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        startindex2 = startindex1 !1
        endindex2 = endindex1 !m
        do j = i, n
!           dx = A(1:m,i)-A(1:m,j)
            dx = A(1:d,i) - A(1:d,j)
            dist = sqrt(sum(dx**2))
            call sOPIP(dist,OPc,IPc)
            K(startindex1:endindex1,&
              startindex2:endindex2) = &
                matmul(matmul(&
                transpose(dA(1:d,startindex1:endindex1)),&
               (&
                matmul(reshape(dx,(/ d, 1 /)),&
                       reshape(dx,(/ 1, d /))) * OPc + &
                Id * IPc &
               )&
                             ),dA(1:d,startindex2:endindex2))

            K(startindex2:endindex2,&
              startindex1:endindex1) = transpose(&
                K(startindex1:endindex1,&
                  startindex2:endindex2)        )

            startindex2 = startindex2 + m
            endindex2 = endindex2 + m
        end do
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do




!   sVVFwJ_counter = sVVFwJ_counter + 1
!   if (modulo(sVVFwJ_counter,1200) == 1) then
!   do i = 1, n
!       write(sVVFwJ_filename,&
!             FMT="('K',I0.2,'_',I0.7,'.dat')") &
!           i, sVVFwJ_counter / 12
!       open(6666,file=&
!                 "/home/kazuumi/rsun_lts/kazuumi"//&
!                 "/venus-NEW/SEP01_test3_again_10/"//&
!                 sVVFwJ_filename)
!       do j = 1, m
!           write(6666,FMT="(15(ES14.7,1x))") &
!               K((i-1)*m+j,1+(i-1)*m:i*m) 
!       end do
!       close(6666)
!   end do

!   write(sVVFwJ_filename,&
!         FMT="('Q',I0.2,'_',I0.7,'.dat')") &
!       0, sVVFwJ_counter / 12
!   open(6666,file=&
!             "/home/kazuumi/rsun_lts/kazuumi"//&
!             "/venus-NEW/SEP01_test3_again_10/"//&
!             sVVFwJ_filename)
!   do i = 1, n
!       write(6666,FMT="(10(ES14.7,1x))") A(1:d,i)
!   end do
!   close(6666)
!   end if




    ! Add the regularization
    do i = 1, m*n
        K(i,i) = K(i,i) - lam
    end do

    call QRDECOMP(K,m*n,m*n,Q,R)
    Qt = transpose(Q)
    L = transpose(R)

    do i = 1, m
        call BackSubstitute(L,m*n,m*n,&
                Wtmp(i,1:m*n),Omega(i,1:m*n),&
                error_flag)
        if (error_flag) exit
    end do

    if (error_flag) then
        print *, "ERROR in sVVF"
        Omega = 0.0d0
        do i = 1, m
            Omega(i,i) = 1.0d0
        end do
    else
        Omega = matmul(Wtmp,Qt)
    end if

end subroutine sVVFwJ


subroutine storeRKHSInverseKernel(A,dA,m,d,n,F)
    implicit none
    integer,intent(in) :: m, d, n
    real(dp), dimension(d,n), intent(in) :: A
    real(dp), dimension(d,m*n), intent(in) :: dA
    real(dp), dimension(m,n), intent(in) :: F

    real(dp),dimension(d,d) :: Id
    real(dp),dimension(m*n,m*n) :: K, Q, R
!   real(dp),dimension(m*n,m*n) :: Qt, L
    real(dp),dimension(d) :: dx
    real(dp) :: dist, OPc, IPc

    real(dp),parameter :: lam = 1.0d-15

    integer :: startindex1,endindex1
    integer :: startindex2,endindex2

    logical :: error_flag

    integer,save :: sVVFwJ_counter = 0
    character(15) :: sVVFwJ_filename

    integer :: i, j

    print *, "this is Kazuumi 0"

    if (RKHSarrays_are_not_allocated) then
        RKHSarrays_are_not_allocated = .false.
        allocate(RKHS_A(10,1000),&
                 RKHS_dA(10,12000),&
                 RKHS_F(12000),&
                 RKHS_Qt(12000,12000),&
                 RKHS_L(12000,12000))
    end if

    print *, "this is Kazuumi 1"

    Id = 0.0d0
    do i = 1, d
        Id(i,i) = 1.0d0
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        startindex2 = startindex1 !1
        endindex2 = endindex1 !m
        do j = i, n
!           dx = A(1:m,i)-A(1:m,j)
            dx = A(1:d,i) - A(1:d,j)
            dist = sqrt(sum(dx**2))
            call sOPIP(dist,OPc,IPc)
            K(startindex1:endindex1,&
              startindex2:endindex2) = &
                matmul(matmul(&
                transpose(dA(1:d,startindex1:endindex1)),&
               (&
                matmul(reshape(dx,(/ d, 1 /)),&
                       reshape(dx,(/ 1, d /))) * OPc + &
                Id * IPc &
               )&
                             ),dA(1:d,startindex2:endindex2))

            K(startindex2:endindex2,&
              startindex1:endindex1) = transpose(&
                K(startindex1:endindex1,&
                  startindex2:endindex2)        )

            startindex2 = startindex2 + m
            endindex2 = endindex2 + m
        end do
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do

    print *, "this is Kazuumi 2"

    ! Add the regularization
    do i = 1, m*n
        K(i,i) = K(i,i) - lam
    end do

    call QRDECOMP(K,m*n,m*n,Q,R)
    RKHS_Qt(1:m*n,1:m*n) = transpose(Q)
    RKHS_L(1:m*n,1:m*n) = transpose(R)

    RKHS_A(1:d,1:n) = A
    RKHS_dA(1:d,1:m*n) = dA

    print *, "this is Kazuumi 3"

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        RKHS_F(startindex1:endindex1) = F(1:m,i)
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do
!   RKHS_F(1:m*n) = reshape(F,(/m*n/))

    return
end subroutine storeRKHSInverseKernel

!subroutine getRKHSTarget(m,d,n,a0,da0,Omega)
subroutine getRKHSTarget(m,d,n,a0,da0,F)
    implicit none
    integer,intent(in) :: m, d, n
    real(dp), dimension(d), intent(in) :: a0
    real(dp), dimension(d,m), intent(in) :: da0
!   real(dp), dimension(m,m*n), intent(out) :: Omega
    real(dp), dimension(m,m*n) :: Omega
    real(dp), dimension(m), intent(out) :: F

    real(dp), dimension(m,m*n) :: Wtmp
    real(dp),dimension(d,d) :: Id
    real(dp),dimension(m*n,m*n) :: K, Q, R
!   real(dp),dimension(m*n,m*n) :: Qt, L
    real(dp),dimension(d) :: dx
    real(dp) :: dist, OPc, IPc

    real(dp),parameter :: lam = 1.0d-15

    integer :: startindex1,endindex1
    integer :: startindex2,endindex2

    logical :: error_flag

    integer,save :: sVVFwJ_counter = 0
    character(15) :: sVVFwJ_filename

    integer :: i, j

    Id = 0.0d0
    do i = 1, d
        Id(i,i) = 1.0d0
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        dx = RKHS_A(1:d,i) - a0
        dist = sqrt(sum(dx**2))
        call sOPIP(dist,OPc,IPc)
        Omega(1:m,startindex1:endindex1) = &
            matmul(matmul(transpose(da0),&
           (&
            matmul(reshape(dx,(/ d, 1 /)),&
                   reshape(dx,(/ 1, d /))) * OPc + &
            Id * IPc &
           )&
                         ),RKHS_dA(1:d,startindex1:endindex1))
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m
    end do

    do i = 1, m
        call BackSubstitute(RKHS_L(1:m*n,1:m*n),m*n,m*n,&
                Wtmp(i,1:m*n),Omega(i,1:m*n),&
                error_flag)
        if (error_flag) exit
    end do

    if (error_flag) then
        print *, "ERROR in sVVF"
        Omega = 0.0d0
        do i = 1, m
            Omega(i,i) = 1.0d0
        end do
    else
        Omega = matmul(Wtmp,RKHS_Qt(1:m*n,1:m*n))
    end if

    F = reshape(matmul(Omega,reshape(RKHS_F(1:m*n),(/m*n,1/))),(/m/))

    return
end subroutine getRKHSTarget

subroutine storeRKHSJacobian(A,dA,m,d,n,F)
    implicit none
    integer,intent(in) :: m, d, n
    real(dp), dimension(d,n), intent(in) :: A
    real(dp), dimension(d,m*n), intent(in) :: dA
    real(dp), dimension(m,n), intent(in) :: F

    if (RKHSarrays_are_not_allocated) then
        RKHSarrays_are_not_allocated = .false.
        allocate(RKHS_A(10,1000),&
                 RKHS_dA(10,12000),&
                 RKHS_F(12000))

        RKHS_A(1:d,1:n) = A
        RKHS_dA(1:d,1:m*n) = dA

        !open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/RKHSalpha.dat")
        !open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/RKHSalphaMulBySqrt2.dat")

        !open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/RKHSalpha-700train20k.dat")
        !open(6666,file="/home/kazuumi/rsun_lts/kazuumi/venus-NEW/RKHSalpha-700train20kMulByYstd.dat")

        open(6666,file="/home/kazuumi/rsun_lts/kazuumi/alpha9421774.dat")
        read(6666,FMT=*) RKHS_F(1:m*n)
        close(6666)
    end if

    return
end subroutine storeRKHSJacobian


subroutine getRKHSTarget_fromalpha(m,d,n,a0,da0,F)
    implicit none
    integer,intent(in) :: m, d, n
    real(dp), dimension(d), intent(in) :: a0
    real(dp), dimension(d,m), intent(in) :: da0
    real(dp), dimension(m,m*n) :: Omega
    real(dp), dimension(m), intent(out) :: F

    real(dp), dimension(m,m*n) :: Wtmp
    real(dp),dimension(d,d) :: Id
    real(dp),dimension(m*n,m*n) :: K, Q, R
    real(dp),dimension(d) :: dx
    real(dp) :: dist, OPc, IPc

    real(dp),parameter :: lam = 1.0d-15

    integer :: startindex1,endindex1
    integer :: startindex2,endindex2

    logical :: error_flag
    real(dp) :: dum1

    integer :: i, j

    dum1 = 1.0d24

    Id = 0.0d0
    do i = 1, d
        Id(i,i) = 1.0d0
    end do

    startindex1 = 1
    endindex1 = m
    do i = 1, n
        dx = RKHS_A(1:d,i) - a0
        dist = sqrt(sum(dx**2))
        call sOPIP(dist,OPc,IPc)
        Omega(1:m,startindex1:endindex1) = &
            matmul(matmul(transpose(da0),&
           (&
            matmul(reshape(dx,(/ d, 1 /)),&
                   reshape(dx,(/ 1, d /))) * OPc + &
            Id * IPc &
           )&
                         ),RKHS_dA(1:d,startindex1:endindex1))
        startindex1 = startindex1 + m
        endindex1 = endindex1 + m

        dum1 = min(dum1,dist)
    end do

    print *, "this is Kazuumi 1"
    print *, "min Jacobian dist (target, library):", dum1

    F = reshape(matmul(Omega,reshape(RKHS_F(1:m*n),(/m*n,1/))),(/m/))

    print *, "this is Kazuumi 2"

    return
end subroutine getRKHSTarget_fromalpha




! Return the scalar component for the two
! components of the kernel, the:
!   (dx)(dx^T) component and the
!   I_mn component
! of the outer product (OP) and inner
! product (IP) components
!
! Because the kernel is (assumed) stationary,
! this only takes the distance between the
! two as an argument
subroutine sOPIP(d,OPc,IPc)
    implicit none
    real(dp),intent(in) :: d
    real(dp),intent(out) :: OPc, IPc

    real(dp) :: s1 = 52.0d0 !1.0d0
    real(dp) :: dthreshold = 1.0d-18

    real(dp) :: d2,d3,invd,invd2,invd3
    real(dp) :: dum1,dum2,dum3,dum4,dum5,dum6

    real(dp),parameter :: sqrt5 = sqrt(5.0d0)

    integer :: kID = 2

    integer :: i, j

    select case(kID)

    ! The Matern kernel (as in the sGDML paper)
    case(1)

        if (d < dthreshold) then
            OPc = 0.0d0
            IPc = 0.0d0
            return
        end if

        d2 = d*d
        d3 = d2*d

        invd = 1.0d0 / d
        invd2 = invd**2
        invd3 = invd2*invd

        dum1 = 5.0d0 / s1
        dum2 = sqrt(dum1)
        dum3 = exp(-dum2 * d)

        dum4 = sqrt(2.0d0) * dum1

        dum5 = -(dum2*invd2 + invd3)
        dum6 = 1.0d0 / 24.0d0
  
        OPc = dum6 * (&
               (2)*(1)*(dum4**2) * (0 -dum2*( 2*invd3 + dum5*d2)) + &
               (6)*(2)*(dum4**1) * (-invd3 -dum2*( invd2 + dum5*d )) + &
               (24)*(1)*(dum4**0) *(0 -dum2*( 0 + dum5 )) &
              ) * dum3
        IPc = dum6 * (&
               (2)*(1)*(dum4**2) * d2 + &
               (6)*(2)*(dum4**1) * d + &
               (24)*(1)*(dum4**0) *1 &
              ) * dum3 * (-dum2) * invd

    ! The Matern kernel (as in the sGDML script
    !                    without permutation)
    case(2)

        dum1 = s1**2
        dum2 = (5.0d0 / (3 * (dum1**2))) * exp(-sqrt5 * d / s1)

        OPc = 5 * dum2
        IPc = - (dum1 + s1 * sqrt5 * d) * dum2

    end select
    
    return
end subroutine sOPIP



subroutine testCLS()
    implicit none
    integer, parameter :: m = 4
    integer, parameter :: p = 1
    integer, parameter :: n = 4
    character(1) :: mtext
    character(1) :: ptext
    character(1) :: ntext
    real(dp), dimension(m,n) :: A
    real(dp), dimension(p,n) :: C
    real(dp), dimension(m) :: b
    real(dp), dimension(p) :: d
    real(dp), dimension(n) :: x

    real(dp), dimension(m+p,n) :: S,R
    real(dp), dimension(m+p,m+p) :: Q

    integer :: i, j

    write(mtext,FMT="(I1)") m
    write(ptext,FMT="(I1)") p
    write(ntext,FMT="(I1)") n

    do i = 1, m
        do j = 1, n
!
! for GNU fortran:
!           A(i,j) = rand()
! for intel fortran:
            call random_number(A(i,j))
!
!
        end do
    end do

    do i = 1, p
        do j = 1, n
!
! for GNU fortran:
!           C(i,j) = rand()
! for intel fortran:
            call random_number(C(i,j))
!
!
        end do
    end do

    C = 1.0d0

    b = 0.0d0
    d = 1.0d0

    print *, ""
    print *, "A"
    print *, ""

    do i = 1, m
        write(6,FMT="("//ntext//"(F10.6))") A(i,:)
    end do

    print *, ""
    print *, "C"
    print *, ""

    do i = 1, p
        write(6,FMT="("//ntext//"(F10.6))") C(i,:)
    end do

    call CLS2(A,m,n,C,p,d,b,x)

    print *, ""
    print *, "x"
    print *, ""

    do i = 1, n
        write(6,FMT="(F10.6)") x(i)
    end do

    b = reshape(matmul(A,reshape(x,(/n,1/))),(/m/))

    print *, ""
    print *, "b = A x"
    print *, ""

    do i = 1, m
        write(6,FMT="(F10.6)") b(i)
    end do

    print *, ""
    print *, "||b|| = ", sqrt(sum(b**2))
    print *, ""

end subroutine testCLS

!This is for LOWER (LEFT) TRIANGULAR matrix A
!with equation of form Ax = b
subroutine BackSubstitute(A,m,n,x,b,error_flag)
    implicit none
    integer,intent(in) :: m, n
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(n), intent(out) :: x
    real(dp), dimension(m), intent(in) :: b
    logical,intent(out) :: error_flag
    
    real(dp) :: division_tolerance = 1.0d-30 ! 1.0d-15
    real(dp) :: other_terms

    integer :: i, j

    error_flag = .false.
    x = 0.0d0

    do i = 1, m
        if (abs(A(i,i)) < division_tolerance) then
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""
            print *, "ERROR in back substitution... EXITING"
            print *, ""
            print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            print *, ""

            error_flag = .true.
            exit
        end if

        other_terms = 0.0d0
        do j = 1, i-1
            other_terms = other_terms + x(j)*A(i,j)
        end do

        x(i) = (b(i) - other_terms) / A(i,i)
    end do

end subroutine BackSubstitute

!This is for UPPER (RIGHT) TRIANGULAR matrix A
!with equation of form Ax = b
subroutine ForwardSubstitute(A,m,n,x,b)
    implicit none
    integer,intent(in) :: m, n
    real(dp), dimension(m,n), intent(in) :: A
    real(dp), dimension(n), intent(out) :: x
    real(dp), dimension(m), intent(in) :: b
    
    real(dp) :: other_terms, dividor
    integer :: extra_pivots

    integer :: i, j

    extra_pivots = m+1
    x = 0.0d0

    do i = 1, m
        dividor = A(m-i+1,m-i+1)

        if (abs(dividor) < 1.0d-15) then
            if (extra_pivots <= n) then
                !need to add this
                extra_pivots = extra_pivots + 1
            end if

            other_terms = 0.0d0
            do j = (m-i+1), m-1
                other_terms = other_terms + x(j)*A(m-i+1,j)
            end do

            if (abs(other_terms - b(m-i+1)) < 1.d-10) then
                dividor = 1.0d0
            else
                print *, ""
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, ""
                print *, "ERROR in forward substitution... EXITING"
                print *, ""
                print *, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                print *, ""

                x = 0.0d0
                exit
            end if
        end if

        other_terms = 0.0d0
        do j = (m-i+1), m
            other_terms = other_terms + x(j)*A(m-i+1,j)
        end do

        x(m-i+1) = (b(m-i+1) - other_terms) / dividor
    end do

end subroutine ForwardSubstitute



subroutine getFlooredSqrt(i1,i2)
    implicit none
    integer,intent(in) :: i1
    integer,intent(out) :: i2

    if (i1 < 2) then
        i2 = i1
        return
    end if

    do i2 = 2, i1
        if (i2*i2 > i1) exit
    end do
    
    i2 = i2 - 1
    return

end subroutine getFlooredSqrt

subroutine ZtoN(z,n)
    implicit none
    integer,intent(in) :: z
    integer,intent(out) :: n

    if (z <= 0) then
        n = 2 * abs(z)
    else
        n = 2*z - 1
    end if

    return

end subroutine ZtoN

subroutine NNtoN(n1,n2,n)
    implicit none
    integer,intent(in) :: n1,n2
    integer,intent(out) :: n
    integer :: tempn

    tempn = n1 + n2
    
    n = n1 + ((tempn)*(tempn+1))/2

    return

end subroutine NNtoN

subroutine NtoZ(n,z)
    implicit none
    integer,intent(in) :: n
    integer,intent(out) :: z

    if (modulo(n,2)==0) then
        z = - n / 2
    else
        z = (n + 1) / 2
    end if

    return

end subroutine NtoZ

subroutine NtoNN(n,n1,n2)
    implicit none
    integer,intent(in) :: n
    integer,intent(out) :: n1,n2
    integer :: tempn

    call getFlooredSqrt(8*n+1,tempn)
    tempn = (tempn - 1)/2
    
    n1 = n - ((tempn)*(tempn+1))/2
    n2 = ((tempn)*(tempn+3))/2 - n

    return

end subroutine NtoNN


subroutine getFlattened(Nvar,x,xflat)
    implicit none
    integer,intent(in) :: Nvar
    integer,dimension(Nvar),intent(in) :: x
    integer,dimension(Nvar) :: xPOS
    integer :: tempx,i
    integer,intent(out) :: xflat

    do i = 1, Nvar
        call ZtoN(x(i),xPOS(i))
    end do

    tempx = xPOS(1)

    do i = 2, Nvar
        call NNtoN(tempx,xPOS(i),xflat)
        tempx = xflat
    end do

    xflat = xflat + 1

    return

end subroutine getFlattened

subroutine getExpanded(Nvar,x,xexpanded)
    implicit none
    integer,intent(in) :: Nvar, x
    integer :: tempx1,tempx2,i
    integer,dimension(Nvar),intent(out) :: xexpanded

    tempx1 = x - 1
    do i = Nvar, 2, -1
        call NtoNN(tempx1,tempx2,xexpanded(i))
        tempx1 = tempx2
    end do

    xexpanded(1) = tempx1

    do i = 1, Nvar
        tempx1 = xexpanded(i)
        call NtoZ(tempx1,xexpanded(i))
    end do

    return

end subroutine getExpanded


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! A more-or-less direct copy and paste of faolane's
! Fortran90 Munkre's Algorithm implementation. A few
! minor adjustments have been made to:
!  (1) adhere to my coding style
!  (2) accept a cost matrix C of type real(dp)
!  (3) assume all costs are non-negative and
!             that we are minimizing cost
!  (4) remove all the extra stuff
!
! Here is a link to the github page this came from:
!    https://github.com/faolane/LAP
! It's a very nice code! It is fairly well-commented
! and the flow of the algorithm is clear from how it
! is written.
!

subroutine munkres(C,n,jSol)
   ! Implementation of the Munkres algorithm (also referred to as the Hungarian
   ! algorithm). J. Munkres, Journal of the SIAM 5, 1 (1957)
   ! The following implementation is based on
   ! http://csclab.murraystate.edu/%7Ebob.pilgrim/445/munkres.html
   implicit none

   ! norm for matrix representation C(j,i) : j = columns, i = rows
   integer,intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
   real(dp), dimension(n,n), intent(in) :: C   ! cost matrix
   real(dp), dimension(n,n) :: CC   ! cost matrix

   real(dp) :: minCC
   real(dp) :: thresholdCC = 1.0d-9

   integer, dimension(n), intent(out) :: jSol ! solution indices

   ! following variables are use only for Munkres algo
   integer, dimension(n,n) :: M    ! mask matrix
   integer, dimension(n) :: rowCover, colCover !cover row and cols
   integer :: pathRow0 = 0, pathCol0 = 0  ! starting point for path finding part

   integer :: step, i, j, tmp
   logical :: done

   integer :: colCount
   logical :: exit_flag, starInRow
   integer :: row, col
   integer :: pathCount
   integer, dimension(2*n+1,2) :: path

   CC = C

   done = .false.
   step = 1
   tmp = 0

   do i = 1, n
      M(:,i) = 0
      rowCover(i) = 0
      colCover(i) = 0
   enddo

   do while(.not. done)
      select case(step)
      case(1)
 
          do i = 1, n
             minCC = CC(1,i)
             do j = 1, n
                if (CC(j,i) < minCC) minCC = CC(j,i)
             enddo
             CC(:,i) = CC(:,i) - minCC
          enddo

          step = 2

      case(2)

          do i = 1, n
             do j = 1, n
                if (CC(j,i) < thresholdCC .and. rowCover(i) == 0 .and. colCover(j) == 0) then
                   M(j,i) = 1
                   rowCover(i) = 1
                   colCover(j) = 1
                endif
             enddo
          enddo
          ! uncovers
          do i = 1, n
             rowCover(i) = 0
             colCover(i) = 0
          enddo

          step = 3

      case(3)

          colCount = 0
          do i = 1, n
             do j = 1, n
                ! if starred and column is uncovered
                if (M(j,i) == 1 .and. colCover(j) == 0) then
                   colCover(j) = 1
                   colCount = colCount + 1
                endif
             enddo
          enddo

          if (colCount == n) then
             step = 0
          else
             step = 4
          endif

      case(4)

          done = .false.

          do while (.not. done)
             ! find an uncovered zero
             starInRow = .false.

             row = 0
             col = 0
             exit_flag = .false.
             do i = 1, n
                do j = 1, n
                   if (CC(j,i) < thresholdCC .and. rowCover(i) == 0 .and. colCover(j) == 0) then
                      row = i
                      col = j
                      exit_flag = .true.
                      exit
                   endif
                enddo
                if (exit_flag) exit
             enddo
             if (row == 0) then !no zero uncoverred left
                done = .true.
                step = 6
             else
                M(col,row) = 2 !primed zero
                ! search if there is a starred zero in the same row
                do j = 1, n
                   if (M(j,row) == 1) then
                      starInRow = .true.
                      col = j
                   endif
                enddo
                if (starInRow) then ! if there is a starred zero in line
                   rowCover(row) = 1
                   colCover(col) = 0
                else ! if no starred zero in line
                   done = .true.
                   step = 5
                   pathRow0 = row
                   pathCol0 = col
                endif
             endif
          enddo

          done = .false.

      case(5)

          pathCount = 1

          path(pathCount,1) = pathRow0
          path(pathCount,2) = pathCol0

          done = .false.

          do while (.not. done)
             ! search for a starred zero in column
             row = 0
             col = path(pathCount,2)
             do i = 1, n
                if (M(col,i) == 1) row = i
             enddo
             if (row /= 0) then ! update path
                pathCount = pathCount + 1
                path(pathCount,1) = row
                path(pathCount,2) = path(pathCount-1,2)
             else
                done = .true.
             endif
             if (.not. done) then
                ! search for a prime zero in row
                do j = 1, n
                   if (M(j,row) == 2) col = j
                enddo
                ! update path
                pathCount = pathCount + 1
                path(pathCount,1) = path(pathCount-1,1)
                path(pathCount,2) = col
             endif
          enddo

          ! augment path
          do i = 1, pathCount
             if(M(path(i,2),path(i,1)) == 1) then
                M(path(i,2),path(i,1)) = 0
             else
                M(path(i,2),path(i,1)) = 1
             endif
          enddo

          ! clear covers and erase primes
          do i = 1, n
             rowCover(i) = 0
             colCover(i) = 0
             do j = 1, n
                if (M(j,i) == 2) M(j,i) = 0
             enddo
          enddo

          step = 3

          done = .false.

      case(6)

          minCC = huge(CC(1,1))

          do i = 1, n
             do j = 1, n
                if (rowCover(i) == 0 .and. colCover(j) == 0 .and. CC(j,i) < minCC) then
                   minCC = CC(j,i)
                endif
             enddo
          enddo
          do i = 1, n
             do j = 1, n
                if (rowCover(i) == 1) CC(j,i) = CC(j,i) + minCC
                if (colCover(j) == 0) CC(j,i) = CC(j,i) - minCC
             enddo
          enddo

          step = 4

      case default ! done
         do i = 1, n
            do j = 1, n
               if (M(j,i) == 1) jSol(i) = j
            enddo
         enddo
         done = .true.
      end select
   enddo

end subroutine munkres


! This returns the next combination for
! indexes of N choose M. If there are no
! more combinations then exit_flag will
! return .true. This assumes the first
! combination is (1,2,...,m).
subroutine nextCombination(m,n,c,exit_flag)
implicit none
integer,intent(in) :: m, n
integer,dimension(m),intent(inout) :: c

integer :: i, j, k
integer :: maxIndex

logical :: exit_flag

exit_flag = .false.

maxIndex = n
do i = m, 1, -1
    if (c(i) < maxIndex) then
        exit_flag = .true.

        k = c(i) + 1
        c(i) = k
        do j = i+1, m
            k = k + 1
            c(j) = k
        end do
    end if

    maxIndex = maxIndex - 1
    if (exit_flag) exit
end do

exit_flag = (.not. (exit_flag))

return
end subroutine nextCombination

subroutine chaseTwiddle_test()
implicit none
integer,parameter :: m = 4
integer,parameter :: n = 8
integer,dimension(n+2) :: p
integer :: x, y, z

integer :: i, j, k
logical :: exit_flag

call chaseInitTwiddle(m,n,p)

i = 1
do
    call chaseTwiddle(m,n,p,x,y,z,exit_flag)
    print *, i
    if (exit_flag) exit
    i = i + 1
end do

return
end subroutine chaseTwiddle_test

subroutine chaseInitTwiddle(m,n,p)
implicit none
integer,intent(in) :: m, n
integer,dimension(n+2),intent(out) :: p

integer :: i, j, k

p(1) = n+2 ! n+1
p(2:(n-m)+1) = 1 ! 0
i = (n-m)+2
do j = 1, m
    p(i) = j+1 !j
    i = i + 1
end do
p(n+2) = -1 !-2

if (m == 0) p(1) = 2 !1

return
end subroutine chaseInitTwiddle

subroutine chaseTwiddle(m,n,p,x,y,z,exit_flag)
implicit none
integer,intent(in) :: m, n
integer,dimension(n+2),intent(inout) :: p
integer,intent(out) :: x, y, z
logical,intent(out) :: exit_flag

integer :: i, j, k

exit_flag = .false.

j = 2
do
    if (p(j) > 1) exit !if (p(j) > 0) exit
    j = j + 1
end do

if (p(j-1) == 1) then !if (p(j-1) == 0) then
    do i = j-1, 2, -1
        p(i) = 0 !p(i) = -1
    end do
    p(j) = 1 !p(j) = 0

    x = 1
    z = 1
    p(2) = 2 !p(2) = 1
    y = j

else
    if (j > 2) p(j-1) = 1 !if (j > 2) p(j-1) = 0
    do
        j = j + 1
        if (p(j) <= 1) exit !if (p(j) <= 0) exit
    end do

    k = j - 1
    i = j
    do 
        if (p(i) /= 1) exit !if (p(i) /= 0) exit
        p(i) = 0 !p(i) = -1
        i = i + 1
    end do

    if (p(i) == 0) then !if (p(i) == -1) then
        p(i) = p(k)
        z = p(k)
        x = i
        y = k
        p(k) = 0 !p(k) = -1
    else
        if (i == p(1)) then
            exit_flag = .true.
        else
            p(j) = p(i)
            z = p(i)
            p(i) = 1 !p(i) = 0
            x = j
            y = i
        end if
    end if

end if

return
end subroutine chaseTwiddle

end module FUNCTIONS
