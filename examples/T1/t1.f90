module t1common
    implicit none

    real(8), allocatable :: a(:,:)
    real(8), allocatable :: b(:), x0(:), xexact(:)
    integer :: neq

contains
    subroutine t1init()
        character(10) :: rep
        character(7) :: field
        character(19) :: symm
        integer :: rows, cols, nnz, nnzm, i
        integer, allocatable :: indx(:), jndx(:)
        real(8), allocatable :: rval(:)
        complex :: cval(1)
        integer :: ival(1)

        ! Read FS 183 6 matrix from Matrix Market collection
        open(11, file='fs_183_6.mtx', status='OLD')
        call mminfo(11, rep, field, symm, rows, cols, nnzm)
        if (rep /= 'coordinate' .or. field /= 'real' .or. symm /= 'general' &
            & .or. rows /= cols) &
            & stop 1

        allocate(indx(nnzm), jndx(nnzm), rval(nnzm))
        
        call mmread(11, rep, field, symm, rows, cols, nnz, nnzm, &
            & indx, jndx, ival, rval, cval)

        close(11)

        neq = rows
        allocate(a(neq,neq), b(neq), x0(neq), xexact(neq))
        ! convert to dense
        a(:,:) = 0.d0
        do i = 1, nnz
            a(indx(i), jndx(i)) = rval(i)
        end do
            
        ! RHS vector according to 10.1137/070707373

        xexact(:) = 1.d0
        call dgemv('N', neq, neq, 1.d0, a, neq, xexact, 1, 0.d0, b, 1)
        x0(:) = 0.d0
    end subroutine t1init

    subroutine t1f(n, xcur, fcur, rpar, ipar, itrmf)
        integer, intent(in) :: n
        real(8), intent(in) :: xcur(n)
        real(8), intent(out) :: fcur(n)
        real(8), intent(inout) :: rpar(*)
        integer, intent(inout) :: ipar(*), itrmf
        
        call dcopy(n, b, 1, fcur, 1)
        call dgemv('N', n, n, 1.d0, a, n, xcur, 1, -1.d0, fcur, 1)
        itrmf = 0
    end subroutine t1f

    subroutine t1jv(n, xcur, fcur, ijob, v, z, rpar, ipar, iinf, riinf, itrmjv)
        integer, intent(in) :: n, ijob, iinf(3)
        real(8), intent(in) :: xcur(n), fcur(n), v(n), riinf(2)
        real(8), intent(out) :: z(n)

        real(8), intent(inout) :: rpar(*)
        integer, intent(inout) :: ipar(*), itrmjv

        if (ijob == 1) then
            ! no preconditioner
            itrmjv = 2
            return
        end if
        call dgemv('N', n, n, 1.d0, a, n, v, 1, 0.d0, z, 1)
        itrmjv = 0
    end subroutine t1jv
end module t1common

program t1
    use t1common
    implicit none

    real(8), allocatable :: x(:), rwork(:)
    real(8) :: rpar(1), ftol, rinp(9)
    integer :: info(6), ipar(1), inp(12), lrw, iterm

    real(8), parameter :: stptol = 0.d0
    !real(8), parameter :: ltol = 1.d-12
    !integer, parameter :: kdmax = 100
    !integer, parameter :: ikrysl = 3
    !integer, parameter :: iresup = 1
    integer :: ikrysl, iresup, kdmax
    real(8) :: ltol

    real(8), external :: ddot, dnrm2

    call t1init

    allocate(x(neq))

    call dcopy(neq, x0, 1, x, 1)

    write (*,*) 'Please enter: tol kdmax ikrysl iresup'
    read (*,*) ltol, kdmax, ikrysl, iresup
    if (ikrysl == 1 .or. ikrysl == 2) then
        ! provide memory for BiCG methods
        kdmax = 9
    end if

    inp(1) = 0
    inp(2) = 1
    inp(3) = ikrysl
    inp(4) = kdmax
    inp(5) = 0
    inp(6) = 10000
    inp(7) = iresup
    inp(8) = 0
    inp(9) = -1
    inp(10) = 3
    inp(11) = 4
    inp(12) = 0

    call nitdflts(rinp)
    rinp(6) = ltol

    lrw = neq*(kdmax+5)+kdmax*(kdmax+6)+1
    allocate(rwork(lrw))

    call t1f(neq, x, rwork, rpar, ipar, iterm)

    ftol = (1.d0+ltol)*ltol*dnrm2(neq, rwork, 1)

    call nitsol(neq, x, t1f, t1jv, ftol, stptol, inp, rinp, info, &
        & rwork, rpar, ipar, iterm, ddot, dnrm2)

    write (*,*) 'iterm: ', iterm
    ! check the solution against the exact one
    call daxpy(neq, -1.d0, xexact, 1, x, 1)
    write (*,*) 'Error: ', dnrm2(neq, x, 1)/dnrm2(neq, xexact, 1)

end program t1
