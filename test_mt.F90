program test_mt
    use mt19937_mod, only: MT19937_RANDOM_SEED, MT19937_RANDOM_NUMBER
    USE iso_c_binding, ONLY : int32 => c_int32_t, real32 => c_float, real64 => c_double
    implicit none
    real(kind=real64) :: r(10000), r2(10000)
    integer :: u, ios

    open(newunit=u, file='data.txt', action='read', iostat=ios, &
        status='old')
    if (ios /= 0) stop "FAILED TO OPEN data.txt"

    CALL MT19937_RANDOM_SEED([1])
    CALL MT19937_RANDOM_NUMBER(r)

    read(u, '(F15.12)') r2

    if ( all(abs(r - r2) < 0.00000000001_real64) ) then
        print*, "all identical"
    else
        print*, "Failures"
    end if

    close(u)

end program test_mt
