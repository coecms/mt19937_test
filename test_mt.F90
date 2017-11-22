program test_mt
    use mt19937_mod, only: MT19937_RANDOM_SEED, MT19937_RANDOM_NUMBER
    USE iso_c_binding, ONLY : int32 => c_int32_t, real32 => c_float, real64 => c_double
    implicit none
    real(kind=real64) :: r(10000)
    integer :: u, ios

    open(newunit=u, file='data.txt', action='write', iostat=ios, &
        status='unknown')
    if (ios /= 0) stop "FAILED TO OPEN data.txt"

    CALL MT19937_RANDOM_SEED([1_int32])
    CALL MT19937_RANDOM_NUMBER(r)

    write(u, '(F8.5)') r

    close(u)

end program test_mt
