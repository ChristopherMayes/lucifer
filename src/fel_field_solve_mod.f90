!+
! Module fel_field_solve_mod
!
! The passes of the averaged field solve over one plane of the field: the multiply by
! the step's propagator, the division by the transform's points and the addition of the
! source. fel_field_step (fel_track_mod) calls them between its transforms.
!
! The rounding of every operation here is fixed by the source. Each complex product is
! written out in its real and imaginary parts, a product is rounded before it is added
! or subtracted, and this file is compiled with -ffp-contract=off (CMakeLists.txt), so
! the compiler forms no fused multiply-add from them. The same statements written as
! complex array expressions in a file compiled at -O3 with contraction left on do not
! have one rounding: the compiler multiplies two complex arrays with separate
! instructions or with fused ones as the statements around them change, and the two
! differ in the last bit. A run's digits then moved when any line of the routine that
! held those expressions was edited. With the rounding written here, they do not.
!
! The arrays are explicit shape and carry no pointer or target attribute, so the
! compiler takes them as distinct and each loop is one pass over memory. The caller
! passes a plane of n by n points as n*n elements.
!-

module fel_field_solve_mod

use precision_def

implicit none

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_solve_copy (n, w, e)
!
! Routine to copy a plane into the work buffer, w = e.
!
! Input:
!   n      -- integer: The plane's points.
!   e(n)   -- complex(rp): The plane.
!
! Output:
!   w(n)   -- complex(rp): Its copy.
!-

subroutine fel_solve_copy (n, w, e)

integer n, i
complex(rp) w(n), e(n)

!

do i = 1, n
  w(i) = e(i)
enddo

end subroutine fel_solve_copy

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_solve_filter (n, w, sigmoid)
!
! Routine to multiply a transform by the source filter's sigmoid, w = w * sigmoid.
!
! Input:
!   n           -- integer: The plane's points.
!   w(n)        -- complex(rp): A transform, the source's or under the filter's
!                    self-test the field's.
!   sigmoid(n)  -- real(rp): The sigmoid, in the transform's order.
!
! Output:
!   w(n)        -- complex(rp): The product.
!-

subroutine fel_solve_filter (n, w, sigmoid)

integer n, i
complex(rp) w(n)
real(rp) sigmoid(n)

!

do i = 1, n
  w(i) = cmplx(real(w(i), rp) * sigmoid(i), aimag(w(i)) * sigmoid(i), rp)
enddo

end subroutine fel_solve_filter

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_solve_propagate (n, w, k)
!
! Routine to multiply a transformed plane by the step's propagator, w = w * k.
!
! Input:
!   n      -- integer: The plane's points.
!   w(n)   -- complex(rp): The plane's transform.
!   k(n)   -- complex(rp): The propagator.
!
! Output:
!   w(n)   -- complex(rp): The product.
!-

subroutine fel_solve_propagate (n, w, k)

integer n, i
complex(rp) w(n), k(n)
real(rp) wr, wi

!

do i = 1, n
  wr = real(w(i), rp);  wi = aimag(w(i))
  w(i) = cmplx(wr * real(k(i), rp) - wi * aimag(k(i)), wr * aimag(k(i)) + wi * real(k(i), rp), rp)
enddo

end subroutine fel_solve_propagate

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_solve_propagate_add (n, w, k, s, p)
!
! Routine to multiply a transformed plane by the step's propagator and add the
! transformed source, w = w * k + p * s. With the source filter on the source is added
! in Fourier space, before the one inverse transform.
!
! Input:
!   n      -- integer: The plane's points.
!   w(n)   -- complex(rp): The plane's transform.
!   k(n)   -- complex(rp): The propagator.
!   s(n)   -- complex(rp): The source's transform.
!   p      -- complex(rp): The factor on the source: 2 for a field of one polarization,
!               twice the plane's polarization component for a field of two.
!
! Output:
!   w(n)   -- complex(rp): The plane's transform after the step.
!-

subroutine fel_solve_propagate_add (n, w, k, s, p)

integer n, i
complex(rp) w(n), k(n), s(n), p
real(rp) wr, wi, sr, si, pr, pi_

!

pr = real(p, rp);  pi_ = aimag(p)

do i = 1, n
  wr = real(w(i), rp);  wi = aimag(w(i))
  sr = real(s(i), rp);  si = aimag(s(i))
  w(i) = cmplx((wr * real(k(i), rp) - wi * aimag(k(i))) + (pr * sr - pi_ * si), &
               (wr * aimag(k(i)) + wi * real(k(i), rp)) + (pr * si + pi_ * sr), rp)
enddo

end subroutine fel_solve_propagate_add

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_solve_normalize (n, e, w, rn)
!
! Routine to write a plane back from its inverse transform, e = w / rn, rn the
! transform's points, which an unnormalized transform pair multiplies by.
!
! Input:
!   n      -- integer: The plane's points.
!   w(n)   -- complex(rp): The inverse transform.
!   rn     -- real(rp): The transform's points.
!
! Output:
!   e(n)   -- complex(rp): The plane.
!-

subroutine fel_solve_normalize (n, e, w, rn)

integer n, i
complex(rp) e(n), w(n)
real(rp) rn

!

do i = 1, n
  e(i) = cmplx(real(w(i), rp) / rn, aimag(w(i)) / rn, rp)
enddo

end subroutine fel_solve_normalize

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!+
! Subroutine fel_solve_normalize_add (n, e, w, rn, s, p)
!
! Routine to write a plane back from its inverse transform and add the source,
! e = w / rn + p * s. With the source filter off the source is added in real space.
!
! Input:
!   n      -- integer: The plane's points.
!   w(n)   -- complex(rp): The inverse transform.
!   rn     -- real(rp): The transform's points.
!   s(n)   -- complex(rp): The source.
!   p      -- complex(rp): The factor on the source: 2 for a field of one polarization,
!               twice the plane's polarization component for a field of two.
!
! Output:
!   e(n)   -- complex(rp): The plane after the step.
!-

subroutine fel_solve_normalize_add (n, e, w, rn, s, p)

integer n, i
complex(rp) e(n), w(n), s(n), p
real(rp) rn, sr, si, pr, pi_

!

pr = real(p, rp);  pi_ = aimag(p)

do i = 1, n
  sr = real(s(i), rp);  si = aimag(s(i))
  e(i) = cmplx(real(w(i), rp) / rn + (pr * sr - pi_ * si), aimag(w(i)) / rn + (pr * si + pi_ * sr), rp)
enddo

end subroutine fel_solve_normalize_add

end module
