!+
! Module fel_field_mod
!
! Types for the field set. Each field of the set carries its own slippage state and
! escape bank.
!
! The types have a module of their own so fel_track_mod and fel_device_mod can both take
! the set as an array of fel_field_struct. Passing a component section instead, such as
! ffield(:)%wf, makes the compiler pack a temporary of allocatable arrays and deep-free
! it on return. fel_track_mod re-exports these types.
!-

module fel_field_mod

use wavefront_mod

implicit none

!+
! Struct fel_slip_struct
!
! Slippage state of one field record: Genesis's Field::first and Field::accuslip, and
! the two run facts they need. One per wavefront.
!
! The default is steady state. With timerun false, fel_apply_slippage does nothing, as
! in Genesis, and with first = 0, fel_field_index is the identity.
!-

type fel_slip_struct
  logical :: timerun = .false.  ! Time-dependent run? False: slippage is a no-op.
  integer :: first = 0          ! Rotation offset of the field record, 0-based (Field::first).
  real(rp) :: accuslip = 0      ! Accumulated slippage [radiation wavelengths] (Field::accuslip).
  integer :: sample = 1         ! Wavelengths per slice (Control::sample), the beam's
                                !   n_wavelength. An integer by construction.
  real(rp) :: u_escaped = 0     ! Energy [J] slippage carried out of the window, summed
                                !   over every zero-filled slice. The ledger's escape
                                !   column: E_beam + U_window + U_escaped closes without
                                !   wakes.
end type

!+
! Struct fel_bank_struct
!
! Scratch carrier for the slices one fel_apply_slippage call carries out of the window.
! The caller passes it for the light itself (keep_escaped_field), and not only its
! energy. It is reset and refilled each call.
!
! The caller drains it at once to file, so it holds one grid plane per rotation of the
! call: about 1 inside an undulator and about 10 over an interlude.
!-

type fel_bank_struct
  complex(wf_rp), allocatable :: plane(:,:,:)  ! Transmitted planes, in transmission order.
  complex(wf_rp), allocatable :: plane_y(:,:,:)  ! Ey planes, allocated when Ey is live.
  integer :: n = 0                             ! How many this call transmitted.
end type

!+
! Struct fel_field_struct
!
! One radiation field of the walk's field set: its harmonic number, its wavefront and
! its own slippage state and escape bank. The walk carries an ordered set of these, as
! Genesis carries vector<Field*>, and the fundamental is always entry 1.
!
! The ponderomotive phase, the phi0 advance and the slippage schedule are all defined
! against the fundamental. A harmonic field couples through fc(h) and the phase
! h*theta, and diffracts at its own wavelength: wf%wavelength = fundamental / harm.
! Every field shares the time window. So the slippage state of every field advances in
! fundamental wavelengths, Genesis's one Control::sample, and the records rotate in
! lockstep.
!-

type fel_field_struct
  integer :: harm = 1               ! Harmonic number h (Field::harm).
  type (wavefront_struct) :: wf
  type (fel_slip_struct) :: slip
  type (fel_bank_struct) :: bank
end type

end module fel_field_mod
