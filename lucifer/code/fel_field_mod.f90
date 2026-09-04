!+
! Module fel_field_mod
!
! The field set's types: one radiation field of the walk's set with its own slippage
! state and escape bank. They sit in a module of their own so that both the tracker
! (fel_track_mod) and the device seam (fel_device_mod) can take the set as an array
! of fel_field_struct, since passing a component section such as ffield(:)%wf, a
! wavefront with allocatable arrays inside, makes the compiler pack a temporary and
! deep-free it on return. fel_track_mod re-exports these types, so its users see them
! as before.
!-

module fel_field_mod

use wavefront_mod

implicit none

!+
! Struct fel_slip_struct
!
! The slippage state of one field record: Genesis's Field::first and Field::accuslip,
! plus the two run facts they are meaningless without. One per wavefront. The default is
! the steady state: timerun false makes fel_apply_slippage a no-op (as Genesis's)
! and first = 0 makes fel_field_index the identity.
!-

type fel_slip_struct
  logical :: timerun = .false.  ! Time-dependent run? False: slippage is a no-op.
  integer :: first = 0          ! Rotation offset of the field record, 0-based (Field::first).
  real(rp) :: accuslip = 0      ! Accumulated slippage [radiation wavelengths] (Field::accuslip).
  real(rp) :: sample = 1        ! Slice spacing / radiation wavelength (Control::sample).
  real(rp) :: u_escaped = 0     ! Energy [J] transmitted out of the window by slippage
                                ! (summed over every zero-filled slice). The TD energy
                                ! ledger's escape column, so E_beam + U_window + U_escaped
                                ! closes in a wake-free run.
end type

!+
! Struct fel_bank_struct
!
! Scratch carrier for the field slices one fel_apply_slippage call transmits out of the
! window: the caller passes it when it wants the light itself, not just its banked
! energy (keep_escaped_field). Reset and refilled per call. The caller drains it
! immediately (streams to file), so peak memory is a handful of grid planes -- one per
! rotation of the call, ~1 inside undulators, ~10 over an interlude.
!-

type fel_bank_struct
  complex(wf_rp), allocatable :: plane(:,:,:)  ! Transmitted planes, in transmission order.
  complex(wf_rp), allocatable :: plane_y(:,:,:)  ! Ey planes, allocated when Ey is live.
  integer :: n = 0                             ! How many this call transmitted.
end type

!+
! Struct fel_field_struct
!
! One radiation field of the walk's field set: the harmonic number, the wavefront
! record, and that field's own slippage state and escape bank. The walk carries an
! ordered set of these (Genesis's vector<Field*>), with the fundamental always
! entry 1. The ponderomotive phase, the phi0 advance and the slippage schedule are
! all defined against the fundamental. A harmonic field couples through fc(h) and the
! phase h*theta and diffracts at its own wavelength. A single-entry set is the
! pre-harmonic walk, bit for bit. wf%wavelength = fundamental wavelength / harm. Every
! field shares the time window, so the slippage state advances in fundamental-wavelength
! units for all of them (Genesis's one Control::sample) and the records rotate in
! lockstep.
!-

type fel_field_struct
  integer :: harm = 1               ! Harmonic number h (Field::harm).
  type (wavefront_struct) :: wf
  type (fel_slip_struct) :: slip
  type (fel_bank_struct) :: bank
end type

end module fel_field_mod
