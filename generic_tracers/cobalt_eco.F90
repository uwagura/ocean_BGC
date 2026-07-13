!> The cobalt_eco module contains subroutines for the diel vertical migration (DVM)
!! component of the COBALT biogeochemical model. The DVM-specific computation is
!! activated via the do_dvm flag in generic_COBALT_nml, but DVM tracer fields are
!! always registered to maintain a consistent NUM_ZOO=5 data layout.
module COBALT_eco

  use cobalt_types
  use MOM_file_parser, only : get_param, param_file_type
  use g_tracer_utils, only : g_tracer_add, g_tracer_type
  use g_tracer_utils, only : g_send_data, g_tracer_get_common, g_tracer_get_pointer
  use g_tracer_utils, only : g_diag_type
  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use time_manager_mod, only: time_type

  implicit none; private

  public cobalt_eco_add_tracers
  public eco_cobalt_reg_diag
  public eco_cobalt_send_diag
  public dvm_add_params
  public dvm_alloc_arrays
  public dvm_dealloc_arrays
  public dvm_migration

  contains

  !> Read all Diel Vertical Migration (DVM) namelist parameters into the zooplankton
  !! and COBALT derived types. Extracted from generic_COBALT user_add_params (Stage 2);
  !! called only when do_dvm is true.
  subroutine dvm_add_params(param_file, zoo, cobalt)
    type(param_file_type),                   intent(in)    :: param_file
    type(zooplankton), dimension(NUM_ZOO),   intent(inout) :: zoo
    type(generic_COBALT_type),               intent(inout) :: cobalt

    call get_param(param_file, "generic_COBALT", "q_p_2_n_vmmdz", zoo(4)%q_p_2_n, "Medium migrating zooplankton P:N", &
                   units="mol P mol N-1", default= 1.0/18.0)
    call get_param(param_file, "generic_COBALT", "q_p_2_n_vmlgz", zoo(5)%q_p_2_n, "Large migrating zooplankton P:N", &
                   units="mol P mol N-1", default= 1.0/16.0)
   call get_param(param_file, "generic_COBALT", "imax_vmmdz", zoo(4)%imax, &
                   "max ingestion rate for medium migrating zooplankton @ 0 deg. C", units="day-1", default=0.57, scale=I_sperd)
    call get_param(param_file, "generic_COBALT", "imax_vmlgz", zoo(5)%imax, &
                   "max ingestion rate for large migrating zooplankton @ 0 deg. C", units="day-1", default= 0.23, scale=I_sperd)
    call get_param(param_file, "generic_COBALT", "ki_vmmdz", zoo(4)%ki, "half-sat for ingestion by medium migrating zooplankton", &
                   units="mol N kg-1", default=1.25e-6)
    call get_param(param_file, "generic_COBALT", "ki_vmlgz", zoo(5)%ki, "half-sat for ingestion by large migrating zooplankton", &
                   units="mol N kg-1", default=1.25e-6)
    call get_param(param_file, "generic_COBALT", "ktemp_vmmdz", zoo(4)%ktemp, &
                   "exponential temperature dependence of medium migrating zooplankton rates", units="deg. C-1", default=0.063)
    call get_param(param_file, "generic_COBALT", "ktemp_vmlgz", zoo(5)%ktemp, &
                   "exponential temperature dependence of large migrating zooplankton rates", units="deg. C-1", default=0.063)
    call get_param(param_file, "generic_COBALT", "dvm_I_thresh_smz", zoo(1)%dvm_I_thresh, "Irradiance threshold for small zooplankton DVM", &
                   units="W m-2", default=0.0001)
    call get_param(param_file, "generic_COBALT", "dvm_I_thresh_mdz", zoo(2)%dvm_I_thresh, "Irradiance threshold for medium zooplankton DVM", &
                   units="W m-2", default=0.0001)
    call get_param(param_file, "generic_COBALT", "dvm_I_thresh_lgz", zoo(3)%dvm_I_thresh, "Irradiance threshold for large zooplankton DVM", &
                   units="W m-2", default=0.0001)
    call get_param(param_file, "generic_COBALT", "dvm_I_thresh_vmmdz", zoo(4)%dvm_I_thresh, "Irradiance threshold for medium migrating zooplankton DVM", &
                   units="W m-2", default=0.0001)
    call get_param(param_file, "generic_COBALT", "dvm_I_thresh_vmlgz", zoo(5)%dvm_I_thresh, "Irradiance threshold for large migrating zooplankton DVM", &
                   units="W m-2", default=0.0001)
    call get_param(param_file, "generic_COBALT", "swim_max_smz", zoo(1)%swim_max, "Maximum swimming speed for small zooplankton", &
                   units="m s-1", default=0.0)
    call get_param(param_file, "generic_COBALT", "swim_max_mdz", zoo(2)%swim_max, "Maximum swimming speed for medium zooplankton", &
                   units="m s-1", default=0.0)
    call get_param(param_file, "generic_COBALT", "swim_max_lgz", zoo(3)%swim_max, "Maximum swimming speed for large zooplankton", &
                   units="m s-1", default=0.0)
    call get_param(param_file, "generic_COBALT", "swim_max_vmmdz", zoo(4)%swim_max, "Maximum swimming speed for medium migrating zooplankton", &
                   units="m s-1", default=0.0) ! 0.06
    call get_param(param_file, "generic_COBALT", "swim_max_vmlgz", zoo(5)%swim_max, "Maximum swimming speed for large migrating zooplankton", &
                   units="m s-1", default=0.0) ! 0.08
    call get_param(param_file, "generic_COBALT", "swim_ref_smz", zoo(1)%swim_ref, "Reference swimming speed for small zooplankton", &
                   units="m s-1", default=0.0)
    call get_param(param_file, "generic_COBALT", "swim_ref_mdz", zoo(2)%swim_ref, "Reference swimming speed for medium zooplankton", &
                   units="m s-1", default=0.0)
    call get_param(param_file, "generic_COBALT", "swim_ref_lgz", zoo(3)%swim_ref, "Reference swimming speed for large zooplankton", &
                   units="m s-1", default=0.0)
    call get_param(param_file, "generic_COBALT", "swim_ref_vmmdz", zoo(4)%swim_ref, "Reference swimming speed for medium migrating zooplankton", &
                   units="m s-1", default=0.12)  ! 0.12
    call get_param(param_file, "generic_COBALT", "swim_ref_vmlgz", zoo(5)%swim_ref, "Reference swimming speed for large migrating zooplankton", &
                   units="m s-1", default=0.16)  ! 0.16
    call get_param(param_file, "generic_COBALT", "k_I_dvm_smz", zoo(1)%k_I_dvm, "Half-saturation irradiance for small zooplankton", &
                   units="W m-2", default=0.1)
    call get_param(param_file, "generic_COBALT", "k_I_dvm_mdz", zoo(2)%k_I_dvm, "Half-saturation irradiance for medium zooplankton", &
                   units="W m-2", default=0.1)
    call get_param(param_file, "generic_COBALT", "k_I_dvm_lgz", zoo(3)%k_I_dvm, "Half-saturation irradiance for large zooplankton", & 
                   units="W m-2", default=0.1)
    call get_param(param_file, "generic_COBALT", "k_I_dvm_vmmdz", zoo(4)%k_I_dvm, "Half-saturation irradiance for medium migrating zooplankton", &
                   units="W m-2", default=0.1)
    call get_param(param_file, "generic_COBALT", "k_I_dvm_vmlgz", zoo(5)%k_I_dvm, "Half-saturation irradiance for large migrating zooplankton", &
                   units="W m-2", default=0.1)
    call get_param(param_file, "generic_COBALT", "swim_stop_o2_smz", zoo(1)%swim_stop_o2, "Oxygen level to stop swimming for small zooplankton", &
                   units="mol O2", default=60.0e-6)
    call get_param(param_file, "generic_COBALT", "swim_stop_o2_mdz", zoo(2)%swim_stop_o2, "Oxygen level to stop swimming for medium zooplankton", &
                   units="mol O2", default=60.0e-6)
    call get_param(param_file, "generic_COBALT", "swim_stop_o2_lgz", zoo(3)%swim_stop_o2, "Oxygen level to stop swimming for large zooplankton", & 
                   units="mol O2", default=60.0e-6)
    call get_param(param_file, "generic_COBALT", "swim_stop_o2_vmmdz", zoo(4)%swim_stop_o2, "Oxygen level to stop swimming for medium migrating zooplankton", &
                   units="mol O2", default=60.0e-6)
    call get_param(param_file, "generic_COBALT", "swim_stop_o2_vmlgz", zoo(5)%swim_stop_o2, "Oxygen level to stop swimming for large migrating zooplankton", &
                   units="mol O2", default=60.0e-6)
    call get_param(param_file, "generic_COBALT", "smz_ipa_vmmdz", zoo(1)%ipa_vmmdz, &
                   "innate availability of medium migrating zooplankton to small zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "smz_ipa_vmlgz", zoo(1)%ipa_vmlgz, &
                   "innate availability of large migrating zooplankton to small zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "mdz_ipa_vmmdz", zoo(2)%ipa_vmmdz, &
                   "innate availability of medium migrating zooplankton to medium zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "mdz_ipa_vmlgz", zoo(2)%ipa_vmlgz, &
                   "innate availability of large migrating zooplankton to medium zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_smp", zoo(4)%ipa_smp, &
                   "innate availability of small phytoplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.4)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_mdp", zoo(4)%ipa_mdp, &
                   "innate availability of medium phytoplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_lgp", zoo(4)%ipa_lgp, &
                   "innate availability of large phytoplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_diaz",zoo(4)%ipa_diaz, &
                   "innate availability of diazotrophs to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.75)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_smz", zoo(4)%ipa_smz, &
                   "innate availability of small zooplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_mdz", zoo(4)%ipa_mdz, &
                   "innate availability of medium zooplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_lgz", zoo(4)%ipa_lgz, &
                   "innate availability of large zooplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_vmmdz", zoo(4)%ipa_vmmdz, &
                   "innate availability of medium migrating zooplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_vmlgz", zoo(4)%ipa_vmlgz, &
                   "innate availability of large migrating zooplankton to medium migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_bact", zoo(4)%ipa_bact, &
                   "innate availability of bacteria to medium migrating zooplankton feeding (0-1)", units="none", default=0.0)
    call get_param(param_file, "generic_COBALT", "vmmdz_ipa_det", zoo(4)%ipa_det, &
                   "innate availability of detritus to medium migrating zooplankton feeding (0-1)", units="none", default=0.0)
       call get_param(param_file, "generic_COBALT", "lgz_ipa_vmmdz", zoo(3)%ipa_vmmdz, &
                   "innate availability of medium migrating zooplankton to large zooplankton feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "lgz_ipa_vmlgz", zoo(3)%ipa_vmlgz, &
                   "innate availability of large migrating zooplankton to large zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_smp", zoo(5)%ipa_smp, &
                   "innate availability of small phytoplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_mdp", zoo(5)%ipa_mdp, &
                   "innate availability of medium phytoplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=0.4)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_lgp", zoo(5)%ipa_lgp, &
                   "innate availability of large phytoplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_diaz",zoo(5)%ipa_diaz, &
                   "innate availability of diazotrophs to large migrating zooplankton feeding (0-1)", units="none", &
                   default=0.4)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_smz", zoo(5)%ipa_smz, &
                   "innate availability of small zooplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_mdz", zoo(5)%ipa_mdz, &
                   "innate availability of medium zooplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_lgz", zoo(5)%ipa_lgz, &
                   "innate availability of large zooplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_vmmdz", zoo(5)%ipa_vmmdz, &
                   "innate availability of medium migrating zooplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_vmlgz", zoo(5)%ipa_vmlgz, &
                   "innate availability of large migrating zooplankton to large migrating zooplankton feeding (0-1)", units="none", &
                   default=0.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_bact",zoo(5)%ipa_bact, &
                   "innate availability of bacteria to large migrating zooplankton feeding (0-1)", units="none", default=0.0)
    call get_param(param_file, "generic_COBALT", "vmlgz_ipa_det", zoo(5)%ipa_det, &
                   "innate availability of detritus to large migrating zooplankton feeding (0-1)", units="none", default=0.0)
    call get_param(param_file, "generic_COBALT", "nswitch_vmmdz", zoo(4)%nswitch, &
                   "prey switching parameter 1 for medium migrating zooplankton", units="none", default=2.0)
    call get_param(param_file, "generic_COBALT", "nswitch_vmlgz", zoo(5)%nswitch, &
                   "prey switching parameter 1 for large migrating zooplankton", units="none", default=2.0)
    call get_param(param_file, "generic_COBALT", "mswitch_vmmdz", zoo(4)%mswitch, &
                   "prey switching parameter 2 for medium migrating zooplankton", units="none", default=2.0)
    call get_param(param_file, "generic_COBALT", "mswitch_vmlgz", zoo(5)%mswitch, &
                   "prey switching parameter 2 for large migrating zooplankton", units="none", default=2.0)
    call get_param(param_file, "generic_COBALT", "gge_max_vmmdz",zoo(4)%gge_max, &
                   "maximum gross growth efficiency for medium migrating zooplankton", units="none", default=0.4)
    call get_param(param_file, "generic_COBALT", "gge_max_vmlgz",zoo(5)%gge_max, &
                   "maximum gross growth efficiency for large migrating zooplankton", units="none", default=0.4)
    call get_param(param_file, "generic_COBALT", "bresp_vmmdz", zoo(4)%bresp, &
                   "basal respiration rate for medium migrating zooplankton", units="day-1", default=0.008,scale=I_sperd)
    call get_param(param_file, "generic_COBALT", "bresp_vmlgz", zoo(5)%bresp, &
                   "basal respiration rate for large migrating zooplankton", units="day-1", default=0.0032, scale=I_sperd)
    call get_param(param_file, "generic_COBALT", "phi_det_vmmdz", zoo(4)%phi_det, &
                   "fraction of ingestion by medium migrating zooplankton to detritus", units="none", default=0.15)
    call get_param(param_file, "generic_COBALT", "phi_det_vmlgz", zoo(5)%phi_det, &
                   "fraction of ingestion by large migrating zooplankton to detritus", units="none", default=0.30)
    call get_param(param_file, "generic_COBALT", "phi_ldon_vmmdz", zoo(4)%phi_ldon, &
                   "fraction of N ingestion by medium migrating zooplankton to labile dissolved organic nitrogen", &
                   units="none", default=0.625*(0.30-zoo(4)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_ldon_vmlgz", zoo(5)%phi_ldon, &
                   "fraction of N ingestion by large migrating zooplankton to labile dissolved organic nitrogen", &
                   units="none", default=0.625*(0.30-zoo(5)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_ldop_vmmdz", zoo(4)%phi_ldop, &
                   "fraction of P ingestion by medium migrating zooplankton to labile dissolved organic phosphorus", &
                   units="none", default=0.575*(0.30-zoo(4)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_ldop_vmlgz", zoo(5)%phi_ldop, &
                   "fraction of P ingestion by large migrating zooplankton to labile dissolved organic phosphorus", &
                   units="none", default=0.575*(0.30-zoo(5)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_srdon_vmmdz", zoo(4)%phi_srdon, &
                   "fraction of N ingestion by medium migrating zooplankton to semi-refractory dissolved organic nitrogen", &
                   units="none", default=0.075*(0.30-zoo(4)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_srdon_vmlgz", zoo(5)%phi_srdon, &
                   "fraction of N ingestion by large migrating zooplankton to semi-refractory dissolved organic nitrogen", &
                   units="none", default=0.075*(0.30-zoo(5)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_srdop_vmmdz", zoo(4)%phi_srdop, &
                   "fraction of P ingestion by medium migrating zooplankton to semi-refractory dissolved organic phosphorus", &
                   units="none", default=0.125*(0.30-zoo(4)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_srdop_vmlgz", zoo(5)%phi_srdop, &
                   "fraction of P ingestion by large migrating zooplankton to semi-refractory dissolved organic phosphorus", &
                   units="none", default=0.125*(0.30-zoo(5)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_sldon_vmmdz", zoo(4)%phi_sldon, &
                   "fraction of N ingestion by medium migrating zooplankton to semi-labile dissolved organic nitrogen", &
                   units="none", default=0.3*(0.30-zoo(4)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_sldon_vmlgz", zoo(5)%phi_sldon, &
                   "fraction of N ingestion by large migrating zooplankton to semi-labile dissolved organic nitrogen", &
                   units="none", default=0.3*(0.30-zoo(5)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_sldop_vmmdz", zoo(4)%phi_sldop, &
                   "fraction of P ingestion by medium migrating zooplankton to semi-labile dissolved organic phosphorus", &
                   units="none", default=0.3*(0.30-zoo(4)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_sldop_vmlgz", zoo(5)%phi_sldop, &
                   "fraction of P ingestion by large migrating zooplankton to semi-labile dissolved organic phosphorus", &
                   units="none", default=0.3*(0.30-zoo(5)%phi_det))
    call get_param(param_file, "generic_COBALT", "phi_det_si_vmmdz", zoo(4)%phi_det_si, &
                   "fraction of silica ingestion by medium migrating zooplankton to si detritus", units="none", default=0.15)
    call get_param(param_file, "generic_COBALT", "phi_det_si_vmlgz", zoo(5)%phi_det_si, &
                   "fraction of silica ingestion by large migrating zooplankton to si detritus", units="none", default=0.30)
    call get_param(param_file, "generic_COBALT", "hp_ipa_vmmdz", cobalt%hp_ipa_vmmdz, &
                   "innate availability of medium migrating zooplankton to higher predator feeding (0-1)", units="none", &
                   default=1.0)
    call get_param(param_file, "generic_COBALT", "hp_ipa_vmlgz", cobalt%hp_ipa_vmlgz, &
                   "innate availability of large migrating zooplankton to higher predator feeding (0-1)", units="none", &
                   default=1.0)

  end subroutine dvm_add_params

  !> Allocate all Diel Vertical Migration (DVM) state and initialise it to zero:
  !! the cobalt-level dvm component (12 source/sink flux arrays; the tracer
  !! indices are scalars set later during tracer registration) and the zoo(4:5)
  !! gut/metabolite/migration fields. Extracted verbatim from generic_COBALT
  !! user_allocate_arrays (Stage 3); called only when do_dvm is true.
  !! NOTE: the jnvmlgz_met allocation re-zeroes cobalt%dvm%jnvmlgz rather than
  !! jnvmlgz_met; this reproduces an mpoupon quirk verbatim to preserve
  !! bit-for-bit answers (jnvmlgz_met is fully written before it is read).
  subroutine dvm_alloc_arrays(cobalt, zoo, isd, ied, jsd, jed, nk)
    type(generic_COBALT_type),             intent(inout) :: cobalt
    type(zooplankton), dimension(NUM_ZOO), intent(inout) :: zoo
    integer,                               intent(in)    :: isd, ied, jsd, jed, nk

    integer :: n

    allocate(cobalt%dvm)
    allocate(cobalt%dvm%jnvmmdz(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jnvmmdz=0.0
    allocate(cobalt%dvm%jnvmlgz(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jnvmlgz=0.0
    allocate(cobalt%dvm%jnvmmdz_gut(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jnvmmdz_gut=0.0
    allocate(cobalt%dvm%jnvmlgz_gut(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jnvmlgz_gut=0.0
    allocate(cobalt%dvm%jpvmmdz_gut(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jpvmmdz_gut=0.0
    allocate(cobalt%dvm%jpvmlgz_gut(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jpvmlgz_gut=0.0
    allocate(cobalt%dvm%jfevmmdz_gut(isd:ied, jsd:jed, 1:nk))       ; cobalt%dvm%jfevmmdz_gut=0.0
    allocate(cobalt%dvm%jfevmlgz_gut(isd:ied, jsd:jed, 1:nk))       ; cobalt%dvm%jfevmlgz_gut=0.0
    allocate(cobalt%dvm%jsivmmdz_gut(isd:ied, jsd:jed, 1:nk))       ; cobalt%dvm%jsivmmdz_gut=0.0
    allocate(cobalt%dvm%jsivmlgz_gut(isd:ied, jsd:jed, 1:nk))       ; cobalt%dvm%jsivmlgz_gut=0.0
    allocate(cobalt%dvm%jnvmmdz_met(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jnvmmdz_met=0.0
    allocate(cobalt%dvm%jnvmlgz_met(isd:ied, jsd:jed, 1:nk))        ; cobalt%dvm%jnvmlgz=0.0
    do n = NUM_BASE_ZOO+1, NUM_ZOO
       allocate(zoo(n)%lim_nut_n_ingestion(isd:ied,jsd:jed,nk))   ; zoo(n)%lim_nut_n_ingestion   = 0.0
       allocate(zoo(n)%jmetabo_n(isd:ied,jsd:jed,nk))      ; zoo(n)%jmetabo_n      = 0.0
       allocate(zoo(n)%f_gut_n(isd:ied,jsd:jed,nk))        ; zoo(n)%f_gut_n        = 0.0
       allocate(zoo(n)%f_gut_p(isd:ied,jsd:jed,nk))        ; zoo(n)%f_gut_p        = 0.0
       allocate(zoo(n)%f_gut_fe(isd:ied,jsd:jed,nk))       ; zoo(n)%f_gut_fe       = 0.0
       allocate(zoo(n)%f_gut_si(isd:ied,jsd:jed,nk))       ; zoo(n)%f_gut_si       = 0.0
       allocate(zoo(n)%f_met_n(isd:ied,jsd:jed,nk))        ; zoo(n)%f_met_n        = 0.0
       allocate(zoo(n)%jclear_gut_n(isd:ied,jsd:jed,nk))   ; zoo(n)%jclear_gut_n   = 0.0
       allocate(zoo(n)%jprod_gut_n(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_gut_n    = 0.0
       allocate(zoo(n)%jclear_gut_p(isd:ied,jsd:jed,nk))   ; zoo(n)%jclear_gut_p   = 0.0
       allocate(zoo(n)%jprod_gut_p(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_gut_p    = 0.0
       allocate(zoo(n)%jclear_gut_fe(isd:ied,jsd:jed,nk))  ; zoo(n)%jclear_gut_fe  = 0.0
       allocate(zoo(n)%jprod_gut_fe(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_gut_fe   = 0.0
       allocate(zoo(n)%jclear_gut_si(isd:ied,jsd:jed,nk))  ; zoo(n)%jclear_gut_si  = 0.0
       allocate(zoo(n)%jprod_gut_si(isd:ied,jsd:jed,nk))   ; zoo(n)%jprod_gut_si   = 0.0
       allocate(zoo(n)%jclear_met_n(isd:ied,jsd:jed,nk))   ; zoo(n)%jclear_met_n   = 0.0
       allocate(zoo(n)%jprod_met_n(isd:ied,jsd:jed,nk))    ; zoo(n)%jprod_met_n    = 0.0
       allocate(zoo(n)%vmove_met(isd:ied,jsd:jed,nk))      ; zoo(n)%vmove_met      = 0.0
       allocate(zoo(n)%vmove_gut(isd:ied,jsd:jed,nk))      ; zoo(n)%vmove_gut      = 0.0
       allocate(zoo(n)%vmove_gut_p(isd:ied,jsd:jed,nk))    ; zoo(n)%vmove_gut_p    = 0.0
       allocate(zoo(n)%vmove_gut_fe(isd:ied,jsd:jed,nk))   ; zoo(n)%vmove_gut_fe   = 0.0
       allocate(zoo(n)%vmove_gut_si(isd:ied,jsd:jed,nk))   ; zoo(n)%vmove_gut_si   = 0.0
    enddo
  end subroutine dvm_alloc_arrays

  !> Deallocate all DVM state allocated by dvm_alloc_arrays. Mirror of that
  !! routine; called only when do_dvm is true.
  subroutine dvm_dealloc_arrays(cobalt, zoo)
    type(generic_COBALT_type),             intent(inout) :: cobalt
    type(zooplankton), dimension(NUM_ZOO), intent(inout) :: zoo

    integer :: n

    deallocate(cobalt%dvm%jnvmmdz)
    deallocate(cobalt%dvm%jnvmlgz)
    deallocate(cobalt%dvm%jnvmmdz_gut)
    deallocate(cobalt%dvm%jnvmlgz_gut)
    deallocate(cobalt%dvm%jpvmmdz_gut)
    deallocate(cobalt%dvm%jpvmlgz_gut)
    deallocate(cobalt%dvm%jfevmmdz_gut)
    deallocate(cobalt%dvm%jfevmlgz_gut)
    deallocate(cobalt%dvm%jsivmmdz_gut)
    deallocate(cobalt%dvm%jsivmlgz_gut)
    deallocate(cobalt%dvm%jnvmmdz_met)
    deallocate(cobalt%dvm%jnvmlgz_met)
    deallocate(cobalt%dvm)
    do n = NUM_BASE_ZOO+1, NUM_ZOO
       deallocate(zoo(n)%lim_nut_n_ingestion)
       deallocate(zoo(n)%jmetabo_n)
       deallocate(zoo(n)%f_gut_n)
       deallocate(zoo(n)%f_gut_p)
       deallocate(zoo(n)%f_gut_fe)
       deallocate(zoo(n)%f_gut_si)
       deallocate(zoo(n)%f_met_n)
       deallocate(zoo(n)%jclear_gut_n)
       deallocate(zoo(n)%jprod_gut_n)
       deallocate(zoo(n)%jclear_gut_p)
       deallocate(zoo(n)%jprod_gut_p)
       deallocate(zoo(n)%jclear_gut_fe)
       deallocate(zoo(n)%jprod_gut_fe)
       deallocate(zoo(n)%jclear_gut_si)
       deallocate(zoo(n)%jprod_gut_si)
       deallocate(zoo(n)%jclear_met_n)
       deallocate(zoo(n)%jprod_met_n)
       deallocate(zoo(n)%vmove_met)
       deallocate(zoo(n)%vmove_gut)
       deallocate(zoo(n)%vmove_gut_p)
       deallocate(zoo(n)%vmove_gut_fe)
       deallocate(zoo(n)%vmove_gut_si)
    enddo
  end subroutine dvm_dealloc_arrays

  !> Compute diel vertical migration (DVM) swimming velocities (vmove*) for the
  !! migrating zooplankton groups. Extracted verbatim from generic_COBALT
  !! update_from_source (Stage 4): the "smart migration" prep (prey-weighted
  !! normalised cumulative biomass) followed by the vmove assignment loop. In
  !! the original routine these two blocks were ~560 lines apart, but both
  !! depend only on read-only state (f_n, rho_dzt), so co-locating them here is
  !! arithmetically identical. The former update_from_source allocatable
  !! temporaries are now automatic locals. Called only when do_dvm is true; the
  !! subsequent g_tracer_set_values('vmove',...) calls stay in generic_COBALT.
  subroutine dvm_migration(zoo, phyto, cobalt, rho_dzt, ilb, jlb, isc, iec, jsc, jec, nk)
    type(zooplankton),   dimension(NUM_ZOO),   intent(inout) :: zoo
    type(phytoplankton), dimension(NUM_PHYTO), intent(in)    :: phyto
    type(generic_COBALT_type),                 intent(in)    :: cobalt
    real, dimension(ilb:,jlb:,:),              intent(in)    :: rho_dzt
    integer,                                   intent(in)    :: ilb, jlb, isc, iec, jsc, jec, nk

    ! smart-migration temporaries (were allocatables in update_from_source)
    real, dimension(isc:iec,jsc:jec,nk) :: vmmd_rho_dzt, vmlg_rho_dzt
    real, dimension(isc:iec,jsc:jec,nk) :: vmmd_prey_rho_dzt, vmlg_prey_rho_dzt
    real, dimension(isc:iec,jsc:jec,nk) :: vmmd_norm_cum, vmlg_norm_cum
    real, dimension(isc:iec,jsc:jec,nk) :: vmmd_prey_norm_cum, vmlg_prey_norm_cum
    real, dimension(isc:iec,jsc:jec)    :: vmmd_int, vmlg_int, vmmd_prey_int, vmlg_prey_int
    integer :: i, j, k, n
    real    :: swim

    ! --- smart-migration prep: prey-weighted normalised cumulative biomass ---
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; !{
        vmmd_rho_dzt(i,j,k) = rho_dzt(i,j,k) * zoo(4)%f_n(i,j,k)
        vmlg_rho_dzt(i,j,k) = rho_dzt(i,j,k) * zoo(5)%f_n(i,j,k)
        vmmd_prey_rho_dzt(i,j,k) =  rho_dzt(i,j,k) * ( phyto(DIAZO)%f_n(i,j,k) + \
                                    phyto(LARGE)%f_n(i,j,k) + zoo(1)%f_n(i,j,k) )
        vmlg_prey_rho_dzt(i,j,k) = rho_dzt(i,j,k) * ( phyto(DIAZO)%f_n(i,j,k) + \
                                   phyto(LARGE)%f_n(i,j,k) + zoo(2)%f_n(i,j,k) + zoo(4)%f_n(i,j,k) )
    enddo; enddo; enddo;  !}  i, j, k

    do j = jsc, jec ; do i = isc, iec; !{
        vmmd_int(i,j) = 0
        vmlg_int(i,j) = 0
        vmmd_prey_int(i,j) = 0
        vmlg_prey_int(i,j) = 0
    enddo; enddo;  !}  j, i

    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec; !{
        vmmd_int(i,j)  =  vmmd_int(i,j)  +   vmmd_rho_dzt(i,j,k)
        vmlg_int(i,j)  =  vmlg_int(i,j)  +   vmlg_rho_dzt(i,j,k) 
        vmmd_prey_int(i,j)  =  vmmd_prey_int(i,j)  +   vmmd_prey_rho_dzt(i,j,k) 
        vmlg_prey_int(i,j)  =  vmlg_prey_int(i,j)  +   vmlg_prey_rho_dzt(i,j,k)
    enddo; enddo; enddo;  !}  i, j, k 

    do j = jsc, jec ; do i = isc, iec; !{
        vmmd_norm_cum(i,j,1)  = vmmd_rho_dzt(i,j,1)    /   vmmd_int(i,j)
        vmlg_norm_cum(i,j,1)    = vmlg_rho_dzt(i,j,1)      /   vmlg_int(i,j)
        vmmd_prey_norm_cum(i,j,1)  = vmmd_prey_rho_dzt(i,j,1)  /   vmmd_prey_int(i,j)
        vmlg_prey_norm_cum(i,j,1)    = vmlg_prey_rho_dzt(i,j,1)    /   vmlg_prey_int(i,j)
    enddo; enddo;  !}  j, i

    do k = 2, nk ; do j = jsc, jec ; do i = isc, iec; !{
        ! Normalization and cumulative sum
        vmmd_norm_cum(i,j,k)  = vmmd_norm_cum(i,j,k-1) + vmmd_rho_dzt(i,j,k)    /   vmmd_int(i,j)
        vmlg_norm_cum(i,j,k)    = vmlg_norm_cum(i,j,k-1) + vmlg_rho_dzt(i,j,k)      /   vmlg_int(i,j)
        vmmd_prey_norm_cum(i,j,k)  =  vmmd_prey_norm_cum(i,j,k-1)  + vmmd_prey_rho_dzt(i,j,k)  /   vmmd_prey_int(i,j)
        vmlg_prey_norm_cum(i,j,k)    = vmlg_prey_norm_cum(i,j,k-1) + vmlg_prey_rho_dzt(i,j,k)    /   vmlg_prey_int(i,j)
    enddo; enddo; enddo;  !}  i, j, k

    !
    ! 3.2.4 Zooplankton migration: vmove assignment
    !
    do k = 1, nk ; do j = jsc, jec ; do i = isc, iec   !{
       do n = 2, NUM_ZOO !{
           
           swim = zoo(n)%swim_max * abs( LOG(zoo(n)%dvm_I_thresh/(epsln+cobalt%irr_inst(i,j,k))) / 0.0232) / &
                   (50.0 + abs( LOG(zoo(n)%dvm_I_thresh/(epsln+cobalt%irr_inst(i,j,k))) / 0.0232))
           
           ! Upward swimming (during night)
           if ( cobalt%irr_inst(i,j,k) .lt. zoo(n)%dvm_I_thresh ) then
              
                ! Medium migratory zooplankton (4)
                if ( n .eq. 4 .and. do_dvm ) then
                   if (vmmd_prey_norm_cum(i,j,k) .gt.  vmmd_norm_cum(i,j,k)) then
                     zoo(n)%vmove(i,j,k) = -swim ! Upward
                     zoo(n)%vmove_met(i,j,k) = -swim ! Upward  
                     zoo(n)%vmove_gut(i,j,k) = -swim ! Upward
                     zoo(n)%vmove_gut_p(i,j,k)  = -swim ! Upward
                     zoo(n)%vmove_gut_fe(i,j,k) = -swim ! Upward 
                     zoo(n)%vmove_gut_si(i,j,k) = -swim ! Upward   
                  else
                     zoo(n)%vmove(i,j,k) = swim  ! Downward
                     zoo(n)%vmove_met(i,j,k) = swim ! Downward  
                     zoo(n)%vmove_gut(i,j,k) = swim ! Downward     
                     zoo(n)%vmove_gut_p(i,j,k)  = swim ! Downward 
                     zoo(n)%vmove_gut_fe(i,j,k) = swim ! Downward
                     zoo(n)%vmove_gut_si(i,j,k) = swim ! Downward 
                  endif
             
                ! Large migratory zooplankton (5)
                else if ( n .eq. 5 .and. do_dvm ) then
                   if (vmlg_prey_norm_cum(i,j,k) .gt.  vmlg_norm_cum(i,j,k)) then
                     zoo(n)%vmove(i,j,k)        = -swim ! Upward
                     zoo(n)%vmove_met(i,j,k)    = -swim ! Upward  
                     zoo(n)%vmove_gut(i,j,k)    = -swim ! Upward
                     zoo(n)%vmove_gut_p(i,j,k)  = -swim ! Upward
                     zoo(n)%vmove_gut_fe(i,j,k) = -swim ! Upward
                     zoo(n)%vmove_gut_si(i,j,k) = -swim ! Upward
                  else
                     zoo(n)%vmove(i,j,k)        = swim ! Downward
                     zoo(n)%vmove_met(i,j,k)    = swim ! Downward
                     zoo(n)%vmove_gut(i,j,k)    = swim ! Downward
                     zoo(n)%vmove_gut_p(i,j,k)  = swim ! Downward
                     zoo(n)%vmove_gut_fe(i,j,k) = swim ! Downward
                     zoo(n)%vmove_gut_si(i,j,k) = swim ! Downward
                  endif   
             
               ! Others (1,2,3)
               else
                  zoo(n)%vmove(i,j,k) = -swim ! Upward
               endif

           ! Downward swimming (during day)
           else
                
                ! Medium migratory zooplankton (4)
                if ( n .eq. 4 .and. do_dvm ) then
                    ! Enough oxygen
                    if (cobalt%f_o2(i,j,k) .gt. zoo(n)%swim_stop_o2) then
                         zoo(n)%vmove(i,j,k) = swim        ! Downward
                         zoo(n)%vmove_met(i,j,k) = swim    ! Downward
                         zoo(n)%vmove_gut(i,j,k) = swim    ! Downward
                         zoo(n)%vmove_gut_p(i,j,k)  = swim ! Downward
                         zoo(n)%vmove_gut_fe(i,j,k) = swim ! Downward
                         zoo(n)%vmove_gut_si(i,j,k) = swim ! Downward
                    
                    ! Not enough oxygen
                    else
                         zoo(n)%vmove(i,j,k) = 0.0        ! No swimming
                         zoo(n)%vmove_met(i,j,k) = 0.0    ! No swimming
                         zoo(n)%vmove_gut(i,j,k) = 0.0    ! No swimming
                         zoo(n)%vmove_gut_p(i,j,k)  = 0.0 ! No swimming
                         zoo(n)%vmove_gut_fe(i,j,k) = 0.0 ! No swimming
                         zoo(n)%vmove_gut_si(i,j,k) = 0.0 ! No swimming
                    endif
                                 

                ! Large migratory zooplankton (5)
                else if ( n .eq. 5 .and. do_dvm ) then
                   ! Enough oxygen
                   if (cobalt%f_o2(i,j,k) .gt. zoo(n)%swim_stop_o2) then
                        zoo(n)%vmove(i,j,k) = swim        ! Downward
                        zoo(n)%vmove_met(i,j,k) = swim    ! Downward
                        zoo(n)%vmove_gut(i,j,k) = swim    ! Downward
                        zoo(n)%vmove_gut_p(i,j,k)  = swim ! Downward
                        zoo(n)%vmove_gut_fe(i,j,k) = swim ! Downward
                        zoo(n)%vmove_gut_si(i,j,k) = swim ! Downward
               
                   ! Not enough oxygen
                   else
                        zoo(n)%vmove(i,j,k) = 0.0        ! No swimming
                        zoo(n)%vmove_met(i,j,k) = 0.0    ! No swimming
                        zoo(n)%vmove_gut(i,j,k) = 0.0    ! No swimming
                        zoo(n)%vmove_gut_p(i,j,k)  = 0.0 ! No swimming
                        zoo(n)%vmove_gut_fe(i,j,k) = 0.0 ! No swimming
                        zoo(n)%vmove_gut_si(i,j,k) = 0.0 ! No swimming
                   endif

               ! Others (1,2,3)
               else
                  ! Enough oxygen
                  if (cobalt%f_o2(i,j,k) .gt. zoo(n)%swim_stop_o2) then 
                       zoo(n)%vmove(i,j,k) = swim ! Downward
              
                  ! Not enough oxygen
                  else
                       zoo(n)%vmove(i,j,k) = 0.0 ! No swimming
                  endif
               endif 
          endif
 
       enddo !} n
    enddo; enddo; enddo; !} i, j, k
  end subroutine dvm_migration


  !> Register DVM-specific prognostic tracers for migrating zooplankton groups.
  !! These tracers are always registered (NUM_ZOO=5 is fixed) but their DVM-specific
  !! computation (swim, gut clearance, bioenergetics) is gated by do_dvm.
  subroutine cobalt_eco_add_tracers(tracer_list, package_name)
    type(g_tracer_type), pointer :: tracer_list
    character(len=*), intent(in) :: package_name

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nvmmdz',            &
         longname   = 'medium migrating Zooplankton Nitrogen (Biomass)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nvmmdz_gut',            &
         longname   = 'medium migrating Zooplankton Nitrogen (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'pvmmdz_gut',            &
         longname   = 'medium migrating Zooplankton Phosphorus (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'sivmmdz_gut',            &
         longname   = 'medium migrating Zooplankton Silicon (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'fevmmdz_gut',            &
         longname   = 'medium migrating Zooplankton Iron (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nvmmdz_met',            &
         longname   = 'medium migrating Zooplankton Nitrogen (Metabolites)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nvmlgz',            &
         longname   = 'large migrating Zooplankton Nitrogen (Biomass)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nvmlgz_gut',            &
         longname   = 'large migrating Zooplankton Nitrogen (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'pvmlgz_gut',            &
         longname   = 'large migrating Zooplankton Phosphorus (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'sivmlgz_gut',            &
         longname   = 'large migrating Zooplankton Silicon (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'fevmlgz_gut',            &
         longname   = 'large migrating Zooplankton Iron (Gut)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

    call g_tracer_add(tracer_list,package_name,&
         name       = 'nvmlgz_met',            &
         longname   = 'large migrating Zooplankton Nitrogen (Metabolites)', &
         units      = 'mol/kg',         &
         prog       = .true.,           &
         move_vertical = .true.)

  end subroutine cobalt_eco_add_tracers

  !> Register DVM-specific diagnostic fields for migrating zooplankton groups (zoo(4) and zoo(5)).
  !! Called from cobalt_reg_diagnostics only when do_dvm is true.
  subroutine eco_cobalt_reg_diag(diag_list, axes, init_time, zoo)
    type(g_diag_type), pointer                         :: diag_list
    integer,                             intent(in)    :: axes(3)
    type(time_type),                     intent(in)    :: init_time
    type(zooplankton), dimension(NUM_ZOO), intent(inout) :: zoo

    ! local
    type(vardesc) :: vardesc_temp

    vardesc_temp = vardesc("jzloss_n_vmMdz","Medium-sized migrating zooplankton nitrogen loss to zooplankton",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jzloss_n_vmLgz","Large migratingzooplankton nitrogen loss to zooplankton",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jzloss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jzloss_p_vmMdz","Medium-sized migrating zooplankton phosphorus loss to zooplankton",&
                           'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jzloss_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jzloss_p_vmLgz","Large migrating zooplankton phosphorus loss to zooplankton",&
                           'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jzloss_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jhploss_n_vmMdz","Medium-sized migrating zooplankton nitrogen loss to higher predators",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jhploss_n_vmLgz","Large migrating zooplankton nitrogen loss to higher predators",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jhploss_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jhploss_p_vmMdz","Medium-sized migrating zooplankton phosphorus loss to higher predators",&
                           'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jhploss_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jhploss_p_vmLgz","Large migrating zooplankton phosphorus loss to higher predators",&
                           'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jhploss_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
     vardesc_temp = vardesc("jingest_n_vmMdz","Ingestion of nitrogen by medium-sized migrating zooplankton",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_n_vmLgz","Ingestion of nitrogen by large migrating zooplankton",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jingest_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_p_vmMdz","Ingestion of phosphorous by medium-sized migrating zooplankton",&
                           'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_p_vmLgz","Ingestion of phosphorous by large migrating zooplankton",&
                           'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jingest_p = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
     vardesc_temp = vardesc("jingest_sio2_vmMdz","Ingestion of sio2 by medium-sized migrating zooplankton",&
                           'h','L','s','mol SiO2 kg-1 s-1','f')
    zoo(4)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_sio2_vmLgz","Ingestion of sio2 by large migrating zooplankton",&
                           'h','L','s','mol SiO2 kg-1 s-1','f')
    zoo(5)%id_jingest_sio2 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_fe_vmMdz","Ingestion of Fe by medium-sized migrating zooplankton",&
                           'h','L','s','mol Fe kg-1 s-1','f')
    zoo(4)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_fe_vmLgz","Ingestion of Fe by large migrating zooplankton",&
                           'h','L','s','mol Fe kg-1 s-1','f')
    zoo(5)%id_jingest_fe = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ndet_vmMdz","Production of nitrogen detritus by medium migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ndet_vmLgz","Production of nitrogen detritus by large migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jprod_ndet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_pdet_vmMdz","Production of phosphorous detritus by medium migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_pdet_vmLgz","Production of phosphorous detritus by large migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jprod_pdet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1) 
    vardesc_temp = vardesc("jprod_sidet_vmMdz","Production of opal detritus by medium migrating zooplankton",&
                   'h','L','s','mol SiO2 kg-1 s-1','f')
    zoo(4)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sidet_vmLgz","Production of opal detritus by large migrating zooplankton",&
                   'h','L','s','mol SiO2 kg-1 s-1','f')
    zoo(5)%id_jprod_sidet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sio4_vmMdz","Production of sio4 through grazing/dissolution",&
                   'h','L','s','mol SiO4 kg-1 s-1','f')
    zoo(4)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sio4_vmLgz","Production of sio4 through grazing/dissolution",&
                   'h','L','s','mol SiO4 kg-1 s-1','f')
    zoo(5)%id_jprod_sio4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_fedet_vmMdz","Production of iron detritus by medium migrating zooplankton",&
                   'h','L','s','mol Fe kg-1 s-1','f')
    zoo(4)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_fedet_vmLgz","Production of iron detritus by large migrating zooplankton",&
                   'h','L','s','mol Fe kg-1 s-1','f')
    zoo(5)%id_jprod_fedet = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ldon_vmMdz","Production of labile dissolved organic nitrogen by medium migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ldon_vmLgz","Production of labile dissolved organic nitrogen by large migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jprod_ldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ldop_vmMdz","Production of labile dissolved organic phosphorous by medium migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ldop_vmLgz","Production of labile dissolved organic phosphorous by large migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jprod_ldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_srdon_vmMdz","Production of semi-refractory dissolved organic nitrogen by medium migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_srdon_vmLgz","Production of semi-refractory dissolved organic nitrogen by large migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jprod_srdon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_srdop_vmMdz","Production of semi-refractory dissolved organic phosphorous by medium migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_srdop_vmLgz","Production of semi-refractory dissolved organic phosphorous by large migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jprod_srdop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sldon_vmMdz","Production of semi-labile dissolved organic nitrogen by medium migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sldon_vmLgz","Production of semi-labile dissolved organic nitrogen by large migrating zooplankton",&
                   'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jprod_sldon = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sldop_vmMdz","Production of semi-labile dissolved organic phosphorous by medium migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(4)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_sldop_vmLgz","Production of semi-labile dissolved organic phosphorous by large migrating zooplankton",&
                   'h','L','s','mol P kg-1 s-1','f')
    zoo(5)%id_jprod_sldop = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
     vardesc_temp = vardesc("jprod_fed_vmMdz","Production of dissolved iron by medium migrating zooplankton",&
                   'h','L','s','mol Fe kg-1 s-1','f')
    zoo(4)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_fed_vmLgz","Production of dissolved iron by large migrating zooplankton",&
                   'h','L','s','mol Fe kg-1 s-1','f')
    zoo(5)%id_jprod_fed = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
     vardesc_temp = vardesc("jprod_po4_vmMdz","Production of phosphate by medium migrating zooplankton",&
                   'h','L','s','mol PO4 kg-1 s-1','f')
    zoo(4)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_po4_vmLgz","Production of phosphate by large migrating zooplankton",&
                   'h','L','s','mol PO4 kg-1 s-1','f')
    zoo(5)%id_jprod_po4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_nh4_vmMdz","Production of ammonia by medium migrating zooplankton",&
                   'h','L','s','mol NH4 kg-1 s-1','f')
    zoo(4)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_nh4_vmLgz","Production of ammonia by large migrating zooplankton",&
                   'h','L','s','mol NH4 kg-1 s-1','f')
    zoo(5)%id_jprod_nh4 = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_nvmmdz","Production of new biomass (nitrogen) by medium-sized migrating zooplankton",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(4)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_nvmlgz","Production of new biomass (nitrogen) by large migrating zooplankton",&
                           'h','L','s','mol N kg-1 s-1','f')
    zoo(5)%id_jprod_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("o2lim_vmMdz","Oxygen limitation of medium-sized migrating zooplankton",'h','L','s','dimensionless','f')
    zoo(4)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("o2lim_vmlgz","Oxygen limitation of large migrating zooplankton",'h','L','s','dimensionless','f')
    zoo(5)%id_o2lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("temp_lim_vmMdz","Temperature limitation of medium-sized migrating zooplankton",'h','L','s','dimensionless','f')
    zoo(4)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("temp_lim_vmlgz","Temperature limitation of large migrating zooplankton",'h','L','s','dimensionless','f')
    zoo(5)%id_temp_lim = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("vmove_vmMdz","Vertically migrating medium zooplankton movement",'h','L','s','m s-1','f')
    zoo(4)%id_vmove = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("vmove_vmLgz","Vertically migrating large zooplankton movement",'h','L','s','m s-1','f')
    zoo(5)%id_vmove = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("AE_vmMdz","Assimilation efficiency of migrating medium zooplankton",'h','L','s','dimensionless','f')
    zoo(4)%id_assim_eff = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("AE_vmLgz","Assimilation efficiency of migrating large zooplankton",'h','L','s','dimensionless','f')
    zoo(5)%id_assim_eff = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jmetabo_nvmMdz","Production of nh4 by vertically migrating medium zooplankton metabolism, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(4)%id_jmetabo_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jmetabo_nvmLgz","Production of nh4 by vertically migrating large zooplankton metabolism, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(5)%id_jmetabo_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_nvmmdz_100","Medium migrating zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_nvmlgz_100","Large migrating zooplankton nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(5)%id_jprod_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
   vardesc_temp = vardesc("jingest_n_nvmmdz_100","Medium migrating zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jingest_n_nvmlgz_100","Large migrating zooplankton nitrogen ingestion integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(5)%id_jingest_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jzloss_nvmmdz_100","Medium migrating zooplankton nitrogen loss to zooplankton integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jzloss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jhploss_nvmmdz_100","Medium migrating zooplankton nitrogen loss to higher preds. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jhploss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jhploss_nvmlgz_100","Large migrating zooplankton nitrogen loss to higher preds. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(5)%id_jhploss_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
   vardesc_temp = vardesc("jprod_ndet_nvmmdz_100","Medium migrating zooplankton nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_ndet_nvmlgz_100","Large migrating zooplankton nitrogen detritus prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(5)%id_jprod_ndet_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jprod_don_nvmmdz_100","Medium migrating zooplankton dissolved org. nitrogen prod. integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jprod_don_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jremin_n_nvmmdz_100","Medium migrating zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(4)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jremin_n_nvmlgz_100","Large migrating zooplankton nitrogen remineralization integral in upper 100m",'h','1','s','mol m-2 s-1','f')
    zoo(5)%id_jremin_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("nvmmdz_100","Medium migrating zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(4)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("nvmlgz_100","Large migrating zooplankton nitrogen biomass in upper 100m",'h','1','s','mol m-2','f')
    zoo(5)%id_f_n_100 = register_diag_field(package_name, vardesc_temp%name, axes(1:2),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("vmove_Smz","Small zooplankton movement",'h','L','s','m s-1','f')
    zoo(1)%id_vmove = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("vmove_Mdz","medium zooplankton movement",'h','L','s','m s-1','f')
    zoo(2)%id_vmove = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("vmove_Lgz","Large zooplankton movement",'h','L','s','m s-1','f')
    zoo(3)%id_vmove = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("AE_Smz","Assimilation efficiency of small zooplankton",'h','L','s','dimensionless','f')
    zoo(1)%id_assim_eff = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("AE_Mdz","Assimilation efficiency of medium zooplankton",'h','L','s','dimensionless','f')
    zoo(2)%id_assim_eff = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("AE_Lgz","Assimilation efficiency of large zooplankton",'h','L','s','dimensionless','f')
    zoo(3)%id_assim_eff = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
          init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jmetabo_nSmz","Production of nh4 by small zooplankton metabolism, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(1)%id_jmetabo_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jmetabo_nMdz","Production of nh4 by medium zooplankton metabolism, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(2)%id_jmetabo_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)
    vardesc_temp = vardesc("jmetabo_nLgz","Production of nh4 by large zooplankton metabolism, layer integral",&
                           'h','L','s','mol N m-2 s-1','f')
    zoo(3)%id_jmetabo_n = register_diag_field(package_name, vardesc_temp%name, axes(1:3),&
         init_time, vardesc_temp%longname,vardesc_temp%units, missing_value = missing_value1)


  end subroutine eco_cobalt_reg_diag

  !> Send DVM-specific diagnostics for migrating zooplankton groups (zoo(4) and zoo(5)).
  !! Called from cobalt_send_diagnostics only when do_dvm is true.
  !! When post_vertdiff=.true., computes and sends zoo(4)/zoo(5) 100m integrals and
  !! flux diagnostics. When post_vertdiff=.false., sends source-term rate diagnostics.
  subroutine eco_cobalt_send_diag(tracer_list, model_time, grid_tmask, rho_dzt, &
                                   ilb, jlb, tau, zoo, cobalt, post_vertdiff)
    type(g_tracer_type), pointer                          :: tracer_list
    type(time_type),                         intent(in)   :: model_time
    real, dimension(ilb:,jlb:,:),            intent(in)   :: rho_dzt
    integer,                                 intent(in)   :: ilb, jlb, tau
    type(zooplankton), dimension(NUM_ZOO),   intent(inout) :: zoo
    type(generic_COBALT_type),               intent(inout) :: cobalt
    logical,                                 intent(in), optional :: post_vertdiff
    real, dimension(:,:,:), pointer                       :: grid_tmask

    ! local variables
    integer :: isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau
    integer :: n, i, j, k, k_100
    logical :: used
    logical :: is_post_vertdiff
    real :: drho_dzt
    integer, dimension(:,:), pointer :: mask_coast, grid_kmt
    real, dimension(:,:), allocatable :: rho_dzt_100

    call g_tracer_get_common(isc, iec, jsc, jec, isd, ied, jsd, jed, nk, ntau, &
         grid_tmask=grid_tmask, grid_mask_coast=mask_coast, grid_kmt=grid_kmt)

    is_post_vertdiff = .false.
    if (present(post_vertdiff)) is_post_vertdiff = post_vertdiff

    select case (is_post_vertdiff)

      case (.true.)
        !
        ! Compute zoo(4)/zoo(5) 100m depth-integral fields
        !
        allocate(rho_dzt_100(isc:iec, jsc:jec))
        do j = jsc, jec ; do i = isc, iec  !{
          rho_dzt_100(i,j) = rho_dzt(i,j,1)
          zoo(4)%f_n_100(i,j) = (cobalt%dvm%p_nvmmdz(i,j,1,tau) + cobalt%dvm%p_nvmmdz_gut(i,j,1,tau) + &
            cobalt%dvm%p_nvmmdz_met(i,j,1,tau)) * rho_dzt(i,j,1)
          zoo(5)%f_n_100(i,j) = (cobalt%dvm%p_nvmlgz(i,j,1,tau) + cobalt%dvm%p_nvmlgz_gut(i,j,1,tau) + &
            cobalt%dvm%p_nvmlgz_met(i,j,1,tau)) * rho_dzt(i,j,1)
        enddo; enddo !} i,j

        do j = jsc, jec ; do i = isc, iec  !{
          k_100 = 1
          do k = 2, grid_kmt(i,j)  !{
            if (rho_dzt_100(i,j) .lt. cobalt%Rho_0 * 100.0) then
              k_100 = k
              rho_dzt_100(i,j) = rho_dzt_100(i,j) + rho_dzt(i,j,k)
              zoo(4)%f_n_100(i,j) = zoo(4)%f_n_100(i,j) + &
                (cobalt%dvm%p_nvmmdz(i,j,k,tau) + cobalt%dvm%p_nvmmdz_gut(i,j,k,tau) + &
                 cobalt%dvm%p_nvmmdz_met(i,j,k,tau)) * rho_dzt(i,j,k)
              zoo(5)%f_n_100(i,j) = zoo(5)%f_n_100(i,j) + &
                (cobalt%dvm%p_nvmlgz(i,j,k,tau) + cobalt%dvm%p_nvmlgz_gut(i,j,k,tau) + &
                 cobalt%dvm%p_nvmlgz_met(i,j,k,tau)) * rho_dzt(i,j,k)
            endif
          enddo  !} k

          if (k_100 .gt. 1 .and. k_100 .lt. grid_kmt(i,j)) then
            drho_dzt = cobalt%Rho_0 * 100.0 - rho_dzt_100(i,j)
            zoo(4)%f_n_100(i,j) = zoo(4)%f_n_100(i,j) + &
              (cobalt%dvm%p_nvmmdz(i,j,k_100,tau) + cobalt%dvm%p_nvmmdz_gut(i,j,k_100,tau) + &
               cobalt%dvm%p_nvmmdz_met(i,j,k_100,tau)) * drho_dzt
            zoo(5)%f_n_100(i,j) = zoo(5)%f_n_100(i,j) + &
              (cobalt%dvm%p_nvmlgz(i,j,k_100,tau) + cobalt%dvm%p_nvmlgz_gut(i,j,k_100,tau) + &
               cobalt%dvm%p_nvmlgz_met(i,j,k_100,tau)) * drho_dzt
          endif
        enddo; enddo !} i,j
        deallocate(rho_dzt_100)

        ! Send zoo(4)/zoo(5) 100m biomass integrals
        do n = NUM_BASE_ZOO+1, NUM_ZOO  !{
          used = g_send_data(zoo(n)%id_f_n_100, zoo(n)%f_n_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        enddo !} n

        ! Send zoo(4)/zoo(5) 100m flux integral diagnostics
        do n = NUM_BASE_ZOO+1, NUM_ZOO  !{
          used = g_send_data(zoo(n)%id_jprod_n_100, zoo(n)%jprod_n_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
          used = g_send_data(zoo(n)%id_jingest_n_100, zoo(n)%jingest_n_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
          used = g_send_data(zoo(n)%id_jremin_n_100, zoo(n)%jremin_n_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
          used = g_send_data(zoo(n)%id_jzloss_n_100, zoo(n)%jzloss_n_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
          used = g_send_data(zoo(n)%id_jprod_don_100, zoo(n)%jprod_don_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
          used = g_send_data(zoo(n)%id_jhploss_n_100, zoo(n)%jhploss_n_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
          used = g_send_data(zoo(n)%id_jprod_ndet_100, zoo(n)%jprod_ndet_100, &
            model_time, rmask = grid_tmask(:,:,1), is_in=isc, js_in=jsc, ie_in=iec, je_in=jec)
        enddo !} n

      case (.false.)
        !
        ! Send zoo(4)/zoo(5) source-term rate diagnostics
        !
        do n = NUM_BASE_ZOO+1, NUM_ZOO  !{
          used = g_send_data(zoo(n)%id_jzloss_n, zoo(n)%jzloss_n, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jhploss_n, zoo(n)%jhploss_n, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jzloss_p, zoo(n)%jzloss_p, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jhploss_p, zoo(n)%jhploss_p, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jingest_n, zoo(n)%jingest_n, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jingest_p, zoo(n)%jingest_p, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jingest_sio2, zoo(n)%jingest_sio2, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jingest_fe, zoo(n)%jingest_fe, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_ndet, zoo(n)%jprod_ndet, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_pdet, zoo(n)%jprod_pdet, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_ldon, zoo(n)%jprod_ldon, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_ldop, zoo(n)%jprod_ldop, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_sldon, zoo(n)%jprod_sldon, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_sldop, zoo(n)%jprod_sldop, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_srdon, zoo(n)%jprod_srdon, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_srdop, zoo(n)%jprod_srdop, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_fed,  zoo(n)%jprod_fed, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_fedet, zoo(n)%jprod_fedet, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_sidet, zoo(n)%jprod_sidet, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_sio4, zoo(n)%jprod_sio4, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_po4,  zoo(n)%jprod_po4, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_nh4,  zoo(n)%jprod_nh4, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_jprod_n, zoo(n)%jprod_n, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_o2lim, zoo(n)%o2lim, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
          used = g_send_data(zoo(n)%id_temp_lim, zoo(n)%temp_lim, &
            model_time, rmask = grid_tmask, is_in=isc, js_in=jsc, ks_in=1,ie_in=iec, je_in=jec, ke_in=nk)
        enddo !} n

    end select

  end subroutine eco_cobalt_send_diag

end module COBALT_eco
