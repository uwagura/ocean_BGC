!> The cobalt_eco module contains subroutines for the diel vertical migration (DVM)
!! component of the COBALT biogeochemical model. The DVM-specific computation is
!! activated via the do_dvm flag in generic_COBALT_nml, but DVM tracer fields are
!! always registered to maintain a consistent NUM_ZOO=5 data layout.
module COBALT_eco

  use cobalt_types
  use g_tracer_utils, only : g_tracer_add, g_tracer_type

  implicit none; private

  public cobalt_eco_add_tracers

  contains

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

end module COBALT_eco
