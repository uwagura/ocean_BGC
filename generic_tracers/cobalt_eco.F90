!> The cobalt_eco module contains subroutines for the diel vertical migration (DVM)
!! component of the COBALT biogeochemical model. The DVM-specific computation is
!! activated via the do_dvm flag in generic_COBALT_nml, but DVM tracer fields are
!! always registered to maintain a consistent NUM_ZOO=5 data layout.
module COBALT_eco

  use cobalt_types
  use g_tracer_utils, only : g_tracer_add, g_tracer_type
  use g_tracer_utils, only : g_send_data, g_tracer_get_common, g_tracer_get_pointer
  use g_tracer_utils, only : g_diag_type
  use g_tracer_utils, only : register_diag_field=>g_register_diag_field
  use time_manager_mod, only: time_type

  implicit none; private

  public cobalt_eco_add_tracers
  public eco_cobalt_reg_diag
  public eco_cobalt_send_diag

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
          zoo(4)%f_n_100(i,j) = (cobalt%p_nvmmdz(i,j,1,tau) + cobalt%p_nvmmdz_gut(i,j,1,tau) + &
            cobalt%p_nvmmdz_met(i,j,1,tau)) * rho_dzt(i,j,1)
          zoo(5)%f_n_100(i,j) = (cobalt%p_nvmlgz(i,j,1,tau) + cobalt%p_nvmlgz_gut(i,j,1,tau) + &
            cobalt%p_nvmlgz_met(i,j,1,tau)) * rho_dzt(i,j,1)
        enddo; enddo !} i,j

        do j = jsc, jec ; do i = isc, iec  !{
          k_100 = 1
          do k = 2, grid_kmt(i,j)  !{
            if (rho_dzt_100(i,j) .lt. cobalt%Rho_0 * 100.0) then
              k_100 = k
              rho_dzt_100(i,j) = rho_dzt_100(i,j) + rho_dzt(i,j,k)
              zoo(4)%f_n_100(i,j) = zoo(4)%f_n_100(i,j) + &
                (cobalt%p_nvmmdz(i,j,k,tau) + cobalt%p_nvmmdz_gut(i,j,k,tau) + &
                 cobalt%p_nvmmdz_met(i,j,k,tau)) * rho_dzt(i,j,k)
              zoo(5)%f_n_100(i,j) = zoo(5)%f_n_100(i,j) + &
                (cobalt%p_nvmlgz(i,j,k,tau) + cobalt%p_nvmlgz_gut(i,j,k,tau) + &
                 cobalt%p_nvmlgz_met(i,j,k,tau)) * rho_dzt(i,j,k)
            endif
          enddo  !} k

          if (k_100 .gt. 1 .and. k_100 .lt. grid_kmt(i,j)) then
            drho_dzt = cobalt%Rho_0 * 100.0 - rho_dzt_100(i,j)
            zoo(4)%f_n_100(i,j) = zoo(4)%f_n_100(i,j) + &
              (cobalt%p_nvmmdz(i,j,k_100,tau) + cobalt%p_nvmmdz_gut(i,j,k_100,tau) + &
               cobalt%p_nvmmdz_met(i,j,k_100,tau)) * drho_dzt
            zoo(5)%f_n_100(i,j) = zoo(5)%f_n_100(i,j) + &
              (cobalt%p_nvmlgz(i,j,k_100,tau) + cobalt%p_nvmlgz_gut(i,j,k_100,tau) + &
               cobalt%p_nvmlgz_met(i,j,k_100,tau)) * drho_dzt
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
