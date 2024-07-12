from Code.Codebase import Case_Maker

AVL_file  = 'TitanWing.avl'
Case_file = 'Nu_Run_File'
Case_name = 'Cruise'

Case_vcv = [
    'a a 0'
,   'b b 0'
,   'r r 0'
,   'p p 0'
,   'y y 0'
,   'd1 d1 0'
]

Case_flight_condition = {
    'alpha'     : 0.0
,   'beta'      : 0.0
,   'pb/2V'     : 0.0
,   'qc/2V'     : 0.0
,   'rb/2V'     : 0.0
,   'CL'        : 0.0
,   'CDo'       : 0.0
,   'bank'      : 0.0
,   'elevation' : 0.0
,   'heading'   : 0.0
,   'Mach'      : 0.0
,   'velocity'  : 444.5331
,   'density'   : 23.77e-4
,   'grav.acc.' : 32.17
,   'turn_rad.' : 0.0
,   'load_fac.' : 0.0
}

Case_mass_properties = {
    'X_cg'      : -1
,   'Y_cg'      : -1
,   'Z_cg'      : -1
,   'mass'      : 32451
,   'Ixx'       : 1.0
,   'Iyy'       : 1.0
,   'Izz'       : 1.0
,   'Ixy'       : 0.0
,   'Iyz'       : 0.0
,   'Izx'       : 0.0
}

Case_viscous_properties = {
    'visc CL_a' : 0.0
,   'visc CL_u' : 0.0
,   'visc CM_a' : 0.0
,   'visc CM_u' : 0.0
}

#Don't need to edit this
cm = Case_Maker(avl_file=AVL_file, run_file=Case_file)
cm._open_existing_dot_run()
return_code = cm.insert_new_run_case(case_name=Case_name,
                                     vcv=Case_vcv,
                                     flight=Case_flight_condition,
                                     mass=Case_mass_properties,
                                     viscous=Case_viscous_properties)
if return_code == 100: cm._write_the_dot_run_file()
else: pass
