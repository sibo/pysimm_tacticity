from pysimm.apps import mc_md
from pysimm import system

frame = system.read_lammps('irmof1_drei.lmps')
gas1 = system.read_lammps('ch4.lmps')

mc_props = {'rigid_type': False,
            'max_ins': 1000,
            'Chemical_Potential_Info': -30.09,
            'Temperature_Info': 300,
            'Run_Type': {'steps': 250},
            'CBMC_Info': {'rcut_cbmc': 2.0},
            'Simulation_Length_Info': {'run': 10000,
                                       'coord_freq': 10000,
                                       'prop_freq': 500},
            'VDW_Style': {'cut_val': 9.0},
            'Charge_Style': {'cut_val': 9.0},
            'Property_Info': {'prop1': 'energy_total',
                              'prop2': 'pressure'}}

md_props = {'ensemble': 'npt',
            'temp': 300,
            'pressure': 15,
            'timestep': 0.25,
            'cutoff': 9.0,
            'length': 10000,
            'thermo': 2500,
            'dump': 2500,
            'print_to_screen': False}

sim_result = mc_md.mc_md(gas1, frame, mcmd_niter=3, sim_folder='results', mc_props=mc_props, md_props=md_props)
sim_result.write_lammps('MOFplusME.lmps')
sim_result.write_xyz('MOFplusME.xyz')
