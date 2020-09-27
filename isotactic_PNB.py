from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk_tacticity
from pysimm.models.monomers.dreiding.NbTMS_H2_tacticity import monomer as NbTMS
from pysimm.apps.equilibrate import equil

def run(test=False):
    # we'll make a polystyrene monomer from the pysimm models database
    A = NbTMS(isomer="exoexo")
    # we'll instantiate a Dreiding forcefield object for use later
    f = forcefield.Dreiding()
    
    # the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
    A.apply_charges(f, charges='gasteiger')
    
    # the buckingham potential isn't great at small distances, and therefore we use the LJ potential while growing the polymer
    A.pair_style = 'lj/cut'
    
    # run the polymer random walk tacticity method with 10 total repeat units
    polymer = random_walk_tacticity(A, 20, forcefield=f,capped=True,tacticity='syndiotactic',rotation=180,errorCheck=False,sim=0)
    writeFileFormats(polymer,"polymer",unwrap=True)

    #quick opt of polymer
    lmps.quick_min(polymer, min_style='fire',etol=1.0e-4,maxiter=100000)

    # write a few different file formats
    polymer.unwrap()
    writeFileFormats(polymer,"polymer_fired")
    
    #pack multiple copies of polymer
    polymers = system.replicate(polymer, 8, density=0.005)
    #polymers = polymer
    writeFileFormats(polymers,"polymers")
    lmps.quick_min(polymers, min_style='fire',etol=1.0e-4,maxiter=100000)
    writeFileFormats(polymers,"polymers_fired")
    
    #quickmd
    nvt_settings = {
        'name': 'nvt_md',
        'print_to_screen': True,
        'ensemble': 'nvt',
        'temperature': {
            'start': 100,
            'stop': 300
        },
        'new_v': True,
        'length': 10000
    }
    npt_settings = {
        'name': 'npt_md',
        'print_to_screen': True,
        'ensemble': 'npt',
        'temperature': 300,
        'new_v': True,
        'pressure': {
            'start': 1000,
            'stop': 1
        },
        'length': 100000,
        'thermo_style': 'custom step temp press density'
    }
    #nvt runs okay, but npt fails...need to add neigh_modify command to reneighbor more often during compression of npt step
    #lmps.quick_md(polymers, debug=True, **nvt_settings)
    #writeFileFormats(polymers,"polymers_nvt")
    #lmps.quick_md(polymers, debug=True, **npt_settings)
    #lmps.quick_md(polymers)
    #writeFileFormats(polymers,"polymers_npt")
    sim = lmps.Simulation(polymers, name='nptAndReneighbor', debug=True)
    sim.add_custom('neigh_modify delay 0')
    sim.add(lmps.Velocity(temperature=1000))
    sim.add_md(length=10000,ensemble='npt',temperature=1000,pressure=5000)
    sim.run()
    writeFileFormats(polymers,"polymers_npt")
    writeFileFormats(polymers,"polymers_npt_unwrapped",unwrap=True)


    #21-step equilibration
    equil(polymers,np=1,pmax=50000)
    writeFileFormats(polymers,"polymers_equil")
    polymers.unwrap()
    writeFileFormats(polymers,"polymers_equil_unwrap")
        
def writeFileFormats(s,name,**kwargs):
    unwrap = kwargs.get('unwrap',False)
    if unwrap:
        s.unwrap()
    s.write_xyz(name + '.xyz')
    s.write_yaml(name + '.yaml')
    s.write_lammps(name + '.lmps')
    s.write_chemdoodle_json(name + '.json')
    if unwrap:
        s.wrap()

if __name__ == '__main__':
    run()
