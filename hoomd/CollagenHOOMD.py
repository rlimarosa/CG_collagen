# Collagen CG modeling in HOOMD blue
import argparse
import datetime
import itertools
import os
import sys
import yaml

import numpy as np
import mdtraj

import hoomd
import hoomd.md as md
import gsd.hoomd

# Dictionary lookups for the types in the simulation
getTypebyName = {}
getTypebyNameBonds = {}
getTypebyNameAngles = {}
getTypebyNameDihedrals = {}
getNamebyType = {}
sigma_per_nm = 2.0

r_buffer_cell = 0.5
timestep_size = 0.001
t_damp = 0.05

sigma_lj = {}
sigma_lj[('BB', 'BB')] = 1.0
sigma_lj[('BB', 'HB')] = 0.65
sigma_lj[('HBa', 'HBd')] = 0.3
sigma_lj[('HBa', 'HBa')] = 0.7
sigma_lj[('HBd', 'HBd')] = 0.7

cutoff_lj = {}
cutoff_lj[('BB', 'BB')] = np.power(2.0, 1/6.0) * sigma_lj[('BB', 'BB')]
cutoff_lj[('BB', 'HB')] = np.power(2.0, 1/6.0) * sigma_lj[('BB', 'HB')]
cutoff_lj[('HBa', 'HBd')] = 0.570
cutoff_lj[('HBa', 'HBa')] = np.power(2.0, 1/6.0) * sigma_lj[('HBa', 'HBa')]
cutoff_lj[('HBd', 'HBd')] = np.power(2.0, 1/6.0) * sigma_lj[('HBd', 'HBd')]

epsilon_lj = {}
epsilon_lj[('BB', 'BB')] = 1.0
epsilon_lj[('BB', 'HB')] = 1.0
epsilon_lj[('HBa', 'HBd')] = 50.4
epsilon_lj[('HBa', 'HBa')] = 1.0
epsilon_lj[('HBd', 'HBd')] = 1.0

# Create a status line maker for our output
class Status():

    def __init__(self, sim):
        self.sim = sim

    @property
    def seconds_remaining(self):
        try:
            return (self.sim.final_timestep - self.sim.timestep) / self.sim.tps
        except ZeroDivisionError:
            return 0

    @property
    def etr(self):
        return str(datetime.timedelta(seconds=self.seconds_remaining))

def parse_args():
    parser = argparse.ArgumentParser(prog='CollagenHOOMD.py')

    # Options
    parser.add_argument('--pdb', required=True, help="Path to PDB file")
    parser.add_argument('--mode',   choices=['cpu','gpu'], required=True, help="Compute mode: cpu or gpu")
    parser.add_argument('--box',    type=float, default=110.0, help="Box edge length in sigma (default: 110)")
    opts = parser.parse_args()

    return opts

def print_snapshot(snap):
    if snap.communicator.rank != 0:
        return

    # Alias the arrays for ease
    pos = snap.particles.position  # (N,3)
    tid = snap.particles.typeid
    types = snap.particles.types

    print(f"\nSnapshot: {len(pos)} particles")
    print(f"{'idx':>5s}  {'type':>4s}   {'x':>8s}   {'y':>8s}   {'z':>8s}")
    print('-'*45)
    for i, (x,y,z) in enumerate(pos):
        tname = types[tid[i]]
        print(f"{i:5d}  {tname:>4s}   {x:8.3f}   {y:8.3f}   {z:8.3f}")

def load_pdb_to_snapshot(pdb_path):
    print("Loading PDB into HOOMD snapshot")
    # Load the PDB file with mdtraj
    traj = mdtraj.load_pdb(pdb_path)
    positions = traj.xyz[0]  # shape = (n_atoms,3)
    # Build the names of the atoms, and corresponding typeid (numeric)
    current_typeid = 0
    for atom in traj.topology.atoms:
        if atom.name not in getTypebyName:
            getTypebyName[atom.name] = current_typeid
            getNamebyType[current_typeid] = atom.name
            current_typeid += 1

    # Build an empty hoomd snapshot
    snap = hoomd.Snapshot()
    # This is only accessible on rank 0 for the positions...
    if snap.communicator.rank == 0:
        snap.particles.N = traj.n_atoms
        snap.particles.types = list(getTypebyName.keys())
        snap.particles.position[:] = positions
        # Loop over particles to assign identity
        current_atom_idx = 0
        for atom in traj.topology.atoms:
            snap.particles.typeid[current_atom_idx] = getTypebyName[atom.name]
            if atom.name == 'BBP' or atom.name == 'BBO' or atom.name == 'BBG':
                snap.particles.mass[current_atom_idx] = 3.0
            elif atom.name == 'HBP' or atom.name == 'HBG':
                snap.particles.mass[current_atom_idx] = 1.0
            else:
                print("Encountered an atom name we don't know, exiting")
                sys.exit(1)
            current_atom_idx += 1

    print("Finished reading PDB for coordinates and identity of particles")
    return snap

def create_triple_helix_from_pdb(pdb_path, box_length):
    print("Creating triple helix from PDB")
    # First create the snapshot associated with the positions and types in the PDB file
    snap = load_pdb_to_snapshot(pdb_path)

    # Bonds are set up in the following way for each 3(5)-mer
    # HBP           HBG
    #  |             |
    # BBP -- BBO -- BBG ->
    nchains = 3
    nmer = 12
    nbeads_per_mer = 5
    nbonds_per_mer = 4
    nangles_per_mer = 3
    ndihedrals_per_mer = 1

    # Linear bonds
    # Intra-mer
    getTypebyNameBonds['BBP_BBO'] = 0
    getTypebyNameBonds['BBO_BBG'] = 1
    getTypebyNameBonds['BBP_HBP'] = 2
    getTypebyNameBonds['BBG_HBG'] = 3
    # Inter-mer
    getTypebyNameBonds['BBG_BBP'] = 4
    # Now figure out how many bonds total we have
    snap.bonds.N = nchains*((nmer * nbonds_per_mer) + (nmer - 1))
    snap.bonds.types = ['BBP_BBO', 'BBO_BBG', 'BBP_HBP', 'BBG_HBG', 'BBG_BBP']

    # Angular bonds
    # Intra-mer
    getTypebyNameAngles['BBP_BBO_BBG'] = 0
    getTypebyNameAngles['HBP_BBP_BBO'] = 1
    getTypebyNameAngles['HBG_BBG_BBO'] = 2
    # Inter-mer
    getTypebyNameAngles['BBO_BBG_BBP'] = 3
    getTypebyNameAngles['BBG_BBP_BBO'] = 4
    # Now figure out how many angles total we have
    snap.angles.N = nchains*((nmer * nangles_per_mer) + 2*(nmer - 1))
    snap.angles.types = ['BBP_BBO_BBG', 'HBP_BBP_BBO', 'HBG_BBG_BBO',
                         'BBO_BBG_BBP', 'BBG_BBP_BBO']

    # Dihedrals
    # https://hoomd-blue.readthedocs.io/en/latest/hoomd/md/dihedral/periodic.html
    # Inter-mer
    getTypebyNameDihedrals['HBP_BBP_BBG_HBG'] = 0
    # Intra-mer
    getTypebyNameDihedrals['HBG_BBG_BBP_HBP'] = 1
    # Now figure out how many dihedrals we have total
    snap.dihedrals.N = nchains*((nmer * ndihedrals_per_mer) + (nmer - 1))
    snap.dihedrals.types = ['HBP_BBP_BBG_HBG', 'HBG_BBG_BBP_HBP']

    # March down each chain individually, as this will change the total numbers
    ibond = 0
    iangle = 0
    idihedral = 0
    for ichain in range(3):
        print("Assigning bonds to chain {}".format(ichain))
        # March down the -mer chain adding the bonds as we go
        for imer in range(nmer):
            print("  Creating i-mer {}".format(imer))
            # Do a sanity check to make sure that the types from particles matches what we expect in the bonds
            BBPidx = (ichain * nmer * nbeads_per_mer) + (imer * nbeads_per_mer)
            HBPidx = BBPidx + 1
            BBOidx = BBPidx + 2
            BBGidx = BBPidx + 3
            HBGidx = BBPidx + 4
            if snap.particles.types[snap.particles.typeid[BBPidx]] != 'BBP':
                print("ERROR in BBP identity!")
                sys.exit(1)
            if snap.particles.types[snap.particles.typeid[HBPidx]] != 'HBP':
                print("ERROR in HBP identity!")
                sys.exit(1)
            if snap.particles.types[snap.particles.typeid[BBOidx]] != 'BBO':
                print("ERROR in BBO identity!")
                sys.exit(1)
            if snap.particles.types[snap.particles.typeid[BBGidx]] != 'BBG':
                print("ERROR in BBG identity!")
                sys.exit(1)
            if snap.particles.types[snap.particles.typeid[HBGidx]] != 'HBG':
                print("ERROR in HBG identity!")
                sys.exit(1)
            # Now we can assign the bonds
            snap.bonds.typeid[ibond] = getTypebyNameBonds['BBP_BBO']
            snap.bonds.group[ibond] = [BBPidx, BBOidx]
            ibond += 1
            snap.bonds.typeid[ibond] = getTypebyNameBonds['BBO_BBG']
            snap.bonds.group[ibond] = [BBOidx, BBGidx]
            ibond += 1
            snap.bonds.typeid[ibond] = getTypebyNameBonds['BBP_HBP']
            snap.bonds.group[ibond] = [BBPidx, HBPidx]
            ibond += 1
            snap.bonds.typeid[ibond] = getTypebyNameBonds['BBG_HBG']
            snap.bonds.group[ibond] = [BBGidx, HBGidx]
            ibond += 1
            # Angles
            snap.angles.typeid[iangle] = getTypebyNameAngles['BBP_BBO_BBG']
            snap.angles.group[iangle] = [BBPidx, BBOidx, BBGidx]
            iangle += 1
            snap.angles.typeid[iangle] = getTypebyNameAngles['HBP_BBP_BBO']
            snap.angles.group[iangle] = [HBPidx, BBPidx, BBOidx]
            iangle += 1
            snap.angles.typeid[iangle] = getTypebyNameAngles['HBG_BBG_BBO']
            snap.angles.group[iangle] = [HBGidx, BBGidx, BBOidx]
            iangle += 1
            # Dihedrals
            snap.dihedrals.typeid[idihedral] = getTypebyNameDihedrals['HBP_BBP_BBG_HBG']
            snap.dihedrals.group[idihedral] = [HBPidx, BBPidx, BBGidx, HBGidx]
            idihedral += 1
            # If we are not the last imer, add the bond to the next section too
            if imer < nmer - 1:
                # Bonds
                snap.bonds.typeid[ibond] = getTypebyNameBonds['BBG_BBP']
                snap.bonds.group[ibond] = [BBGidx, BBPidx + nbeads_per_mer]
                ibond += 1
                # Angles
                snap.angles.typeid[iangle] = getTypebyNameAngles['BBO_BBG_BBP']
                snap.angles.group[iangle] = [BBOidx, BBGidx, BBPidx + nbeads_per_mer]
                iangle += 1
                snap.angles.typeid[iangle] = getTypebyNameAngles['BBG_BBP_BBO']
                snap.angles.group[iangle] = [BBGidx, BBPidx + nbeads_per_mer, BBOidx + nbeads_per_mer]
                iangle += 1
                # Dihedrals
                snap.dihedrals.typeid[idihedral] = getTypebyNameDihedrals['HBG_BBG_BBP_HBP']
                snap.dihedrals.group[idihedral] = [HBGidx, BBGidx, BBPidx + nbeads_per_mer, HBPidx + nbeads_per_mer]
                idihedral += 1

    snap.configuration.box = hoomd.Box.cube(L=box_length)

    return snap

if __name__ == "__main__":
    opts = parse_args()
    device = (hoomd.device.GPU(notice_level=3) if opts.mode=='gpu'
              else hoomd.device.CPU(notice_level=3))
    sim = hoomd.Simulation(device=device, seed=1)
    print("MPI enabled:", hoomd.version.mpi_enabled)

    # build & inject snapshot
    #snap = load_pdb_to_snapshot(opts.pdb, opts.box)
    snap = create_triple_helix_from_pdb(opts.pdb, opts.box)
    print_snapshot(snap)
    sim.create_state_from_snapshot(snap)

    # Actually try to run the simulation, first, create the cell list
    cell = hoomd.md.nlist.Cell(buffer=r_buffer_cell, exclusions = ['bond'])

    ###############################
    # Bonded, angle, and dihedral interactions
    ###############################
    # Assign bonded interaction strengths
    # Bonds
    linear_bond = md.bond.Harmonic()
    # Intra
    linear_bond.params['BBP_BBO'] = dict(k=1000.0, r0=0.5)
    linear_bond.params['BBO_BBG'] = dict(k=1000.0, r0=0.5)
    linear_bond.params['BBP_HBP'] = dict(k=1000.0, r0=0.375)
    linear_bond.params['BBG_HBG'] = dict(k=1000.0, r0=0.375)
    # Inter
    linear_bond.params['BBG_BBP'] = dict(k=1000.0, r0=0.5)

    # Angles
    angular_bond = md.angle.Harmonic()
    # Intra
    angular_bond.params['BBP_BBO_BBG'] = dict(k=30.0, t0=np.pi)
    angular_bond.params['HBP_BBP_BBO'] = dict(k=300.0, t0=np.pi/2.0)
    angular_bond.params['HBG_BBG_BBO'] = dict(k=300.0, t0=np.pi/2.0)
    # Inter
    angular_bond.params['BBO_BBG_BBP'] = dict(k=30.0, t0=np.pi)
    angular_bond.params['BBG_BBP_BBO'] = dict(k=30.0, t0=np.pi)

    # Dihedrals
    dihedral_bond = md.dihedral.Periodic()
    # Intra
    dihedral_bond.params['HBP_BBP_BBG_HBG'] = dict(k=15.0, d=1, n=1, phi0=2.0*np.pi/3.0)
    # Inter
    dihedral_bond.params['HBG_BBG_BBP_HBP'] = dict(k=15.0, d=-1, n=1, phi0=2.0*np.pi/3.0)


    ###############################
    # Non-bonded interactions
    ###############################
    wca = md.pair.LJ(nlist = cell)
    wca.mode = 'shift'
    # BBO <--> BBO
    wca.params[('BBO', 'BBO')] = dict(epsilon=epsilon_lj[('BB', 'BB')], sigma=sigma_lj[('BB', 'BB')])
    wca.r_cut['BBO', 'BBO'] = cutoff_lj[('BB', 'BB')]
    # BBO <--> BBP
    wca.params[('BBO', 'BBP')] = dict(epsilon=epsilon_lj[('BB', 'BB')], sigma=sigma_lj[('BB', 'BB')])
    wca.r_cut['BBO', 'BBP'] = cutoff_lj[('BB', 'BB')]
    # BBO <--> BBG
    wca.params[('BBO', 'BBG')] = dict(epsilon=epsilon_lj[('BB', 'BB')], sigma=sigma_lj[('BB', 'BB')])
    wca.r_cut['BBO', 'BBG'] = cutoff_lj[('BB', 'BB')]
    # BBO <--> HBP
    wca.params[('BBO', 'HBP')] = dict(epsilon=epsilon_lj[('BB', 'HB')], sigma=sigma_lj[('BB', 'HB')])
    wca.r_cut['BBO', 'HBP'] = cutoff_lj[('BB', 'HB')]
    # BBO <--> HBG
    wca.params[('BBO', 'HBG')] = dict(epsilon=epsilon_lj[('BB', 'HB')], sigma=sigma_lj[('BB', 'HB')])
    wca.r_cut['BBO', 'HBG'] = cutoff_lj[('BB', 'HB')]

    # BBP <--> BBP
    wca.params[('BBP', 'BBP')] = dict(epsilon=epsilon_lj[('BB', 'BB')], sigma=sigma_lj[('BB', 'BB')])
    wca.r_cut['BBP', 'BBP'] = cutoff_lj[('BB', 'BB')]
    # BBP <--> BBG
    wca.params[('BBP', 'BBG')] = dict(epsilon=epsilon_lj[('BB', 'BB')], sigma=sigma_lj[('BB', 'BB')])
    wca.r_cut['BBP', 'BBG'] = cutoff_lj[('BB', 'BB')]
    # BBP <--> HBP
    wca.params[('BBP', 'HBP')] = dict(epsilon=epsilon_lj[('BB', 'HB')], sigma=sigma_lj[('BB', 'HB')])
    wca.r_cut['BBP', 'HBP'] = cutoff_lj[('BB', 'HB')]
    # BBP <--> HBG
    wca.params[('BBP', 'HBG')] = dict(epsilon=epsilon_lj[('BB', 'HB')], sigma=sigma_lj[('BB', 'HB')])
    wca.r_cut['BBP', 'HBG'] = cutoff_lj[('BB', 'HB')]

    # BBG <--> BBG
    wca.params[('BBG', 'BBG')] = dict(epsilon=epsilon_lj[('BB', 'BB')], sigma=sigma_lj[('BB', 'BB')])
    wca.r_cut['BBG', 'BBG'] = cutoff_lj[('BB', 'BB')]
    # BBG <--> HBP
    wca.params[('BBG', 'HBP')] = dict(epsilon=epsilon_lj[('BB', 'HB')], sigma=sigma_lj[('BB', 'HB')])
    wca.r_cut['BBG', 'HBP'] = cutoff_lj[('BB', 'HB')]
    # BBG <--> HBG
    wca.params[('BBG', 'HBG')] = dict(epsilon=epsilon_lj[('BB', 'HB')], sigma=sigma_lj[('BB', 'HB')])
    wca.r_cut['BBG', 'HBG'] = cutoff_lj[('BB', 'HB')]

    # HBP <--> HBP
    wca.params[('HBP', 'HBP')] = dict(epsilon=epsilon_lj[('HBa', 'HBa')], sigma=sigma_lj[('HBa', 'HBa')])
    wca.r_cut['HBP', 'HBP'] = cutoff_lj[('HBa', 'HBa')]
    # HBP <--> HBG
    wca.params[('HBP', 'HBG')] = dict(epsilon=epsilon_lj[('HBa', 'HBd')], sigma=sigma_lj[('HBa', 'HBd')])
    wca.r_cut['HBP', 'HBG'] = cutoff_lj[('HBa', 'HBd')]

    # HBG <--> HBG
    wca.params[('HBG', 'HBG')] = dict(epsilon=epsilon_lj[('HBd', 'HBd')], sigma=sigma_lj[('HBd', 'HBd')])
    wca.r_cut['HBG', 'HBG'] = cutoff_lj[('HBd', 'HBd')]

    ###############################
    # Integrator, Langevin thermostat
    ###############################
    integrator = md.Integrator(dt=timestep_size)
    integrator.forces.append(linear_bond)
    integrator.forces.append(angular_bond)
    integrator.forces.append(dihedral_bond)
    integrator.forces.append(wca)

    langevin = md.methods.Langevin(filter=hoomd.filter.All(), kT = 1.0)
    # XXX This needs double checking for the conversion to damping coefficient
    langevin.gamma['BBP'] = 3.0 / t_damp
    langevin.gamma['BBO'] = 3.0 / t_damp
    langevin.gamma['BBG'] = 3.0 / t_damp
    langevin.gamma['HBP'] = 1.0 / t_damp
    langevin.gamma['HBG'] = 1.0 / t_damp

    integrator.methods.append(langevin)

    # Add everything to the simulation system
    sim.operations.integrator = integrator

    # Thermalize the system
    sim.run(0)
    sim.state.thermalize_particle_momenta(hoomd.filter.All(), kT=1.0)

    ###############################################################################
    # Print information for the main program
    ###############################################################################
    # Keep track of the thermodynamic information
    thermodynamic_properties = md.compute.ThermodynamicQuantities(
            filter = hoomd.filter.All())
    sim.operations.computes.append(thermodynamic_properties)

    logger = hoomd.logging.Logger()
    logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
    logger.add(thermodynamic_properties)
    
    # Display some quantities to a table while running
    output_logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    status = Status(sim)
    output_logger.add(sim, quantities=['timestep', 'tps'])
    output_logger[('Status', 'etr')] = (status, 'etr', 'string')
    output_logger.add(thermodynamic_properties, quantities=['kinetic_temperature', 'pressure'])
    table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(period=1000),
                              logger=output_logger)
    sim.operations.writers.append(table)

    # Set up writing out GSD trajectories
    gsd_writer = hoomd.write.GSD(filename = 'collagen.gsd',
                                    trigger = hoomd.trigger.Periodic(1000),
                                    mode = 'wb',
                                    filter = hoomd.filter.All(),
                                    logger = logger)
    sim.operations.writers.append(gsd_writer)

    ###############################################################################
    # Run the simulation
    ###############################################################################
    print(f"--------")
    sim.run(1e6)
