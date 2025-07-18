import argparse
import itertools
import string
from sys import stdout
from simtk import openmm, unit
from openmm import app
from openmm import app as appmod
from openmm import XmlSerializer
from openmm.app import Topology
import numpy as np

# =============================================================================
# Parse command-line arguments
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--pdb", default="POG20_1trimervs3.pdb")
parser.add_argument("--outprefix", default="POG20_300K")
parser.add_argument("--temp", type=float, default=300.0,help="Simulation temperature in K")
args = parser.parse_args()

# =============================================================================
# Load structure and build box
# =============================================================================
pdb = app.PDBFile(args.pdb)
modeller = app.Modeller(pdb.topology, pdb.positions)

# Set box size
modeller.topology.setPeriodicBoxVectors(np.eye(3) * 50.0 * unit.nanometers)


# =============================================================================
# Build the system with custom parameters
# =============================================================================
topology = modeller.topology
positions = modeller.positions

system = openmm.System()
for atom in topology.atoms():
    name = atom.name
    mass = {"BBP": 3.0, "BBO": 3.0, "BBG": 3.0, "HBP": 1.0, "HBG": 1.0}.get(name, 1.0)
    system.addParticle(mass * unit.amu)


system.setDefaultPeriodicBoxVectors(*topology.getPeriodicBoxVectors())
bond_force = openmm.HarmonicBondForce()
angle_force = openmm.HarmonicAngleForce()
torsion_force = openmm.PeriodicTorsionForce()

bond_params = {
    ("BBP", "BBO"): (0.25, 100),
    ("BBO", "BBG"): (0.25, 100),
    ("BBG", "BBP"): (0.25, 100),
    ("BBG", "HBG"): (0.185, 100),
    ("BBP", "HBP"): (0.185, 100)
    
}
angle_params = [
    ("BBP", "BBO", "BBG", 3.14159, 2),
    ("BBG", "BBP", "BBO", 3.14159, 2),
    ("BBO", "BBG", "BBP", 3.14159, 2),
    ("HBP", "BBP", "BBO", 1.5708, 30),
    ("HBG", "BBG", "BBO", 1.5708, 30)
]

torsion_params = [
    ("HBP", "BBP", "BBG", "HBG", 1, 2.0944, 30),
    ("HBG", "BBG", "BBP", "HBP", 1, -2.0944, 30)
]


# Build bonds, angles, and torsions inside each chain
for chain in topology.chains():
    residues = list(chain.residues())

    for i, res in enumerate(residues):
        atoms = {a.name: a.index for a in res.atoms()}

        # Intra-residue bonds
        for (a1, a2), (length, k) in bond_params.items():
            if a1 in atoms and a2 in atoms:
                bond_force.addBond(atoms[a1], atoms[a2], length * unit.nanometer, k * unit.kilocalorie_per_mole / unit.nanometer**2)

        # Intra-residue angles
        for (a1, a2, a3, theta0, k) in angle_params:
            if a1 in atoms and a2 in atoms and a3 in atoms:
                angle_force.addAngle(atoms[a1], atoms[a2], atoms[a3],
                                     theta0 * unit.radian,
                                     k * unit.kilocalorie_per_mole / unit.radian**2)

        # Intra-residue torsions
        for (a1, a2, a3, a4, periodicity, phase, k) in torsion_params:
            if a1 in atoms and a2 in atoms and a3 in atoms and a4 in atoms:
                torsion_force.addTorsion(atoms[a1], atoms[a2], atoms[a3], atoms[a4],
                                         periodicity,
                                         phase * unit.radian,
                                         k * unit.kilocalorie_per_mole)

        # Inter-residue bond (connect BBG[i] → BBP[i+1])
        if i < len(residues) - 1:
            this_atoms = {a.name: a.index for a in residues[i].atoms()}
            next_atoms = {a.name: a.index for a in residues[i+1].atoms()}
            if "BBG" in this_atoms and "BBP" in next_atoms:
                bond_force.addBond(this_atoms["BBG"], next_atoms["BBP"],
                                   0.25 * unit.nanometer, 100 * unit.kilocalorie_per_mole / unit.nanometer**2)

                
for f in [bond_force, angle_force, torsion_force]:
    f.setForceGroup(0)
    system.addForce(f)
    
# =============================================================================
# Nonbonded interaction - BB–BB
# =============================================================================
# map bead names to type-indices (shared by all forces)
type_map = {"BBP":0, "BBO":1, "BBG":2, "HBG":3, "HBP":4}
M = 5

# ε & σ only for BB–BB
epsilon_bb   = np.zeros((M, M), dtype='float64')
sigma_bb     = np.zeros((M, M), dtype='float64')

epsilon_bb[0,0] = 1.0    
epsilon_bb[0,1] = 1.0   
epsilon_bb[0,2] = 1.0 
epsilon_bb[1,0] = 1.0   
epsilon_bb[2,0] = 1.0 

sigma_bb[0,0] = 1.0
sigma_bb[0,1] = 1.0
sigma_bb[0,2] = 1.0
sigma_bb[1,0] = 1.0
sigma_bb[2,0] = 1.0

epsilon_bb[1,1] = 1.0    
epsilon_bb[1,2] = 1.0 
epsilon_bb[2,1] = 1.0
sigma_bb[1,1] = 1.0   
sigma_bb[1,2] = 1.0
sigma_bb[2,1] = 1.0

epsilon_bb[2,2] = 1.0 
sigma_bb[2,2] =  1.0

eps_bb_list = epsilon_bb.ravel().tolist()
sig_bb_list = sigma_bb.ravel().tolist()

force_bb = openmm.CustomNonbondedForce(
    '4*eps*((sig/r)^12 - (sig/r)^6); '
    'eps=epsilon(type1,type2); sig=sigma(type1,type2)'
)
force_bb.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
force_bb.setCutoffDistance(1.122 * unit.nanometer)


force_bb.addTabulatedFunction('epsilon',
    openmm.Discrete2DFunction(M, M, eps_bb_list)
)
force_bb.addTabulatedFunction('sigma',
    openmm.Discrete2DFunction(M, M, sig_bb_list)
)
force_bb.addPerParticleParameter('type')


for atom in topology.atoms():
    t = type_map.get(atom.name, 0)
    force_bb.addParticle([t])
'''
for atom in topology.atoms():
    force_bb.addParticle([ type_map[atom.name] ])


# gather bonded‐pair list 
bonded_pairs = []
for idx in range(bond_force.getNumBonds()):
    i, j, _, _ = bond_force.getBondParameters(idx)
    bonded_pairs.append((int(i), int(j)))

force_bb.createExclusionsFromBonds(bonded_pairs, 3)
'''

# grab every directly‐bonded pair (atom1, atom2)
bonded_pairs = [(a.index, b.index) for a,b in topology.bonds()]

force_bb.createExclusionsFromBonds(bonded_pairs, 1)
force_bb.setForceGroup(1)
system.addForce(force_bb)



# =============================================================================
# Nonbonded interaction - BB–HB
# =============================================================================
epsilon_bbhb = np.zeros((M, M), dtype='float64')
sigma_bbhb   = np.zeros((M, M), dtype='float64')
epsilon_bbhb[0,3] = 1.0    # BB–HB
epsilon_bbhb[0,4] = 1.0
sigma_bbhb[0,3]   = 0.65
sigma_bbhb[0,4]   = 0.65

epsilon_bbhb[3,0] = 1.0    # BB–HB
epsilon_bbhb[4,0] = 1.0
sigma_bbhb[3,0]   = 0.65
sigma_bbhb[4,0]   = 0.65

epsilon_bbhb[1,3] = 1.0    # BB–HB
epsilon_bbhb[1,4] = 1.0
sigma_bbhb[1,3]   = 0.65
sigma_bbhb[1,4]   = 0.65

epsilon_bbhb[3,1] = 1.0    # BB–HB
epsilon_bbhb[4,1] = 1.0
sigma_bbhb[3,1]   = 0.65
sigma_bbhb[4,1]   = 0.65


epsilon_bbhb[2,3] = 1.0    # BB–HB
epsilon_bbhb[2,4] = 1.0
sigma_bbhb[2,3]   = 0.65
sigma_bbhb[2,4]   = 0.65

epsilon_bbhb[3,2] = 1.0    # BB–HB
epsilon_bbhb[4,2] = 1.0
sigma_bbhb[3,2]   = 0.65
sigma_bbhb[4,2]   = 0.65


eps_bbhb_list = epsilon_bbhb.ravel().tolist()
sig_bbhb_list = sigma_bbhb.ravel().tolist()

force_bbhb = openmm.CustomNonbondedForce(
    '4*eps*((sig/r)^12 - (sig/r)^6); '
    'eps=epsilon(type1,type2); sig=sigma(type1,type2)'
)
force_bbhb.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
force_bbhb.setCutoffDistance(0.729 * unit.nanometer)

force_bbhb.addTabulatedFunction('epsilon',
    openmm.Discrete2DFunction(M, M, eps_bbhb_list)
)
force_bbhb.addTabulatedFunction('sigma',
    openmm.Discrete2DFunction(M, M, sig_bbhb_list)
)
force_bbhb.addPerParticleParameter('type')


for atom in topology.atoms():
    t = type_map.get(atom.name, 0)
    force_bbhb.addParticle([t])
'''
for atom in topology.atoms():
    force_bbhb.addParticle([ type_map[atom.name] ])
'''
force_bbhb.createExclusionsFromBonds(bonded_pairs, 1)
force_bbhb.setForceGroup(2)
system.addForce(force_bbhb)



# =============================================================================
# Nonbonded interaction - HBD–HBD and HBA-HBA
# =============================================================================
epsilon_hbhb = np.zeros((M, M), dtype='float64')
sigma_hbhb   = np.zeros((M, M), dtype='float64')
epsilon_hbhb[3,3] = 1.0    
sigma_hbhb[3,3]   = 0.7
epsilon_hbhb[4,4] = 1.0    
sigma_hbhb[4,4]   = 0.7

eps_hbhb_list = epsilon_hbhb.ravel().tolist()
sig_hbhb_list = sigma_hbhb.ravel().tolist()

force_hbhb = openmm.CustomNonbondedForce(
    '4*eps*((sig/r)^12 - (sig/r)^6); '
    'eps=epsilon(type1,type2); sig=sigma(type1,type2)'
)
force_hbhb.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)
force_hbhb.setCutoffDistance(0.786 * unit.nanometer)

force_hbhb.addTabulatedFunction('epsilon',
    openmm.Discrete2DFunction(M, M, eps_hbhb_list)
)
force_hbhb.addTabulatedFunction('sigma',
    openmm.Discrete2DFunction(M, M, sig_hbhb_list)
)
force_hbhb.addPerParticleParameter('type')


for atom in topology.atoms():
    t = type_map.get(atom.name, 0)
    force_hbhb.addParticle([t])
'''

for atom in topology.atoms():
    force_hbhb.addParticle([ type_map[atom.name] ])

'''
force_hbhb.createExclusionsFromBonds(bonded_pairs, 1)
force_hbhb.setForceGroup(3)
system.addForce(force_hbhb)



# =============================================================================
# Define the HBD–HBA nonbonded interaction with switch
# =============================================================================
epsilon_hbda = np.zeros((M, M), dtype='float64')
sigma_hbda   = np.zeros((M, M), dtype='float64')
epsilon_hbda[3,4] = 269.5
sigma_hbda[3,4]   = 0.3
epsilon_hbda[4,3] = 269.5
sigma_hbda[4,3]   = 0.3

eps_hbda_list = epsilon_hbda.ravel().tolist()
sig_hbda_list = sigma_hbda.ravel().tolist()

force_hbda = openmm.CustomNonbondedForce(
    'eps*((sig/r)^12 - (sig/r)^6); '
    'eps=epsilon(type1,type2); sig=sigma(type1,type2)'
)
force_hbda.setNonbondedMethod(openmm.CustomNonbondedForce.CutoffPeriodic)

# define your cutoff
cutoff = 2.0 * unit.nanometer
force_hbda.setCutoffDistance(cutoff)

# smooth the potential between cutoffs
force_hbda.setUseSwitchingFunction(True)
force_hbda.setSwitchingDistance(cutoff - 1.9*unit.nanometer)

force_hbda.addTabulatedFunction('epsilon',
    openmm.Discrete2DFunction(M, M, eps_hbda_list)
)
force_hbda.addTabulatedFunction('sigma',
    openmm.Discrete2DFunction(M, M, sig_hbda_list)
)
force_hbda.addPerParticleParameter('type')

for atom in topology.atoms():
    t = type_map.get(atom.name, 0)
    force_hbda.addParticle([t])
'''    
for atom in topology.atoms():
    force_hbda.addParticle([ type_map[atom.name] ])
'''
force_hbda.createExclusionsFromBonds(bonded_pairs, 1)
force_hbda.setForceGroup(4)
system.addForce(force_hbda)

# =============================================================================
# Run Simulation
# =============================================================================
print("Periodic vectors:", system.getDefaultPeriodicBoxVectors())
platform = openmm.Platform.getPlatformByName("OpenCL")

#integrator = openmm.LangevinMiddleIntegrator(300*unit.kelvin, 0.01/unit.picosecond, 0.002*unit.picoseconds)
integrator = openmm.LangevinMiddleIntegrator(args.temp * unit.kelvin, 0.13 / unit.picosecond, 0.001546 * unit.picoseconds) 

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)
simulation.minimizeEnergy()
state = simulation.context.getState(getEnergy=True)
print("Minimized energy:", state.getPotentialEnergy())

#simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.context.setVelocitiesToTemperature(args.temp * unit.kelvin)

simulation.reporters.append(app.DCDReporter(f'{args.outprefix}.dcd', 1000, enforcePeriodicBox=True))
simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True, speed=True))

#simulation.step(5000000)  #7.73ns 
#simulation.step(6468305)   #10 ns
simulation.step(12936610)  #20 ns


with open(f'{args.outprefix}.pdb', 'w') as f:
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    topology.setPeriodicBoxVectors(state.getPeriodicBoxVectors())
    app.PDBFile.writeFile(topology, state.getPositions(), f)
