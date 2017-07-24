from openeye import oechem, oequacpac
from simtk.openmm import app
from simtk.openmm import Vec3
from simtk import unit
import itertools

proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS',
                   'LEU', 'MET', 'PRO', 'THR', 'TYR',
                   'ARG', 'ASP', 'GLN', 'GLY', 'ILE',
                   'LYS', 'PHE', 'SER', 'TRP', 'VAL']

rnaResidues = ['A', 'G', 'C', 'U', 'I']
dnaResidues = ['DA', 'DG', 'DC', 'DT', 'DI']


def oemol_to_openmmTop(mol):
    """
    This function converts an OEMol to an openmm topology
    The OEMol coordinates are assumed to be in Angstrom unit

    Parameters:
    -----------
    mol: OEMol molecule
        The molecule to convert

    Return:
    -------
    topology : OpenMM Topology
        The generated OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated with the
        generated topology in Angstrom units
    """
    # OE Hierarchical molecule view
    hv = oechem.OEHierView(mol, oechem.OEAssumption_BondedResidue + oechem.OEAssumption_ResPerceived + oechem.OEAssumption_PDBOrder)

    # Create empty OpenMM Topology
    topology = app.Topology()
    # Dictionary used to map oe atoms to openmm atoms
    oe_atom_to_openmm_at = {}

    for chain in hv.GetChains():

        # Create empty OpenMM Chain
        openmm_chain = topology.addChain(chain.GetChainID())

        for frag in chain.GetFragments():

            for hres in frag.GetResidues():

                # Get OE residue
                oe_res = hres.GetOEResidue()
                # Create OpenMM residue
                openmm_res = topology.addResidue(oe_res.GetName(), openmm_chain)

                for oe_at in hres.GetAtoms():
                    # Select atom element based on the atomic number
                    element = app.element.Element.getByAtomicNumber(oe_at.GetAtomicNum())
                    # Add atom OpenMM atom to the topology
                    openmm_at = topology.addAtom(oe_at.GetName(), element, openmm_res)
                    openmm_at.index = oe_at.GetIdx()
                    # Add atom to the mapping dictionary
                    oe_atom_to_openmm_at[oe_at] = openmm_at

    # Create bonds preserving the bond ordering
    for bond in mol.GetBonds():
        aromatic = None

        # Set the bond aromaticity
        if bond.IsAromatic():
            aromatic = 'Aromatic'

        topology.addBond(oe_atom_to_openmm_at[bond.GetBgn()], oe_atom_to_openmm_at[bond.GetEnd()],
                         type=aromatic, order=bond.GetOrder())

    dic = mol.GetCoords()
    positions = [Vec3(v[0], v[1], v[2]) for k, v in dic.items()] * unit.angstrom

    return topology, positions


def openmmTop_to_oemol(topology, positions):
    """
    This function converts an OpenMM topology in an OEMol

    Parameters:
    -----------
    topology : OpenMM Topology
        The OpenMM topology
    positions : OpenMM Quantity
        The molecule atom positions associated with the
        topology


    Return:
    -------
    oe_mol : OEMol
        The generated OEMol molecule
    """

    # Create an empty OEMol
    oe_mol = oechem.OEMol()

    # Mapping dictionary between openmm atoms and oe atoms
    openmm_atom_to_oe_atom = {}

    # Python set used to identify atoms that are not in protein residues
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues)

    for chain in topology.chains():
        for res in chain.residues():
            # Create an OEResidue
            oe_res = oechem.OEResidue()
            # Set OEResidue name
            oe_res.SetName(res.name)
            # If the atom is not a protein atom then set its heteroatom
            # flag to True
            if res.name not in keep:
                oe_res.SetFragmentNumber(chain.index + 1)
                oe_res.SetHetAtom(True)
            # Set OEResidue Chain ID
            oe_res.SetChainID(chain.id)
            # res_idx = int(res.id) - chain.index * len(chain._residues)
            # Set OEResidue number
            oe_res.SetResidueNumber(int(res.id))

            for openmm_at in res.atoms():
                # Create an OEAtom  based on the atomic number
                oe_atom = oe_mol.NewAtom(openmm_at.element._atomic_number)
                # Set atom name
                oe_atom.SetName(openmm_at.name)
                # Set Symbol
                oe_atom.SetType(openmm_at.element.symbol)
                # Set Atom index
                oe_res.SetSerialNumber(openmm_at.index + 1)
                # Commit the changes
                oechem.OEAtomSetResidue(oe_atom, oe_res)
                # Update the dictionary OpenMM to OE
                openmm_atom_to_oe_atom[openmm_at] = oe_atom

    # Create the bonds
    for bond in topology.bonds():
        at0 = bond[0]
        at1 = bond[1]
        # Read in the bond order from the OpenMM topology
        bond_order = bond.order

        # If bond order info are not present set the bond order to one
        if not bond_order:
            print("WARNING: Bond order info missing between atom indexes: {}-{}".format(at0.index, at1.index))
            bond_order = 1

        # OE atoms
        oe_atom0 = openmm_atom_to_oe_atom[at0]
        oe_atom1 = openmm_atom_to_oe_atom[at1]

        # Set OE atom aromaticity
        if bond.type:
            oe_atom0.SetAromatic(True)
            oe_atom1.SetAromatic(True)

        # Create the bond
        oe_bond = oe_mol.NewBond(oe_atom0, oe_atom1, bond_order)

        if bond.type:
            oe_bond.SetAromatic(True)

    # Set the OEMol positions
    pos = positions.in_units_of(unit.angstrom) / unit.angstrom
    pos = list(itertools.chain.from_iterable(pos))
    oe_mol.SetCoords(pos)

    return oe_mol


def delete_shell(core_mol, del_mol, cut_off, in_out='in'):
    """
    This function deletes molecules present in the passed argument
    del_mol that are far (in_out=out) or close (in_out=in) than the
    selected cutoff distance (in A) from the passed molecules core_mol

    Parameters:
    -----------
    core_mol: OEMol molecule
        The core molecules
    del_mol: OEMol molecule
        The molecules to be deleted if their distances from the core_mol
        molecules are greater or closer that the selected cutoff distance
    cut_off: python float number
        The threshold distance in A used to mark atom for deletion
    in_out: python string
        A flag used to select if delete molecules far or close than
        the cutoff distance from the core_mol

    Return:
    -------
    reset_del: copy of del_mol where atoms have been deleted with
        reset atom indexes
    """

    if in_out not in ['in', 'out']:
        raise ValueError("The passed in_out parameter is not recognized: {}".format(in_out))

    # Copy the passed molecule to delete in
    to_del = oechem.OEMol(del_mol)

    # Create a OE bit vector mask for each atoms of the
    # molecule to delete
    bv = oechem.OEBitVector(to_del.GetMaxAtomIdx())
    bv.NegateBits()

    # Create the Nearest neighbours
    nn = oechem.OENearestNbrs(to_del, cut_off)
    for nbrs in nn.GetNbrs(core_mol):
        # bv.SetBitOff(nbrs.GetBgn().GetIdx())
        for atom in oechem.OEGetResidueAtoms(nbrs.GetBgn()):
            bv.SetBitOff(atom.GetIdx())

    # Invert selection mask
    if in_out == 'in':
        bv.NegateBits()

    pred = oechem.OEAtomIdxSelected(bv)
    for atom in to_del.GetAtoms(pred):
        to_del.DeleteAtom(atom)

    # It is necessary to reset the atom indexes of the molecule with
    # delete atoms to avoid possible mismatching
    reset_del = oechem.OEMol(to_del)

    return reset_del


def check_shell(core_mol, check_mol, cutoff):
    """
    This function checks if at least one atomic distance from the passed
    check_mol molecule to the core_mol molecule is less than the selected
    cutoff distance in A.

    Parameters:
    -----------
    core_mol: OEMol molecule
        The core molecule
    check_mol: OEMol molecule
        The molecule to be checked if inside or outside a shell
        surrounding the core_mole with radius equal to the cutoff
        threshold
    cut_off: python float number
        The threshold distance in A used to mark atom inside or outside
        the shell

    Return:
    -------
    in_out: python boolean
         True if at least one of check_mol atom distance from core_mole
         is less than the selected cutoff threshold
    """

    # Create a OE bit vector mask for each atoms of the
    # molecule to be checked
    bv = oechem.OEBitVector(check_mol.GetMaxAtomIdx())

    # Create the Nearest neighbours
    nn = oechem.OENearestNbrs(check_mol, cutoff)

    # Check neighbours setting the atom bit mask
    for nbrs in nn.GetNbrs(core_mol):
        bv.SetBitOn(nbrs.GetBgn().GetIdx())

    # Create predicate based on the atom bit mask
    pred = oechem.OEAtomIdxSelected(bv)

    # Checking flag
    in_out = False

    # If just one chem_mol atom is inside the cutoff distance return True
    for atom in check_mol.GetAtoms(pred):
        in_out = True
        break

    return in_out


def sanitizeOEMolecule(molecule):
    """
    This function checks if the molecule has coordinates,
    explicit hydrogens, aromaticity missing and not unique
    atom names. If the molecule does not have coordinates
    a fatal error is raised. If the molecule does not have
    hydrogens or aramatic flags are missing then a copy of
    the molecule is fixed, if missing or not unique atom
    names are found then a copy of the molecule is fixed

    Parameters:
    -----------
    molecule: OEMol
        The molecule to be checked

    Return:
    -------
    mol_copy: OEMol
        A copy of the checked molecule with fixed aromaticity,
        hydrogens and unique atom names if they are missing
    """
    mol_copy = molecule.CreateCopy()

    # Check if the molecule has 3D coordinates
    if not oechem.OEGetDimensionFromCoords(mol_copy):
        oechem.OEThrow.Fatal("The molecule coordinates are set to zero")
    # Check if the molecule has hydrogens
    if not oechem.OEHasExplicitHydrogens(mol_copy):
        oechem.OEAddExplicitHydrogens(mol_copy)
    # Check if the molecule has assigned aromaticity
    if not mol_copy.HasPerceived(oechem.OEPerceived_Aromaticity):
        oechem.OEAssignAromaticFlags(mol_copy, oechem.OEAroModelOpenEye)

    # Check for any missing and not unique atom names.
    # If found reassign all of them as Tripos atom names

    atm_list_names = []

    for atom in mol_copy.GetAtoms():
        atm_list_names.append(atom.GetName())

    reassign_names = False

    if len(set(atm_list_names)) != len(atm_list_names):
        reassign_names = True

    if '' in atm_list_names:
        reassign_names = True

    if reassign_names:
        oechem.OETriposAtomNames(mol_copy)

    return mol_copy

def assignELF10charges(molecule, max_confs=800, strictStereo=True):
    """
     This function computes atomic partial charges for an OEMol by
     using the ELF10 method

    Parameters:
    -----------
    molecule : OEMol object
        The molecule that needs to be charged
    max_confs : integer
        The max number of conformers used to calculate the atomic partial charges
    strictStereo : bool
        a flag used to check if atoms need to have assigned stereo chemistry or not

    Return:
    -------
    mol_copy : OEMol
        a copy of the original molecule with assigned atomic partial charges
    """

    mol_copy = molecule.CreateCopy()

    # The passed molecule could have already conformers. If the conformer number
    # does not exceed the max_conf threshold then max_confs conformations will
    # be generated
    if not mol_copy.GetMaxConfIdx() > 200:
        # Generate up to max_confs conformers
        mol_copy = generate_conformers(mol_copy, max_confs=max_confs, strictStereo=strictStereo)

    # Assign MMFF Atom types
    if not oechem.OEMMFFAtomTypes(mol_copy):
        raise RuntimeError("MMFF atom type assignment returned errors")

    # ELF10 charges
    status = oequacpac.OEAssignCharges(mol_copy, oequacpac.OEAM1BCCELF10Charges())

    if not status:
        raise RuntimeError("OEAssignCharges returned error code %d" % status)

    return mol_copy
