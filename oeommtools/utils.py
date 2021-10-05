# (C) 2021 OpenEye Scientific Software Inc. All rights reserved.
#
# TERMS FOR USE OF SAMPLE CODE The software below ("Sample Code") is
# provided to current licensees or subscribers of OpenEye products or
# SaaS offerings (each a "Customer").
# Customer is hereby permitted to use, copy, and modify the Sample Code,
# subject to these terms. OpenEye claims no rights to Customer's
# modifications. Modification of Sample Code is at Customer's sole and
# exclusive risk. Sample Code may require Customer to have a then
# current license or subscription to the applicable OpenEye offering.
# THE SAMPLE CODE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED.  OPENEYE DISCLAIMS ALL WARRANTIES, INCLUDING, BUT
# NOT LIMITED TO, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
# PARTICULAR PURPOSE AND NONINFRINGEMENT. In no event shall OpenEye be
# liable for any damages or liability in connection with the Sample Code
# or its use.

from openeye import oechem
from simtk.openmm import app
from simtk.openmm import Vec3
from simtk import unit
import itertools
import pyparsing as pyp


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
    molecule: OEMol molecule
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
    hv = oechem.OEHierView(mol,
                           oechem.OEAssumption_ResPerceived +
                           oechem.OEAssumption_PDBOrder)

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
                # print(oe_res.GetName(), oe_res.GetResidueNumber(), oe_res.GetChainID())
                for oe_at in hres.GetAtoms():
                    # print("\t", oe_at.GetName(), oe_res.GetSerialNumber(), oe_at.GetIdx())
                    # Select atom element based on the atomic number
                    element = app.element.Element.getByAtomicNumber(oe_at.GetAtomicNum())
                    # Add atom OpenMM atom to the topology
                    openmm_at = topology.addAtom(oe_at.GetName(), element, openmm_res)
                    # Add atom to the mapping dictionary
                    oe_atom_to_openmm_at[oe_at] = openmm_at

    if topology.getNumAtoms() != mol.NumAtoms():
        raise ValueError("OpenMM topology and OEMol number of atoms mismatching: "
                         "OpenMM = {} vs OEMol  = {}".format(topology.getNumAtoms(), mol.NumAtoms()))

    # Count the number of bonds in the openmm topology
    omm_bond_count = 0

    def IsAmideBond(oe_bond):

        # This supporting function checks if the passed bond is an amide bond or not.
        # Our definition of amide bond C-N between a Carbon and a Nitrogen atom is:
        #          O
        #          ║
        #  CA or O-C-N-
        #            |

        # The amide bond C-N is a single bond
        if oe_bond.GetOrder() != 1:
            return False

        atomB = oe_bond.GetBgn()
        atomE = oe_bond.GetEnd()

        # The amide bond is made by Carbon and Nitrogen atoms
        if not (atomB.IsCarbon() and atomE.IsNitrogen() or
                (atomB.IsNitrogen() and atomE.IsCarbon())):
            return False

        # Select Carbon and Nitrogen atoms
        if atomB.IsCarbon():
            C_atom = atomB
            N_atom = atomE
        else:
            C_atom = atomE
            N_atom = atomB

        # Carbon and Nitrogen atoms must have 3 neighbour atoms
        if not (C_atom.GetDegree() == 3 and N_atom.GetDegree() == 3):
            return False

        double_bonds = 0
        single_bonds = 0

        for bond in C_atom.GetBonds():
            # The C-O bond can be single or double.
            if (bond.GetBgn() == C_atom and bond.GetEnd().IsOxygen()) or \
                    (bond.GetBgn().IsOxygen() and bond.GetEnd() == C_atom):
                if bond.GetOrder() == 2:
                    double_bonds += 1
                if bond.GetOrder() == 1:
                    single_bonds += 1
            # The CA-C bond is single
            if (bond.GetBgn() == C_atom and bond.GetEnd().IsCarbon()) or \
                    (bond.GetBgn().IsCarbon() and bond.GetEnd() == C_atom):
                if bond.GetOrder() == 1:
                    single_bonds += 1
        # Just one double and one single bonds are connected to C
        # In this case the bond is an amide bond
        if double_bonds == 1 and single_bonds == 1:
            return True
        else:
            return False

    # Creating bonds
    bond_mapping = {"Single": app.Single,
                    "Double": app.Double,
                    "Triple": app.Triple,
                    "Amide": app.Amide,
                    "Aromatic": app.Aromatic}

    for oe_bond in mol.GetBonds():

        omm_bond_count += 1

        # Set the bond type
        if oe_bond.GetType() is not "":
            if oe_bond.GetType() in bond_mapping:
                omm_bond_type = bond_mapping[oe_bond.GetType()]
            else:
                omm_bond_type = None
        else:
            if oe_bond.IsAromatic():
                oe_bond.SetType("Aromatic")
                omm_bond_type = app.Aromatic
            elif oe_bond.GetOrder() == 2:
                oe_bond.SetType("Double")
                omm_bond_type = app.Double
            elif oe_bond.GetOrder() == 3:
                oe_bond.SetType("Triple")
                omm_bond_type = app.Triple
            elif IsAmideBond(oe_bond):
                oe_bond.SetType("Amide")
                omm_bond_type = app.Amide
            elif oe_bond.GetOrder() == 1:
                oe_bond.SetType("Single")
                omm_bond_type = app.Single
            else:
                omm_bond_type = None

        topology.addBond(oe_atom_to_openmm_at[oe_bond.GetBgn()], oe_atom_to_openmm_at[oe_bond.GetEnd()],
                         type=omm_bond_type, order=oe_bond.GetOrder())

    if omm_bond_count != mol.NumBonds():
        raise ValueError("OpenMM topology and OEMol number of bonds mismatching: "
                         "OpenMM = {} vs OEMol  = {}".format(omm_bond_count, mol.NumBonds()))
    dic = mol.GetCoords()
    positions = [Vec3(v[0], v[1], v[2]) for k, v in dic.items()] * unit.angstrom

    return topology, positions


def openmmTop_to_oemol(topology, positions, verbose=False):
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

    if topology.getNumAtoms() != oe_mol.NumAtoms():
        raise ValueError("OpenMM topology and OEMol number of atoms mismatching: "
                         "OpenMM = {} vs OEMol  = {}".format(topology.getNumAtoms(), oe_mol.NumAtoms()))

    # Count the number of bonds in the openmm topology
    omm_bond_count = 0

    # Create the bonds
    bond_mapping = {app.Single: "Single",
                    app.Double: "Double",
                    app.Triple: "Triple",
                    app.Amide: "Amide",
                    app.Aromatic: "Aromatic"}

    for omm_bond in topology.bonds():

        omm_bond_count += 1

        at0 = omm_bond[0]
        at1 = omm_bond[1]

        oe_bond_order = omm_bond.order

        # If bond order info are not present set the bond order temporary to one
        if not omm_bond.order:
            oe_bond_order = 1

        # OE atoms
        oe_atom0 = openmm_atom_to_oe_atom[at0]
        oe_atom1 = openmm_atom_to_oe_atom[at1]

        # Create the bond
        oe_bond = oe_mol.NewBond(oe_atom0, oe_atom1, oe_bond_order)

        if omm_bond.type:
            if omm_bond.type == app.Aromatic:
                oe_atom0.SetAromatic(True)
                oe_atom1.SetAromatic(True)
                oe_bond.SetAromatic(True)
                oe_bond.SetType(bond_mapping[omm_bond.type])
            elif omm_bond.type in bond_mapping:
                oe_bond.SetType(bond_mapping[omm_bond.type])
            else:
                oe_bond.SetType("")

    if omm_bond_count != oe_mol.NumBonds():
        raise ValueError("OpenMM topology and OEMol number of bonds mismatching: "
                         "OpenMM = {} vs OEMol  = {}".format(omm_bond_count, oe_mol.NumBonds()))

    # Set the OEMol positions
    pos = positions.in_units_of(unit.angstrom) / unit.angstrom
    pos = list(itertools.chain.from_iterable(pos))
    oe_mol.SetCoords(pos)
    oechem.OESetDimensionFromCoords(oe_mol)

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
    mol_copy = oechem.OEMol(molecule)

    # Check if the molecule has 3D coordinates
    if not mol_copy.NumAtoms() == 1:  # Mono-atomic molecules are skipped
        if not oechem.OEGetDimensionFromCoords(mol_copy):
            raise ValueError("The molecule coordinates are set to zero")
    # Check if the molecule has hydrogens
    if not oechem.OEHasExplicitHydrogens(mol_copy):
        oechem.OEAddExplicitHydrogens(mol_copy)
    # Check if the molecule has assigned aromaticity
    if not mol_copy.HasPerceived(oechem.OEPerceived_Aromaticity):
        # oechem.OEAssignAromaticFlags(mol_copy, oechem.OEAroModelOpenEye)
        oechem.OEAssignAromaticFlags(mol_copy, oechem.OEAroModelMDL)
    if not mol_copy.HasPerceived(oechem.OEPerceived_Chiral):
        oechem.OEPerceiveChiral(mol_copy)

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


def strip_water_ions(in_system):
    """
    This function remove waters and ions molecules
    from the input system

    Parameters:
    ----------
    in_system : oechem.OEMol
        The bio-molecular system to clean
    opt: python dictionary
        The system option

    Output:
    -------
    clean_system : oechem.OEMol
        The cleaned system

    """
    # Copy the input system
    system = oechem.OEMol(in_system)

    # Create a bit vector mask
    bv = oechem.OEBitVector(system.GetMaxAtomIdx())
    bv.NegateBits()

    # Create a Hierarchical View of the protein system
    hv = oechem.OEHierView(system, oechem.OEAssumption_BondedResidue +
                           oechem.OEAssumption_ResPerceived)

    # Looping over the system residues
    for chain in hv.GetChains():
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                res = hres.GetOEResidue()

                # Check if a residue is a mono atomic ion
                natoms = 0
                for at in hres.GetAtoms():
                    natoms += 1

                # Set the atom bit mask off
                if oechem.OEGetResidueIndex(res) == oechem.OEResidueIndex_HOH or natoms == 1:
                    # Set Bit mask
                    atms = hres.GetAtoms()
                    for at in atms:
                        bv.SetBitOff(at.GetIdx())

    # Extract the system without waters or ions
    pred = oechem.OEAtomIdxSelected(bv)
    clean_system = oechem.OEMol()
    oechem.OESubsetMol(clean_system, system, pred)

    return clean_system


def select_oemol_atom_idx_by_language(system, mask=''):
    """
    This function selects the atom indexes from the passed oemol molecular complex
    by using  a defined language. The language allows the selection of the ligand,
    protein, waters, mono-atomic ions, excipients, residue numbers and distance selection. Logic
    operators not, or, and, noh, diff, around can be used to refine the selection

    Parameters
    ----------
    system : OEMol of the bio-molecular complex protein-ligand
        The molecular complex

    mask : python string
        A string used to select atoms. A Backus–Naur Form grammar
        (https://en.wikipedia.org/wiki/Backus–Naur_form) is defined by the python
        module pyparsing.
        The defined grammar tokens are: "ligand", "protein", "ca_protein" ,"water",
        "ions", "excipients" and "resid chain1:res_idx1 chain2:res_idx2 ... res_idxn"
        that respectively define the ligand, the protein, carbon alpha protein atoms,
        water molecules, ions, excipients (not protein, ligand, water or ions) and
        residue numbers. The atom selection can be refined by using the following
        operator tokens:

        "not" = invert selection
        "or" = add selections
        "and" = intersect selections
        "diff" = logic difference between selections
        "noh" = remove hydrogens from the selection
        "around" = select atoms inside the cutoff distance from a given selection

    Returns
    -------
    atom_set : python set
        the select atom indexes

    Notes
    -----
        Example of selection string:
        mask = "ligand or protein"
        mask = "not water or not ions"
        mask = "ligand or protein or excipients"
        mask = "noh protein"
        mask = "resid A:17 B:12 17 18"
        mask = "protein diff resid A:1"
        mask = "5.0 around protein"
    """

    def split(system, ligand_res_name='LIG'):
        """
        This function splits the passed molecule in components and tracks the
        mapping between the original molecule and the split components. The
        mapping is created as separated atom component index sets.

        Parameters:
        -----------
        system: OEMol
            The system to split in components. The components are:
                the protein atoms,
                the protein carbon alpha atoms
                the water atoms,
                the ion atoms,
                the excipients atoms
        Returns:
        --------
        dic_set: python dictionary
            The sysetm is splitted in a dictionary with token words as keys
            and for value the related atom set. The token keywords are:
                protein,
                ca_protein,
                ligand,
                water,
                ions,
                excipients,
                system
        """

        # Define Empty sets
        lig_set = set()
        prot_set = set()
        ca_prot_set = set()
        wat_set = set()
        excp_set = set()
        ion_set = set()

        # Atom Bond Set vector used to contains the whole system
        frags = oechem.OEAtomBondSetVector()

        # Define Options for the Filter
        opt = oechem.OESplitMolComplexOptions()

        # The protein filter is set to avoid that multiple
        # chains are separated during the splitting and peptide
        # molecules are recognized as ligands
        pf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Protein)
        peptide = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Peptide)
        protein_filter = oechem.OEOrRoleSet(pf, peptide)
        opt.SetProteinFilter(protein_filter)

        # The ligand filter is set to recognize just the ligand
        lf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Ligand)
        not_protein_filter = oechem.OENotRoleSet(protein_filter)
        ligand_filter = oechem.OEAndRoleSet(lf, not_protein_filter)
        opt.SetLigandFilter(ligand_filter)

        # The water filter is set to recognize just water molecules
        wf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Water)
        opt.SetWaterFilter(wf)

        # Set Category
        cat = oechem.OEMolComplexCategorizer()
        cat.AddLigandName(ligand_res_name)
        opt.SetCategorizer(cat)

        # Define the system fragments
        if not oechem.OEGetMolComplexFragments(frags, system, opt):
            raise ValueError('Unable to generate the system fragments')

        # Set empty OEMol containers
        prot = oechem.OEMol()
        lig = oechem.OEMol()
        wat = oechem.OEMol()
        excp = oechem.OEMol()

        # Split the protein from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(prot, frags, opt, opt.GetProteinFilter(), atommap):
            raise ValueError('Unable to split the Protein')
        # Populate the protein set and the protein carbon alpha set
        pred = oechem.OEIsCAlpha()
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                prot_set.add(sys_idx)
                at = system.GetAtom(oechem.OEHasAtomIdx(sys_idx))
                if pred(at):
                    ca_prot_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Split the ligand from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(lig, frags, opt, opt.GetLigandFilter(), atommap):
            raise ValueError('Unable to split the Ligand')
        # Populate the ligand set
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                lig_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Split the water from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(wat, frags, opt, opt.GetWaterFilter(), atommap):
            raise ValueError('Unable to split the Water')
        # Populate the water set
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                wat_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Split the excipients from the system
        atommap = oechem.OEAtomArray(system.GetMaxAtomIdx())
        if not oechem.OECombineMolComplexFragments(excp, frags, opt, opt.GetOtherFilter(), atommap):
            raise ValueError('Unable to split the Excipients')
        # Populate the excipient set
        for sys_at in system.GetAtoms():
            sys_idx = sys_at.GetIdx()
            at_idx = atommap[sys_idx]
            if at_idx:
                excp_set.add(sys_idx)
                # print(sys_idx, '->', at_idx)

        # Create the mono-atomic ions set
        for exc_idx in excp_set:
            atom = system.GetAtom(oechem.OEHasAtomIdx(exc_idx))
            if atom.GetDegree() == 0:
                ion_set.add(exc_idx)

        # Create the excipients set which are not protein, ligand, waters or ions
        excipients_set = excp_set - ion_set

        # Create the system set
        system_set = prot_set | lig_set | excp_set | wat_set

        if len(system_set) != system.NumAtoms():
            raise ValueError("The total system atom number {} is different "
                             "from its set representation {}".format(system.NumAtoms(), system_set))

        # The dictionary is used to link the token keywords to the created molecule sets
        dic_set = {'ligand': lig_set, 'protein': prot_set, 'ca_protein': ca_prot_set,
                   'water': wat_set,  'ions': ion_set,     'excipients': excipients_set, 'system': system_set}

        return dic_set

    def build_set(ls, dsets):
        """
        This function select the atom indexes

        Parameters:
        -----------
        ls: python list
            the parsed list with tokens and operand tokes for the selection
        dsets: python dictionary
             the dictionary containing the sets for the selection

        Return:
        -------
        atom_set: python set
            the set containing the atom index
        """

        def noh(ls, dsets):
            """
            This function remove hydrogens from the selection
            """
            data_set = build_set(ls[1], dsets)

            noh_set = set()
            pred = oechem.OEIsHydrogen()

            for idx in data_set:
                atom = system.GetAtom(oechem.OEHasAtomIdx(idx))
                if not pred(atom):
                    noh_set.add(idx)

            return noh_set

        def residues(ls):
            """
            This function select residues based on the residue numbers. An example of
            selection can be:
            mask = 'resid A:16 17 19 B:1'
            """
            # List residue atom index to be restrained
            res_atom_set = set()

            # Dictionary of lists with the chain residues selected to be restrained
            # e.g. {chainA:[res1, res15], chainB:[res19, res17]}
            chain_dic = {'': []}

            # Fill out the chain dictionary
            i = 0
            while i < len(ls):
                if ls[i].isdigit():
                    chain_dic[''].append(int(ls[i]))
                    i += 1
                else:
                    try:
                        chain_dic[ls[i]].append(int(ls[i + 2]))
                    except:
                        chain_dic[ls[i]] = []
                        chain_dic[ls[i]].append(int(ls[i + 2]))
                    i += 3

            # Loop over the molecular system to select the atom indexes to be selected
            hv = oechem.OEHierView(system, oechem.OEAssumption_BondedResidue + oechem.OEAssumption_ResPerceived)
            for chain in hv.GetChains():
                chain_id = chain.GetChainID()
                if chain_id not in chain_dic:
                    continue
                for frag in chain.GetFragments():
                    for hres in frag.GetResidues():
                        res_num = hres.GetOEResidue().GetResidueNumber()
                        if res_num not in chain_dic[chain_id]:
                            continue
                        for oe_at in hres.GetAtoms():
                            res_atom_set.add(oe_at.GetIdx())

            return res_atom_set

        def around(dist, ls):
            """
            This function select atom not far than the threshold distance from
            the current selection. The threshold distance is in Angstrom

            selection can be:
            mask = '5.0 around ligand'
            """
            # at = system.GetAtom(oechem.OEHasAtomIdx(idx))

            # Atom set selection
            atom_set_around = set()

            # Create a OE bit vector mask for each atoms
            bv_around = oechem.OEBitVector(system.GetMaxAtomIdx())

            # Set the mask atom
            for at in system.GetAtoms():
                if at.GetIdx() in ls:
                    bv_around.SetBitOn(at.GetIdx())

            # Predicate
            pred = oechem.OEAtomIdxSelected(bv_around)

            # Create the system molecule based on the atom mask
            molecules = oechem.OEMol()
            oechem.OESubsetMol(molecules, system, pred)

            # Create the Nearest neighbours
            nn = oechem.OENearestNbrs(system, float(dist))

            for nbrs in nn.GetNbrs(molecules):
                for atom in oechem.OEGetResidueAtoms(nbrs.GetBgn()):
                    if atom.GetIdx() in ls:
                        continue
                    atom_set_around.add(atom.GetIdx())

            return atom_set_around

        # Start Body of the selection function by language

        # Terminal Literal return the related set
        if isinstance(ls, str):
            return dsets[ls]
        # Not or Noh
        if len(ls) == 2:
            if ls[0] == 'noh':  # Noh case
                return noh(ls, dsets)
            elif ls[0] == 'not':  # Not case
                return dsets['system'] - build_set(ls[1], dsets)
            else:  # Resid case with one index
                return residues(ls[1])

        if len(ls) == 3:
            if ls[1] == 'or':  # Or Case (set union)
                return build_set(ls[0], dsets) | build_set(ls[2], dsets)
            elif ls[1] == 'and':  # And Case (set intersection)
                return build_set(ls[0], dsets) & build_set(ls[2], dsets)
            elif ls[1] == 'diff':  # Diff case (set difference)
                return build_set(ls[0], dsets) - build_set(ls[2], dsets)
            elif ls[1] == 'around':  # Around case
                return around(ls[0], build_set(ls[2], dsets))
            else:
                return residues(ls[1:])  # Resid case with one or two indexes
        else:
            if ls[0] == 'resid':
                return residues(ls[1:])  # Resid case with multiple indexes
            else:
                raise ValueError("The passed list have too many tokens: {}".format(ls))

    # Parse Action-Maker
    def makeLRlike(numterms):
        if numterms is None:
            # None operator can only by binary op
            initlen = 2
            incr = 1
        else:
            initlen = {0: 1, 1: 2, 2: 3, 3: 5}[numterms]
            incr = {0: 1, 1: 1, 2: 2, 3: 4}[numterms]

        # Define parse action for this number of terms,
        # to convert flat list of tokens into nested list
        def pa(s, l, t):
            t = t[0]
            if len(t) > initlen:
                ret = pyp.ParseResults(t[:initlen])
                i = initlen
                while i < len(t):
                    ret = pyp.ParseResults([ret] + t[i:i + incr])
                    i += incr
                return pyp.ParseResults([ret])

        return pa

    # Selection function body

    # Residue number selection
    id = pyp.Optional(pyp.Word(pyp.alphanums) + pyp.Literal(':')) + pyp.Word(pyp.nums)
    resid = pyp.Group(pyp.Literal("resid") + pyp.OneOrMore(id))

    # Real number for around operator selection
    real = pyp.Regex(r"\d+(\.\d*)?").setParseAction(lambda t: float(t[0]))

    # Define the tokens for the BNF grammar
    operand = pyp.Literal("protein") | pyp.Literal("ca_protein") | \
              pyp.Literal("ligand") | pyp.Literal("water") | \
              pyp.Literal("ions") | pyp.Literal("excipients") | resid

    # BNF Grammar definition with parseAction makeLRlike
    expr = pyp.operatorPrecedence(operand,
                                    [
                                        (None, 2, pyp.opAssoc.LEFT, makeLRlike(None)),
                                        (pyp.Literal("not"), 1, pyp.opAssoc.RIGHT, makeLRlike(1)),
                                        (pyp.Literal("noh"), 1, pyp.opAssoc.RIGHT, makeLRlike(1)),
                                        (pyp.Literal("and"), 2, pyp.opAssoc.LEFT, makeLRlike(2)),
                                        (pyp.Literal("or"), 2, pyp.opAssoc.LEFT, makeLRlike(2)),
                                        (pyp.Literal("diff"), 2, pyp.opAssoc.LEFT, makeLRlike(2)),
                                        (real + pyp.Literal("around"), 1, pyp.opAssoc.RIGHT, makeLRlike(2))
                                    ])
    # Parse the input string
    try:
        ls = expr.parseString(mask, parseAll=True)
    except Exception as e:
        raise ValueError("The passed restraint mask is not valid: {}".format(str(e)))

    # Split the system
    dic_sets = split(system)

    # Select atom indexes
    atom_set = build_set(ls[0], dic_sets)

    return atom_set


def split(complex, ligand_res_name='LIG'):
    """
    This function splits the passed system in protein, ligand,
    water and excipients.

    Parameters:
    ----------
    complex : oechem.OEMol
        The bio-molecular complex to split
    ligand_res_name : Python string
        The ligand residue name used to recognize the ligand

    Output:
    -------
    protein : oechem.OEMol
        The split protein
    ligand : oechem.OEMol
        The split ligand
    wat : oechem.OEMol
        The spit water
    cofactors : oechem.OEMol
        The cofactors
    lipids : oechem.OEMol
        The lipids
    metals : oechem.OEMol
        The metals
    excipients : oechem.OEMol
        The excipients
    """

    # Set empty molecule containers
    protein = oechem.OEMol()
    ligand = oechem.OEMol()
    water = oechem.OEMol()
    other = oechem.OEMol()

    # Define the Filter options before the splitting
    opt = oechem.OESplitMolComplexOptions()

    # The protein filter is set to avoid that multiple
    # chains are separated during the splitting and peptide
    # molecules are recognized as ligands
    pf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Protein)
    peptide = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Peptide)
    protein_filter = oechem.OEOrRoleSet(pf, peptide)
    opt.SetProteinFilter(protein_filter)

    # The ligand filter is set to recognize just the ligand
    lf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Ligand)
    not_protein_filter = oechem.OENotRoleSet(protein_filter)
    ligand_filter = oechem.OEAndRoleSet(lf, not_protein_filter)
    opt.SetLigandFilter(ligand_filter)

    # The water filter is set to recognize just water molecules
    wf = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Water)
    opt.SetWaterFilter(wf)

    # Set Category
    cat = oechem.OEMolComplexCategorizer()
    cat.AddLigandName(ligand_res_name)
    opt.SetCategorizer(cat)

    # Splitting the system
    if not oechem.OESplitMolComplex(ligand, protein, water, other, complex, opt):
        raise ValueError('Unable to split the complex')

    cofactors = oechem.OEMol()
    lipids = oechem.OEMol()
    metals = oechem.OEMol()
    excipients = oechem.OEMol()

    # Splitting the subsystem
    if other.NumAtoms():

        opt_other = oechem.OESplitMolComplexOptions()

        lipid_filter = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Lipid)
        opt_other.SetProteinFilter(lipid_filter)

        cofactor_filter = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Cofactor)
        opt_other.SetLigandFilter(cofactor_filter)

        metal_filter = oechem.OEMolComplexFilterFactory(oechem.OEMolComplexFilterCategory_Metal)
        opt_other.SetWaterFilter(metal_filter)

        # Set Category
        cat_other = oechem.OEMolComplexCategorizer()
        opt_other.SetCategorizer(cat_other)

        if not oechem.OESplitMolComplex(cofactors, lipids, metals, excipients, other, opt_other):
            raise ValueError('Unable to split the Subcategories')

    return protein, ligand, water, cofactors, lipids, metals, excipients
