from openeye import oechem
from openeye import oeomega
import numpy as np
import shutil
from simtk import unit
from simtk.openmm import Vec3
import subprocess
import tempfile
import os
import oeommtools
from oeommtools import data_utils
import random
import string

PACKAGE_DIR = os.path.dirname(os.path.dirname(oeommtools.__file__))


def oesolvate(solute, density=1.0, padding_distance=10.0,
              distance_between_atoms=2.5,
              solvents='tip3p', molar_fractions='1.0',
              geometry='box', close_solvent=True,
              salt='[Na+], [Cl-]', salt_concentration=0.0,
              neutralize_solute=True, verbose=False, return_components=False, **kargs):
    """
    This function solvates the passed solute in a cubic box or a sphere by using Packmol. Packmol
    creates an initial point for molecular dynamics simulations by packing molecule in defined regions
    of space. For additional info:
    http://www.ime.unicamp.br/~martinez/packmol/home.shtml

    The geometry volume is estimated by the using the padding parameter and the solute size.
    The number of solvent molecules is calculated by using the specified density and volume.
    Solvent molecules are specified as comma separated smiles strings. The molar fractions
    of each solvent molecule are specified in a similar fashion. By default if the solute is
    charged counter ions are added to neutralize it

    Parameters:
    -----------
    solute: OEMol molecule
        The solute to solvate
    density: float
        The solution density in g/ml
    padding_distance: float
        The largest dimension of the solute (along the x, y, or z axis) is determined (in A), 
        and a cubic box of size (largest dimension)+2*padding is used
    distance_between_atoms: float
        The minimum distance between atoms in A
    solvents: python string
        A comma separated smiles string or keywords for the solvent molecules.
        Special water models can be selected by using the keywords:
        tip3p for TIP3P water model geometry
    molar_fractions: python string
        A comma separated molar fraction string of the solvent molecules
    close_solvent: boolean
        If True solvent molecules will be placed very close to the solute
    salt: python string
        A comma separated string of the dissociated salt in solution
    salt_concentration: float
        Salt concentration in millimolar
    neutralize_solute: boolean
        If True counter-ions will be added to the solution to neutralize the solute
    verbose: Bool
        If True verbose mode is enabled
    return_components: Bool
        If True the added solvent molecule components are also returned as OEMols

    Return:
    -------
    oe_mol: OEMol
        The solvated system. If the selected geometry is a box a SD tag with
        name 'box_vector' is attached the output molecule containing
        the system box vectors.
    oe_mol_components: multiple OEMol
        If the return_components flag is True the added solvent molecules components
        are also returned as additional OEMols in addition to the whole solvated system
        If some components is missing None is return for that component
    """

    def BoundingBox(molecule):
        """
        This function calculates the Bounding Box of the passed
        molecule

        molecule: OEMol

        return: bb (numpy array)
            the calculated bounding box is returned as numpy array:
            [(xmin,ymin,zmin), (xmax,ymax,zmax)]
        """
        coords = [v for k, v in molecule.GetCoords().items()]
        np_coords = np.array(coords)
        min_coord = np_coords.min(axis=0)
        max_coord = np_coords.max(axis=0)
        bb = np.array([min_coord, max_coord])
        return bb

    if shutil.which("packmol") is None:
        raise (IOError("Packmol executable not found"))

    # Extract solvent smiles strings and mole fractions
    solvents = [sm.strip() for sm in solvents.split(',')]
    fractions = [float(mf) for mf in molar_fractions.split(',')]

    # If the smiles string and mole fractions lists have different lengths raise an error
    if len(solvents) != len(fractions):
        raise ValueError("Selected solvent number and selected molar fraction number mismatch: {} vs {}"
                         .format(len(solvents), len(fractions)))

    # Remove smiles string with 0.0 mole fraction
    solvent_smiles = [solvents[i] for i, v in enumerate(fractions) if fractions[i]]
    mol_fractions = [mf for mf in fractions if mf]

    # Mole fractions are non-negative numbers
    if any([v < 0.0 for v in mol_fractions]):
        raise ValueError("Error: Mole fractions are non-negative real numbers")

    # Mole fractions must sum up to 1.0
    if abs(sum(mol_fractions) - 1.0) > 0.001:
        oechem.OEThrow.Error("Error: Mole fractions do not sum up to 1.0")

    if geometry not in ['box', 'sphere']:
        raise ValueError("Error geometry: the supported geometries are box and sphere not {}".format(geometry))

    # Set Units
    density = density * unit.grams/unit.milliliter
    padding_distance = padding_distance * unit.angstrom
    salt_concentration = salt_concentration * unit.millimolar

    solute_copy = oechem.OEMol(solute)

    # Calculate the Solute Bounding Box
    BB_solute = BoundingBox(solute_copy)

    # Estimate of the box cube length
    box_edge = 2.0*padding_distance + np.max(BB_solute[1] - BB_solute[0])*unit.angstrom

    if geometry == 'box':
        # Box Volume
        Volume = box_edge**3
    if geometry == 'sphere':
        Volume = (4.0/3.0) * 3.14159265 * (0.5*box_edge)**3

    # Omega engine is used to generate conformations
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(1)
    omegaOpts.SetStrictStereo(False)
    omega = oeomega.OEOmega(omegaOpts)

    # Create a string code to identify the solute residues. The code ID used is based
    # on the residue number id, the residue name and the chain id:
    # id+resname+chainID
    hv_solute = oechem.OEHierView(solute_copy,
                                  oechem.OEAssumption_BondedResidue +
                                  oechem.OEAssumption_ResPerceived +
                                  oechem.OEAssumption_PDBOrder)

    solute_resid_list = []
    for chain in hv_solute.GetChains():
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                oe_res = hres.GetOEResidue()
                solute_resid_list.append(str(oe_res.GetResidueNumber())+oe_res.GetName()+chain.GetChainID())

    # Solvent component list_names
    solvent_resid_dic_names = dict()

    # Neutralize solute
    ion_sum_wgt_n_ions = 0.0*unit.grams/unit.mole
    if neutralize_solute:
        # Container for the counter-ions
        oe_ions = []
        # Container for the ion smiles strings
        ions_smiles = []
        solute_formal_charge = 0
        for at in solute_copy.GetAtoms():
            solute_formal_charge += at.GetFormalCharge()
        if solute_formal_charge > 0:
            ions_smiles.append("[Cl-]")
        elif solute_formal_charge < 0:
            ions_smiles.append("[Na+]")
        else:
            pass

        # Total number of counter-ions to neutralize the solute
        n_ions = abs(solute_formal_charge)

        # print("Counter ions to add = {} of {}".format(n_ions, ions_smiles[0]))

        # Ions
        if n_ions >= 1:
            for sm in ions_smiles:
                mol = oechem.OEMol()
                if not oechem.OESmilesToMol(mol, sm):
                    raise ValueError("Error counter ions: SMILES string parsing fails for the string: {}".format(sm))

                # Generate conformer
                if not omega(mol):
                    raise ValueError("Error counter ions: Conformer generation fails for the molecule with "
                                     "smiles string: {}".format(sm))

                oe_ions.append(mol)

                if sm == '[Na+]':
                    solvent_resid_dic_names[' NA'] = mol
                else:
                    solvent_resid_dic_names[' CL'] = mol

            ion_sum_wgt = 0.0 * unit.grams / unit.mole
            for ion in oe_ions:
                # Molecular weight
                ion_sum_wgt += oechem.OECalculateMolecularWeight(ion) * unit.grams / unit.mole

            ion_sum_wgt_n_ions = ion_sum_wgt * n_ions

            # Create ions .pdb files
            ions_smiles_pdbs = []
            for i in range(0, len(ions_smiles)):
                pdb_name = os.path.basename(tempfile.mktemp(suffix='.pdb'))
                pdb_name = ions_smiles[i] + '_' + pdb_name
                ions_smiles_pdbs.append(pdb_name)

            for i in range(0, len(ions_smiles)):
                ofs = oechem.oemolostream(ions_smiles_pdbs[i])
                oechem.OEWriteConstMolecule(ofs, oe_ions[i])

    # Add salts to the solution

    # Solvent smiles string parsing
    char_set = string.ascii_uppercase
    salt_sum_wgt_n_salt = 0.0 * unit.grams/unit.mole
    if salt_concentration > 0.0*unit.millimolar:

        salt_smiles = [sm.strip() for sm in salt.split(',')]

        # Container list of oemol salt molecules generated by using smiles strings
        oe_salt = []

        for sm in salt_smiles:
            mol_salt = oechem.OEMol()
            if not oechem.OESmilesToMol(mol_salt, sm):
                raise ValueError("Error salt: SMILES string parsing fails for the string: {}".format(sm))

            # Generate conformer
            if not omega(mol_salt):
                raise ValueError("Error salt: Conformer generation fails for the "
                                 "molecule with smiles string: {}".format(sm))

            # Unique 3 code letter are set as solvent residue names
            solv_id = ''.join(random.sample(char_set * 3, 3))

            # Try to recognize the residue name
            oechem.OEPerceiveResidues(mol_salt)

            for atmol in mol_salt.GetAtoms():
                res = oechem.OEAtomGetResidue(atmol)
                if res.GetName() == 'UNL':
                    res.SetName(solv_id)
                    oechem.OEAtomSetResidue(atmol, res)
                    if solv_id not in solvent_resid_dic_names:
                        solvent_resid_dic_names[solv_id] = mol_salt
                else:
                    if res.GetName() not in solvent_resid_dic_names:
                        solvent_resid_dic_names[res.GetName()] = mol_salt
                    break

            oe_salt.append(mol_salt)

        n_salt = int(round(unit.AVOGADRO_CONSTANT_NA*salt_concentration*Volume.in_units_of(unit.liter)))

        # for i in range(0, len(salt_smiles)):
        #     print("Number of molecules for the salt component {} = {}".format(salt_smiles[i], n_salt))

        salt_sum_wgt = 0.0 * unit.grams/unit.mole
        for salt in oe_salt:
            # Molecular weight
            salt_sum_wgt += oechem.OECalculateMolecularWeight(salt) * unit.grams / unit.mole

        salt_sum_wgt_n_salt = salt_sum_wgt * n_salt

        # Create salt .pdb files
        if n_salt >= 1:
            salt_pdbs = []
            for i in range(0, len(salt_smiles)):
                pdb_name = os.path.basename(tempfile.mktemp(suffix='.pdb'))
                # pdb_name = salt_smiles[i] + '_' + pdb_name
                salt_pdbs.append(pdb_name)

            for i in range(0, len(salt_smiles)):
                ofs = oechem.oemolostream(salt_pdbs[i])
                oechem.OEWriteConstMolecule(ofs, oe_salt[i])

    # Container list of oemol solvent molecules generated by using smiles strings
    oe_solvents = []

    for sm in solvent_smiles:

        if sm == 'tip3p':
            tip3p_fn = os.path.join(PACKAGE_DIR, 'oeommtools', 'data', 'tip3p.pdb')
            ifs = oechem.oemolistream(tip3p_fn)
            mol_sol = oechem.OEMol()

            if not oechem.OEReadMolecule(ifs, mol_sol):
                raise IOError("It was not possible to read the tip3p molecule file")
        else:

            mol_sol = oechem.OEMol()

            if not oechem.OESmilesToMol(mol_sol, sm):
                raise ValueError("Error solvent: SMILES string parsing fails for the string: {}".format(sm))

            # Generate conformer
            if not omega(mol_sol):
                raise ValueError("Error solvent: Conformer generation fails for "
                                 "the molecule with smiles string: {}".format(sm))

        # Unique 3 code letter are set as solvent residue names
        solv_id = ''.join(random.sample(char_set*3, 3))

        # Try to recognize the residue name
        oechem.OEPerceiveResidues(mol_sol)

        for atmol in mol_sol.GetAtoms():
            res = oechem.OEAtomGetResidue(atmol)
            if res.GetName() == 'UNL':
                res.SetName(solv_id)
                oechem.OEAtomSetResidue(atmol, res)
                if solv_id not in solvent_resid_dic_names:
                    solvent_resid_dic_names[solv_id] = mol_sol
            else:
                if res.GetName() not in solvent_resid_dic_names:
                    solvent_resid_dic_names[res.GetName()] = mol_sol
                break

        oe_solvents.append(mol_sol)

    # Sum of the solvent molecular weights
    solvent_sum_wgt_frac = 0.0 * unit.grams/unit.mole

    for idx in range(0, len(oe_solvents)):
        # Molecular weight
        wgt = oechem.OECalculateMolecularWeight(oe_solvents[idx]) * unit.grams/unit.mole
        solvent_sum_wgt_frac += wgt * mol_fractions[idx]

    # Solute molecular weight
    solute_wgt = oechem.OECalculateMolecularWeight(solute_copy)*unit.gram/unit.mole

    # Estimate of the number of each molecular species present in the solution accordingly
    # to their molar fraction fi:
    #
    # ni = fi*(density*volume*NA - wgt_solute - sum_k(wgt_salt_k*nk) - wgt_ion*n_ion)/sum_j(wgt_nj * fj)
    #
    # where ni is the number of molecule of specie i, density the mixture density, volume the
    # mixture volume, wgt_solute the molecular weight of the solute, wgt_salt_k the molecular
    # weight of the salt component k, nk the number of molecule of salt component k, wgt_ion
    # the counter ion molecular weight, n_ions the number of counter ions and wgt_nj the molecular
    # weight of the molecule specie j with molar fraction fj

    div = (unit.AVOGADRO_CONSTANT_NA * density * Volume -
           (solute_wgt + salt_sum_wgt_n_salt + ion_sum_wgt_n_ions))/solvent_sum_wgt_frac

    # Solvent number of monomers
    n_monomers = [int(round(mf*div)) for mf in mol_fractions]

    if not all([nm > 0 for nm in n_monomers]):
        raise ValueError("Error negative number of solvent components: the density could be too low")

    # for i in range(0, len(solvent_smiles)):
    #     print("Number of molecules for the component {} = {}".format(solvent_smiles[i], n_monomers[i]))

    # Packmol Configuration file setting
    if close_solvent:
        header_template = """\n# Mixture\ntolerance {}\nfiletype pdb\noutput {}\nadd_amber_ter\navoid_overlap no"""
    else:
        header_template = """\n# Mixture\ntolerance {}\nfiletype pdb\noutput {}\nadd_amber_ter\navoid_overlap yes"""

    # Templates strings
    solute_template = """\n\n# Solute\nstructure {}\nnumber 1\nfixed 0. 0. 0. 0. 0. 0.\nresnumbers 1\nend structure"""

    if geometry == 'box':
        solvent_template = """\nstructure {}\nnumber {}\ninside box {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f}\
        \nchain {}\nresnumbers 3\nend structure"""
    if geometry == 'sphere':
        solvent_template = """\nstructure {}\nnumber {}\ninside sphere {:0.3f} {:0.3f} {:0.3f} {:0.3f}\
        \nchain {}\nresnumbers 3\nend structure"""

    # Create solvents .pdb files
    solvent_pdbs = []
    for i in range(0, len(solvent_smiles)):
        pdb_name = os.path.basename(tempfile.mktemp(suffix='.pdb'))
        solvent_pdbs.append(pdb_name)

    for i in range(0, len(solvent_smiles)):
        ofs = oechem.oemolostream(solvent_pdbs[i])
        oechem.OEWriteConstMolecule(ofs, oe_solvents[i])

    solute_pdb = 'solute' + '_' + os.path.basename(tempfile.mktemp(suffix='.pdb'))
    ofs = oechem.oemolostream(solute_pdb)

    if solute_copy.GetMaxConfIdx() > 1:
        raise ValueError("Solutes with multiple conformers are not supported")
    else:
        oechem.OEWriteConstMolecule(ofs, solute_copy)

    # Write Packmol header section
    mixture_pdb = 'mixture' + '_' + os.path.basename(tempfile.mktemp(suffix='.pdb'))
    body = header_template.format(distance_between_atoms, mixture_pdb)
    # Write Packmol configuration file solute section
    body += solute_template.format(solute_pdb)

    # The solute is centered inside the box
    xc = (BB_solute[0][0] + BB_solute[1][0]) / 2.
    yc = (BB_solute[0][1] + BB_solute[1][1]) / 2.
    zc = (BB_solute[0][2] + BB_solute[1][2]) / 2.

    # Correct for periodic box conditions to avoid
    # steric clashes at the box edges
    pbc_correction = 1.0 * unit.angstrom

    xmin = xc - ((box_edge - pbc_correction) / 2.) / unit.angstrom
    xmax = xc + ((box_edge - pbc_correction) / 2.) / unit.angstrom
    ymin = yc - ((box_edge - pbc_correction) / 2.) / unit.angstrom
    ymax = yc + ((box_edge - pbc_correction) / 2.) / unit.angstrom
    zmin = zc - ((box_edge - pbc_correction) / 2.) / unit.angstrom
    zmax = zc + ((box_edge - pbc_correction) / 2.) / unit.angstrom

    # Packmol setting for the solvent section
    body += '\n\n# Solvent'
    for i in range(0, len(solvent_smiles)):
        if geometry == 'box':
            body += solvent_template.format(solvent_pdbs[i],
                                            n_monomers[i],
                                            xmin, ymin, zmin,
                                            xmax, ymax, zmax,
                                            "!")
        if geometry == 'sphere':
            body += solvent_template.format(solvent_pdbs[i],
                                            n_monomers[i],
                                            xc, yc, zc,
                                            0.5 * box_edge / unit.angstrom,
                                            "!")

    # Packmol setting for the salt section
    if salt_concentration > 0.0 * unit.millimolar and n_salt >= 1:
        body += '\n\n# Salt'
        for i in range(0, len(salt_smiles)):
            if geometry == 'box':
                body += solvent_template.format(salt_pdbs[i],
                                                int(round(n_salt)),
                                                xmin, ymin, zmin,
                                                xmax, ymax, zmax,
                                                "?")
            if geometry == 'sphere':
                body += solvent_template.format(salt_pdbs[i],
                                                int(round(n_salt)),
                                                xc, yc, zc,
                                                0.5 * box_edge / unit.angstrom,
                                                "?")

    # Packmol setting for the ions section
    if neutralize_solute and n_ions >= 1:
        body += '\n\n# Counter Ions'
        for i in range(0, len(ions_smiles)):
            if geometry == 'box':
                body += solvent_template.format(ions_smiles_pdbs[i],
                                                n_ions,
                                                xmin, ymin, zmin,
                                                xmax, ymax, zmax,
                                                "*")
            if geometry == 'sphere':
                body += solvent_template.format(ions_smiles_pdbs[i],
                                                n_ions,
                                                xc, yc, zc,
                                                0.5 * box_edge / unit.angstrom,
                                                "*")

    # Packmol configuration file
    packmol_filename = os.path.basename(tempfile.mktemp(suffix='.inp'))

    with open(packmol_filename, 'w') as file_handle:
        file_handle.write(body)

    # Call Packmol
    if not verbose:
        mute_output = open(os.devnull, 'w')
        with open(packmol_filename, 'r') as file_handle:
            subprocess.check_call(['packmol'], stdin=file_handle, stdout=mute_output, stderr=mute_output)
    else:
        with open(packmol_filename, 'r') as file_handle:
            subprocess.check_call(['packmol'], stdin=file_handle)

    # Read in the Packmol solvated system
    solvated = oechem.OEMol()

    if os.path.exists(mixture_pdb+'_FORCED'):
        os.rename(mixture_pdb+'_FORCED', mixture_pdb)
        print("Warning: Packing solution is not optimal")

    with oechem.oemolistream(mixture_pdb) as ifs:
        oechem.OEReadMolecule(ifs, solvated)

    # To avoid to change the user oemol starting solute by reading in
    # the generated mixture pdb file and loosing molecule info, the
    # solvent molecules are extracted from the mixture system and
    # added back to the starting solute

    # Extract from the solution system the solvent molecule components
    hv_solvated = oechem.OEHierView(solvated, oechem.OEAssumption_BondedResidue +
                                    oechem.OEAssumption_ResPerceived)

    bv_solvent = oechem.OEBitVector(solvated.GetMaxAtomIdx())
    bv_salt = oechem.OEBitVector(solvated.GetMaxAtomIdx())
    bv_counter_ions = oechem.OEBitVector(solvated.GetMaxAtomIdx())

    for chain in hv_solvated.GetChains():
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():

                if chain.GetChainID() in ["!", "?", "*"]:

                    # Solvent Molecules
                    if chain.GetChainID() == "!":
                        bv = bv_solvent
                    # Salt Molecules
                    elif chain.GetChainID() == "?":
                        bv = bv_salt
                    # Counter Ions
                    elif chain.GetChainID() == "*":
                        bv = bv_counter_ions
                    else:
                        raise ValueError("Unknown Chain ID")

                    atms = hres.GetAtoms()

                    for at in atms:
                        bv.SetBitOn(at.GetIdx())

    solvent_comp = oechem.OEMol()
    salt_comp = oechem.OEMol()
    counter_ions_comp = oechem.OEMol()

    solvent_comp.SetTitle("Solvent")
    salt_comp.SetTitle("Salt")
    counter_ions_comp.SetTitle("Counter_ions")

    comp_list = list([solvent_comp, salt_comp, counter_ions_comp])
    bv_list = list([bv_solvent, bv_salt, bv_counter_ions])

    for bv, comp in zip(bv_list, comp_list):

        pred = oechem.OEAtomIdxSelected(bv)
        if not oechem.OESubsetMol(comp, solvated, pred):
            raise ValueError("Cannot extract the component from the solvated system")

        # Change ion NA, CL names to Na+, Cl-
        for at in comp.GetAtoms():
            res = oechem.OEAtomGetResidue(at)
            if res.GetName() == ' NA':
                res.SetName("Na+")
                oechem.OEAtomSetResidue(at, res)
            elif res.GetName() == ' CL':
                res.SetName("Cl-")
                oechem.OEAtomSetResidue(at, res)
            else:
                pass

    # Add the solvent components to the original solute copy
    solvated_system = solute.CreateCopy()

    # Set Title
    solvated_system.SetTitle(solute_copy.GetTitle())

    for comp in comp_list:
        oechem.OEAddMols(solvated_system, comp)

    if solvated_system.NumAtoms() != \
            solute.NumAtoms() + \
            solvent_comp.NumAtoms() + \
            salt_comp.NumAtoms() + \
            counter_ions_comp.NumAtoms():
        raise ValueError("Solvated system atom mismatch")

    # Cleaning
    to_delete = solvent_pdbs+[packmol_filename, solute_pdb, mixture_pdb]

    if salt_concentration > 0.0*unit.millimolar and n_salt >= 1:
        to_delete += salt_pdbs
    if neutralize_solute and n_ions >= 1:
        to_delete += ions_smiles_pdbs

    for fn in to_delete:
        try:
            os.remove(fn)
        except:
            pass

    # Calculate the solution total density
    total_wgt = oechem.OECalculateMolecularWeight(solvated_system)*unit.gram/unit.mole
    density_mix = (1/unit.AVOGADRO_CONSTANT_NA)*total_wgt/Volume
    print("Computed Solution Density = {}".format(density_mix.in_units_of(unit.gram/unit.milliliter)))
    # Threshold checking
    ths = 0.1 * unit.gram/unit.milliliter
    if not abs(density - density_mix.in_units_of(unit.gram/unit.milliliter)) < ths:
        raise ValueError("Error: the computed density for the solute {} does not match the selected density {} vs {}"
                         .format(solute_copy.GetTitle(), density_mix, density))

    if geometry == 'box':
        # Define the box vector and attached it as SD tag to the solvated system
        # with ID tag: 'box_vectors'
        box_vectors = (Vec3(box_edge/unit.angstrom, 0.0, 0.0),
                       Vec3(0.0, box_edge/unit.angstrom, 0.0),
                       Vec3(0.0, 0.0, box_edge/unit.angstrom))*unit.angstrom

        box_vectors = data_utils.encodePyObj(box_vectors)
        solvated_system.SetData(oechem.OEGetTag('box_vectors'), box_vectors)

    if return_components:
        components = list()
        for comp in comp_list:
            if comp.NumAtoms():
                components.append(comp)
            else:
                components.append(None)
        return tuple([solvated_system] + components)
    else:
        return solvated_system
