from openeye import oechem
from openeye import oeomega
import numpy as np
import shutil
from simtk import unit
from simtk.openmm import Vec3
import subprocess
import tempfile
import os
from oeommtools import data_utils


def oesolvate(solute, density=1.0, padding_distance=10.0,
              solvents='[H]O[H]', molar_fractions='1.0', close_solvent=True,
              salt='[Na+], [Cl-]', salt_concentration=0.0, neutralize_solute=True, **kargs):
    """
    This function solvates the passed solute in a cubic box by using Packmol. Packmol creates
    an initial point for molecular dynamics simulations by packing molecule in defined regions
    of space. For additional info:
    http://www.ime.unicamp.br/~martinez/packmol/home.shtml

    The cubic box volume is estimated by the using the padding parameter and the solute size.
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
        The distance between the solute and the edge of the box in A
    solvents: python string
        A comma separated smiles string of the solvent molecules
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

    Return:
    -------
    oe_mol: OEMol
        The solvated system. A SD tag with name 'box_vector' is attached the
        output molecule containing the system box vectors
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

    # Remove smiles string with 0.0 mole fraction
    solvent_smiles = [solvents[i] for i, v in enumerate(fractions) if fractions[i]]
    mol_fractions = [mf for mf in fractions if mf]

    # Mole fractions are non-negative numbers
    if any([v < 0.0 for v in mol_fractions]):
        oechem.OEThrow.Fatal("Error: Mole fractions are non-negative real numbers")

    # Mole fractions must sum up to 1.0
    if abs(sum(mol_fractions) - 1.0) > 0.001:
        oechem.OEThrow.Fatal("Error: Mole fractions do not sum up to 1.0")

    # If the smiles string and mole fractions lists have different lengths raise an error
    if len(solvent_smiles) != len(mol_fractions):
        oechem.OEThrow.Fatal("Selected solvent number and selected molar fraction number mismatch: {} vs {}"
                             .format(len(solvent_smiles), len(mol_fractions)))

    # Set Units
    density = density * unit.grams/unit.milliliter
    padding_distance = padding_distance * unit.angstrom
    salt_concentration = salt_concentration * unit.millimolar

    # Calculate the Solute Bounding Box
    BB_solute = BoundingBox(solute)

    # Estimate of the box cube length
    box_edge = 2.0*padding_distance + np.max(BB_solute[1] - BB_solute[0])*unit.angstrom

    # Box Volume
    Volume = box_edge**3

    # Omega engine is used to generate conformations
    omegaOpts = oeomega.OEOmegaOptions()
    omegaOpts.SetMaxConfs(1)
    omega = oeomega.OEOmega(omegaOpts)

    # Neutralize solute
    ion_sum_wgt_n_ions = 0.0*unit.grams/unit.mole
    if neutralize_solute:
        # Container for the counter-ions
        oe_ions = []
        # Container for the ion smiles strings
        ions_smiles = []
        solute_formal_charge = 0
        for at in solute.GetAtoms():
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

        if n_ions >= 1:
            for sm in ions_smiles:
                mol = oechem.OEMol()
                if not oechem.OESmilesToMol(mol, sm):
                    oechem.OEThrow.Fatal("Error counter ions: SMILES string parsing fails for the string: {}".format(sm))
                oe_ions.append(mol)
                # Generate conformer
                if not omega(mol):
                    oechem.OEThrow.Fatal(
                        "Error counter ions: Conformer generation fails for the molecule with smiles string: {}".format(sm))

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

    # Add salt to the solution
    salt_sum_wgt_n_salt = 0.0 * unit.grams/unit.mole
    if salt_concentration > 0.0*unit.millimolar:

        salt_smiles = [sm.strip() for sm in salt.split(',')]

        # Container list of oemol salt molecules generated by using smiles strings
        oe_salt = []

        for sm in salt_smiles:
            mol = oechem.OEMol()
            if not oechem.OESmilesToMol(mol, sm):
                oechem.OEThrow.Fatal("Error salt: SMILES string parsing fails for the string: {}".format(sm))
            oe_salt.append(mol)
            # Generate conformer
            if not omega(mol):
                oechem.OEThrow.Fatal(
                    "Error salt: Conformer generation fails for the molecule with smiles string: {}".format(sm))

        n_salt = int(round(unit.AVOGADRO_CONSTANT_NA*salt_concentration*Volume.in_units_of(unit.liter)))

        # for i in range(0, len(salt_smiles)):
        #     print("Number of molecules for the salt component {} = {}".format(salt_smiles[i], n_salt))

        salt_sum_wgt = 0.0 * unit.grams/unit.mole
        for salt in oe_salt:
            # Molecular weight
            salt_sum_wgt += oechem.OECalculateMolecularWeight(salt) * unit.grams / unit.mole

        salt_sum_wgt_n_salt = salt_sum_wgt * n_salt

        # Create solvents .pdb files
        if n_salt >= 1:
            salt_smiles_pdbs = []
            for i in range(0, len(salt_smiles)):
                pdb_name = os.path.basename(tempfile.mktemp(suffix='.pdb'))
                pdb_name = salt_smiles[i] + '_' + pdb_name
                salt_smiles_pdbs.append(pdb_name)

            for i in range(0, len(salt_smiles)):
                ofs = oechem.oemolostream(salt_smiles_pdbs[i])
                oechem.OEWriteConstMolecule(ofs, oe_salt[i])

    # Container list of oemol solvent molecules generated by using smiles strings
    oe_solvents = []

    # Solvent smiles string parsing
    for sm in solvent_smiles:
        mol = oechem.OEMol()

        if not oechem.OESmilesToMol(mol, sm):
            oechem.OEThrow.Fatal("Error solvent: SMILES string parsing fails for the string: {}".format(sm))

        oe_solvents.append(mol)
        # Generate conformer
        if not omega(mol):
            oechem.OEThrow.Fatal("Error solvent: Conformer generation fails for "
                                 "the molecule with smiles string: {}".format(sm))

    # Sum of the solvent molecular weights
    solvent_sum_wgt_frac = 0.0 * unit.grams/unit.mole

    for idx in range(0, len(oe_solvents)):
        # Molecular weight
        wgt = oechem.OECalculateMolecularWeight(oe_solvents[idx]) * unit.grams/unit.mole
        solvent_sum_wgt_frac += wgt * mol_fractions[idx]

    # Solute molecular weight
    solute_wgt = oechem.OECalculateMolecularWeight(solute)*unit.gram/unit.mole

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
        oechem.OEThrow.Fatal("Error negative number of solvent components: the density could be too low")

    # for i in range(0, len(solvent_smiles)):
    #     print("Number of molecules for the component {} = {}".format(solvent_smiles[i], n_monomers[i]))

    # Packmol Configuration files setting
    if close_solvent:
        header_template = """\n# Mixture\ntolerance {}\nfiletype pdb\noutput {}\nadd_amber_ter\navoid_overlap no"""
    else:
        header_template = """\n# Mixture\ntolerance {}\nfiletype pdb\noutput {}\nadd_amber_ter\navoid_overlap yes"""

    # Templates strings
    solute_template = """\n\n# Solute\nstructure {}\nnumber 1\nfixed 0. 0. 0. 0. 0. 0.\nresnumbers 1\nend structure"""
    solvent_template = """\n\n# Solvent\nstructure {}\nnumber {}\ninside box {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f}\nchain !\nend structure"""
    salt_template = """\n\n# Salt\nstructure {}\nnumber {}\ninside box {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f}\nchain !\nend structure"""
    ions_template = """\n\n# Counter-Ions\nstructure {}\nnumber {}\ninside box {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f} {:0.3f}\nchain !\nend structure"""

    # Create solvents .pdb files
    solvent_smiles_pdbs = []
    for i in range(0, len(solvent_smiles)):
        pdb_name = os.path.basename(tempfile.mktemp(suffix='.pdb'))
        pdb_name = solvent_smiles[i] + '_' + pdb_name
        solvent_smiles_pdbs.append(pdb_name)

    for i in range(0, len(solvent_smiles)):
        ofs = oechem.oemolostream(solvent_smiles_pdbs[i])
        oechem.OEWriteConstMolecule(ofs, oe_solvents[i])

    solute_pdb = 'solute'+'_'+os.path.basename(tempfile.mktemp(suffix='.pdb'))
    ofs = oechem.oemolostream(solute_pdb)
    oechem.OEWriteConstMolecule(ofs, solute)

    # Write Packmol header section
    mixture_pdb = 'mixture'+'_'+os.path.basename(tempfile.mktemp(suffix='.pdb'))
    body = header_template.format(2.0, mixture_pdb)
    # Write Packmol configuration file solute section
    body += solute_template.format(solute_pdb)

    # The solute is centered inside the box
    xc = (BB_solute[0][0]+BB_solute[1][0])/2.
    yc = (BB_solute[0][1]+BB_solute[1][1])/2.
    zc = (BB_solute[0][2]+BB_solute[1][2])/2.
    xmin = xc - (box_edge/2.)/unit.angstrom
    xmax = xc + (box_edge/2.)/unit.angstrom
    ymin = yc - (box_edge/2.)/unit.angstrom
    ymax = yc + (box_edge/2.)/unit.angstrom
    zmin = zc - (box_edge/2.)/unit.angstrom
    zmax = zc + (box_edge/2.)/unit.angstrom

    # Packmol setting for the solvent section
    for i in range(0, len(solvent_smiles)):
        body += solvent_template.format(solvent_smiles_pdbs[i],
                                        n_monomers[i],
                                        xmin, ymin, zmin,
                                        xmax, ymax, zmax)
    # Packmol setting for the salt section
    if salt_concentration > 0.0*unit.millimolar and n_salt >= 1:
        for i in range(0, len(salt_smiles)):
            body += salt_template.format(salt_smiles_pdbs[i],
                                         int(round(n_salt)),
                                         xmin, ymin, zmin,
                                         xmax, ymax, zmax)

    # Packmol setting for the ions section
    if neutralize_solute and n_ions >= 1:
        for i in range(0, len(ions_smiles)):
            body += ions_template.format(ions_smiles_pdbs[i],
                                         n_ions,
                                         xmin, ymin, zmin,
                                         xmax, ymax, zmax)

    # Packmol configuration file
    packmol_filename = os.path.basename(tempfile.mktemp(suffix='.inp'))

    with open(packmol_filename, 'w') as file_handle:
        file_handle.write(body)

    # Call Packmol
    mute_output = open(os.devnull, 'w')
    with open(packmol_filename, 'r') as file_handle:
        subprocess.check_call(['packmol'], stdin=file_handle, stdout=mute_output, stderr=mute_output)

    # Read in the Packmol solvated system
    solvated = oechem.OEMol()
    ifs = oechem.oemolistream(mixture_pdb)
    oechem.OEReadMolecule(ifs, solvated)

    # To avoid to change the user oemol starting solute by reading in
    # the generated mixture pdb file and loosing molecule info, the
    # solvent molecules are extracted from the mixture system and
    # added back to the starting solute

    # Create a string code to identify the solute residues. The code ID used is based
    # on the residue number id, the residue name and the chain id:
    # id+resname+chainID
    hv_solute = oechem.OEHierView(solute, oechem.OEAssumption_BondedResidue +
                                  oechem.OEAssumption_ResPerceived)
    solute_resid_list = []
    for chain in hv_solute.GetChains():
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                oe_res = hres.GetOEResidue()
                solute_resid_list.append(str(oe_res.GetResidueNumber())+oe_res.GetName()+chain.GetChainID())

    # Extract from the solution system the solvent molecules
    # by checking the previous solute generated ID: id+resname+chainID
    hv_solvated = oechem.OEHierView(solvated, oechem.OEAssumption_BondedResidue +
                                    oechem.OEAssumption_ResPerceived)
    bv = oechem.OEBitVector(solvated.GetMaxAtomIdx())
    for chain in hv_solvated.GetChains():
        for frag in chain.GetFragments():
            for hres in frag.GetResidues():
                oe_res = hres.GetOEResidue()
                if str(oe_res.GetResidueNumber())+oe_res.GetName()+chain.GetChainID() not in solute_resid_list:
                    atms = hres.GetAtoms()
                    for at in atms:
                        bv.SetBitOn(at.GetIdx())

    pred = oechem.OEAtomIdxSelected(bv)
    components = oechem.OEMol()
    oechem.OESubsetMol(components, solvated, pred)

    # Add the solvent molecules to the solute copy
    solvated_system = solute.CreateCopy()
    oechem.OEAddMols(solvated_system, components)

    # Cleaning
    to_delete = solvent_smiles_pdbs+[packmol_filename, solute_pdb, mixture_pdb]

    if salt_concentration > 0.0*unit.millimolar and n_salt >= 1:
        to_delete += salt_smiles_pdbs
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
        raise oechem.OEThrow.Fatal("Error: the computed density does not match "
                                   "the selected density {} vs {}"
                                   .format(density_mix, density))

    # Define the box vector and attached it as SD tag to the solvated system
    # with ID tag: 'box_vectors'
    box_vectors = (Vec3(box_edge/unit.angstrom, 0.0, 0.0),
                   Vec3(0.0, box_edge/unit.angstrom, 0.0),
                   Vec3(0.0, 0.0, box_edge/unit.angstrom))*unit.angstrom

    box_vectors = data_utils.PackageOEMol.encodePyObj(box_vectors)
    solvated_system.SetData(oechem.OEGetTag('box_vectors'), box_vectors)

    return solvated_system