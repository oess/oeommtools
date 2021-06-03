import unittest
from simtk import unit
from simtk.openmm import app
from openeye import oechem
from oeommtools import utils
from oeommtools import data_utils
from oeommtools import packmol
import pickle
import numpy as np


class ConversionTester(unittest.TestCase):
    """
    Test conversion functions between OE and OpenMM topologies
    """
    def setUp(self):
        pass

    def test_oemol_to_openmmTop(self):
        protein_fn = "tests/data/T4-protein.pdb"
        mol = oechem.OEMol()
        ifs = oechem.oemolistream(protein_fn)
        oechem.OEReadMolecule(ifs, mol)

        top, omm_pos = utils.oemol_to_openmmTop(mol)

        # Assert Atom numbers
        self.assertEqual(top.getNumAtoms(), mol.NumAtoms())

        for (op_at, oe_at) in zip(top.atoms(), mol.GetAtoms()):
            # Assert atom indexes
            self.assertEqual(op_at.index, oe_at.GetIdx())

        oe_pos = [v for k, v in mol.GetCoords().items()]
        # Assert atom positions
        self.assertEqual(oe_pos, omm_pos.in_units_of(unit.angstrom)/unit.angstrom)

        # Assert bond order
        dic_bond_openmm = dict()
        for bond in top.bonds():
            # OpenMM atoms
            at0_idx = bond[0].index
            at1_idx = bond[1].index
            if at0_idx < at1_idx:
                dic_bond_openmm[(at0_idx, at1_idx)] = (bond.order, bond.type)
            else:
                dic_bond_openmm[(at1_idx, at0_idx)] = (bond.order, bond.type)

        dic_bond_oe = dict()
        for bond in mol.GetBonds():
            # OE atoms
            at0_idx = bond.GetBgnIdx()
            at1_idx = bond.GetEndIdx()
            if at0_idx < at1_idx:
                dic_bond_oe[(at0_idx, at1_idx)] = (bond.GetOrder(), bond.GetType())
            else:
                dic_bond_oe[(at1_idx, at0_idx)] = (bond.GetOrder(), bond.GetType())

        self.assertEqual(dic_bond_openmm, dic_bond_oe)

    def test_openmmTop_to_oemol(self):
        protein_fn = 'tests/data/T4-protein.pdb'

        pdb = app.PDBFile(protein_fn)

        oe_mol = utils.openmmTop_to_oemol(pdb.topology, pdb.positions)

        # Assert
        self.assertEqual(pdb.topology.getNumAtoms(), oe_mol.NumAtoms())

        for (op_at, oe_at) in zip(pdb.topology.atoms(), oe_mol.GetAtoms()):
            self.assertEqual(op_at.index, oe_at.GetIdx())

        oe_pos = [v for k, v in oe_mol.GetCoords().items()]
        np.testing.assert_almost_equal(pdb.getPositions(asNumpy=True).in_units_of(unit.angstrom)/unit.angstrom,
                                       np.array(oe_pos), decimal=2)


class SelectionLanguageTester(unittest.TestCase):
    """
    Test Selection Language for a Complex system
    """

    def setUp(self):
        pass

    # Check Selection Language
    def test_selection_language(self):

        fcomplex = "tests/data/pP38_lp38a_2x_complex.pdb"
        # Read OEMol molecule
        mol = oechem.OEMol()
        with oechem.oemolistream(fcomplex) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        res_dic = {}

        mask = 'protein'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'ligand'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'water'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'ions'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'excipients'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'ca_protein'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'protein or ligand'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'noh (protein or ligand)'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'ca_protein or (noh ligand)'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = 'resid A:4 A:5 A:6 A:7'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        mask = '5.0 around ligand'
        ind_set = utils.select_oemol_atom_idx_by_language(mol, mask=mask)
        res_dic[mask] = ind_set

        dic_fname = "tests/data/restraint_test_P38_lp38a_2x.pickle"

        file = open(dic_fname, 'rb')
        res_dic_loaded = pickle.load(file)

        for k in res_dic_loaded:
            if res_dic[k] == res_dic_loaded[k]:
                pass
            else:
                raise ValueError("Restraints checking Errors on mask: {}".format(k))


class SolvatePackmolTester(unittest.TestCase):
    """
    Test Solvation Packmol
    """

    def setUp(self):
        pass

    def test_solvation_packmol(self):
        # Complex file name
        fcomplex = "tests/data/Bace_protein.pdb"

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(fcomplex) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Solvate the system
        solv_complex = packmol.oesolvate(mol, density=1.0,
                                         padding_distance=10.0,
                                         distance_between_atoms=2.0,
                                         solvents='[H]O[H]',
                                         molar_fractions='1.0',
                                         geometry='box',
                                         close_solvent=True,
                                         salt='[Na+], [Cl-]', salt_concentration=100.0,
                                         neutralize_solute=True,
                                         return_components=False)

        prot, lig, wat, cofactors, lipids, metals, excipients = utils.split(solv_complex)

        npa = prot.GetMaxAtomIdx()
        nla = lig.GetMaxAtomIdx()
        noa = excipients.GetMaxAtomIdx()
        nwa = wat.GetMaxAtomIdx()
        nlpa = lipids.GetMaxAtomIdx()
        nma = metals.GetMaxAtomIdx()
        nca = cofactors.GetMaxAtomIdx()

        self.assertEqual(npa, 6044)
        self.assertEqual(nla, 0)
        self.assertEqual(nlpa, 0)
        self.assertEqual(nca, 0)
        self.assertEqual(nma, 0)
        # Ions added to excipients
        self.assertEqual(noa, 94)
        # Water molecules added
        self.assertAlmostEqual(nwa, 48279, delta=90)

        box_vectors = solv_complex.GetData('box_vectors')
        box_vectors = data_utils.decodePyObj(box_vectors)
        box_vectors = box_vectors.in_units_of(unit.nanometers)

        self.assertAlmostEqual(box_vectors[0][0] / unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[0][1] / unit.nanometers, 0.0)
        self.assertEqual(box_vectors[0][2] / unit.nanometers, 0.0)

        self.assertAlmostEqual(box_vectors[1][1] / unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[1][0] / unit.nanometers, 0.0)
        self.assertEqual(box_vectors[1][2] / unit.nanometers, 0.0)

        self.assertAlmostEqual(box_vectors[2][2] / unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[2][0] / unit.nanometers, 0.0)
        self.assertEqual(box_vectors[2][1] / unit.nanometers, 0.0)

    def test_solvation_packmol_tip3p(self):
        # Complex file name
        fcomplex = "tests/data/Bace_protein.pdb"

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(fcomplex) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Solvate the system
        solv_complex = packmol.oesolvate(mol, density=1.0,
                                         padding_distance=10.0,
                                         distance_between_atoms=2.0,
                                         solvents='tip3p',
                                         molar_fractions='1.0',
                                         geometry='box',
                                         close_solvent=True,
                                         salt='[Na+], [Cl-]', salt_concentration=100.0,
                                         neutralize_solute=True, return_components=False)

        prot, lig, wat,  cofactors, lipids, metals, excipients = utils.split(solv_complex)

        npa = prot.GetMaxAtomIdx()
        nla = lig.GetMaxAtomIdx()
        noa = excipients.GetMaxAtomIdx()
        nwa = wat.GetMaxAtomIdx()
        nlpa = lipids.GetMaxAtomIdx()
        nma = metals.GetMaxAtomIdx()
        nca = cofactors.GetMaxAtomIdx()

        self.assertEqual(npa, 6044)
        self.assertEqual(nla, 0)
        self.assertEqual(nlpa, 0)
        self.assertEqual(nca, 0)
        self.assertEqual(nma, 0)

        # Ions added to excipients
        self.assertEqual(noa, 94)
        # Water molecules added
        self.assertAlmostEqual(nwa, 48279, delta=90)

        box_vectors = solv_complex.GetData('box_vectors')
        box_vectors = data_utils.decodePyObj(box_vectors)
        box_vectors = box_vectors.in_units_of(unit.nanometers)

        self.assertAlmostEqual(box_vectors[0][0] / unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[0][1] / unit.nanometers, 0.0)
        self.assertEqual(box_vectors[0][2] / unit.nanometers, 0.0)

        self.assertAlmostEqual(box_vectors[1][1] / unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[1][0] / unit.nanometers, 0.0)
        self.assertEqual(box_vectors[1][2] / unit.nanometers, 0.0)

        self.assertAlmostEqual(box_vectors[2][2] / unit.nanometers, 8.23, delta=0.01)
        self.assertEqual(box_vectors[2][0] / unit.nanometers, 0.0)
        self.assertEqual(box_vectors[2][1] / unit.nanometers, 0.0)

    def test_solvation_packmol_components(self):
        # Complex file name
        fcomplex = "tests/data/Bace_protein.pdb"

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(fcomplex) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Solvate the system
        solv_complex, solvent, salt, counter_ions = packmol.oesolvate(mol, density=1.0,
                                                                      padding_distance=10.0,
                                                                      distance_between_atoms=2.0,
                                                                      solvents='tip3p',
                                                                      molar_fractions='1.0',
                                                                      geometry='box',
                                                                      close_solvent=True,
                                                                      salt='[Na+], [Cl-]', salt_concentration=100.0,
                                                                      neutralize_solute=True,
                                                                      return_components=True)

        prot, lig, wat, cofactors, lipids, metals, excipients = utils.split(solv_complex)

        npa = prot.GetMaxAtomIdx()
        nla = lig.GetMaxAtomIdx()
        noa = excipients.GetMaxAtomIdx()
        nwa = wat.GetMaxAtomIdx()
        nlpa = lipids.GetMaxAtomIdx()
        nma = metals.GetMaxAtomIdx()
        nca = cofactors.GetMaxAtomIdx()

        self.assertEqual(npa, 6044)
        self.assertEqual(nla, 0)
        self.assertEqual(nlpa, 0)
        self.assertEqual(nca, 0)
        self.assertEqual(nma, 0)

        # The Bace system has 19 water molecules
        self.assertEqual(nwa, solvent.NumAtoms() + 19 * 3)

        # The Bace system has the TLA excipient with 14 atoms
        self.assertEqual(noa, salt.NumAtoms() + counter_ions.NumAtoms() + 14)

    def test_solvation_packmol_components_None(self):
        # Complex file name
        fcomplex = "tests/data/Bace_protein.pdb"

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(fcomplex) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        # Solvate the system
        solv_complex, solvent, salt, counter_ions = packmol.oesolvate(mol, density=1.0,
                                                                      padding_distance=10.0,
                                                                      distance_between_atoms=2.0,
                                                                      solvents='tip3p',
                                                                      molar_fractions='1.0',
                                                                      geometry='box',
                                                                      close_solvent=True,
                                                                      salt='[Na+], [Cl-]', salt_concentration=0.0,
                                                                      neutralize_solute=False,
                                                                      return_components=True)

        prot, lig, wat, cofactors, lipids, metals, excipients = utils.split(solv_complex)

        npa = prot.GetMaxAtomIdx()
        nla = lig.GetMaxAtomIdx()

        self.assertEqual(npa, 6044)
        self.assertEqual(nla, 0)

        assert salt is None
        assert counter_ions is None


class RemoveWaterIonsTester(unittest.TestCase):
    """
    Test Remove Water and Ions
    """

    def setUp(self):
        pass

    def test_remove_water_ions(self):
        # Complex file name
        fcomplex = "tests/data/pP38_lp38a_2x_complex.pdb"

        # Read OEMol molecule
        mol = oechem.OEMol()

        with oechem.oemolistream(fcomplex) as ifs:
            oechem.OEReadMolecule(ifs, mol)

        clean_system = utils.strip_water_ions(mol)

        prot, lig, wat, cofactors, lipids, metals, excipients = utils.split(clean_system)

        npa = prot.GetMaxAtomIdx()
        nla = lig.GetMaxAtomIdx()
        noa = excipients.GetMaxAtomIdx()
        nwa = wat.GetMaxAtomIdx()
        nlpa = lipids.GetMaxAtomIdx()
        nma = metals.GetMaxAtomIdx()
        nca = cofactors.GetMaxAtomIdx()

        self.assertEqual(npa, 5629)
        self.assertEqual(nla, 43)
        self.assertEqual(noa, 0)
        self.assertEqual(nwa, 0)
        self.assertEqual(nlpa, 0)
        self.assertEqual(nca, 0)
        self.assertEqual(nma, 0)


class SplitterTester(unittest.TestCase):
    """
    Test Splitter
    """

    def setUp(self):
        pass

    def test_splitter(self):
        # Complex file name
        flask_name = "tests/data/thrombin_flask.oeb"

        # Read OEMol molecule
        flask = oechem.OEMol()

        with oechem.oemolistream(flask_name) as ifs:
            oechem.OEReadMolecule(ifs, flask)

        prot, lig, wat, cofactors, lipids, metals, excipients = utils.split(flask, ligand_res_name='LIG')

        self.assertEqual(prot.NumAtoms(), 4686)
        self.assertEqual(lig.NumAtoms(), 52)
        self.assertEqual(wat.NumAtoms(), 33783)
        self.assertEqual(cofactors.NumAtoms(), 0)
        self.assertEqual(lipids.NumAtoms(), 0)
        self.assertEqual(metals.NumAtoms(), 0)
        self.assertEqual(excipients.NumAtoms(), 25)


if __name__ == "__main__":
        unittest.main()
