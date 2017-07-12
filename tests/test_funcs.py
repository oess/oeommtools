import unittest
from simtk import unit
from simtk.openmm import app
from openeye import oechem
from oeommtools import utils


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
        dic_bond_openmm = {}
        for bond in top.bonds():
            # OpenMM atoms
            at0_idx = bond[0].index
            at1_idx = bond[1].index
            if at0_idx < at1_idx:
                dic_bond_openmm[(at0_idx, at1_idx)] = bond.order
            else:
                dic_bond_openmm[(at1_idx, at0_idx)] = bond.order

        dic_bond_oe = {}
        for bond in mol.GetBonds():
            # OE atoms
            at0_idx = bond.GetBgnIdx()
            at1_idx = bond.GetEndIdx()
            if at0_idx < at1_idx: 
                dic_bond_oe[(at0_idx, at1_idx)] = bond.GetOrder()
            else:
                dic_bond_oe[(at1_idx, at0_idx)] = bond.GetOrder()

        self.assertEqual(dic_bond_openmm, dic_bond_oe)

    def test_openmmTop_to_oemol(self):
        protein_fn ='tests/data/T4-protein.pdb'
        
        pdb = app.PDBFile(protein_fn)

        oe_mol = utils.openmmTop_to_oemol(pdb.topology, pdb.positions)

        # Assert 
        self.assertEqual(pdb.topology.getNumAtoms(), oe_mol.NumAtoms())

        for (op_at, oe_at) in zip(pdb.topology.atoms(), oe_mol.GetAtoms()):
            self.assertEqual(op_at.index, oe_at.GetIdx())

        import numpy as np
        oe_pos = [v for k, v in oe_mol.GetCoords().items()]
        np.testing.assert_almost_equal(pdb.getPositions(asNumpy=True).in_units_of(unit.angstrom)/unit.angstrom,
                                       np.array(oe_pos), decimal=2)
     
        
if __name__ == "__main__":
        unittest.main()
