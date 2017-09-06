[![Build Status](https://travis-ci.org/oess/oeommtools.svg?branch=master)](https://travis-ci.org/oess/oeommtools)

[![Anaconda Badge](https://anaconda.org/nividic/oeommtools/badges/version.svg)](https://anaconda.org/openeye/oeommtools/badges/version.svg)

# OEOMMTools
These are a collection of tools developed to integrate and mix
the OE Toolkit with the OpenMM API


## Prerequisites
* OE Toolkit
* OpenMM

Authors
-------
* Gaetano Calabro' <gcalabro@eyesopen.com>

## Installation

conda install -c openeye/label/Orion oeommtools -c nividic -c omnia

Usage
-----
The following example solvate a benzene molecule in a box of water:

```python
from oeommtools import packmol
from openeye import oechem
from openeye import oeomega

smiles = 'c1ccccc1'
solute = oechem.OEMol()

omegaOpts = oeomega.OEOmegaOptions()
omegaOpts.SetMaxConfs(1)
omega = oeomega.OEOmega(omegaOpts)

if not oechem.OESmilesToMol(solute, smiles):
    oechem.OEThrow.Fatal("SMILES string parsing fails for the string: {}".format(smiles))
if not omega(solute):
    oechem.OEThrow.Fatal("Conformer generation fails for the molecule with smiles string: {}".format(smile))

solvate_mol = packmol.oesolvate_packmol(solute, density=1.0, padding_distance=10.0,
                                        solvents='[H]O[H]', molar_fractions='1.0', close_solvent=True,
                                        salt='[Na+], [Cl-]', salt_concentration=0.0, neutralize_solute=True)

ofs = oechem.oemolostream("solvated.pdb")
oechem.OEWriteConstMolecule(ofs, solvate_mol)
ofs.close()

```

## Issues
* OEOmmTools is in debugging stage and it has been tested on Ubuntu 16.04 and Mac OSX Sierra
* OeOmmTools has been developed in python 3.5

## Disclaimers
* This code is currently in alpha release status. Use at your own risk. We will almost certainly be making changes to the API in the near future.
