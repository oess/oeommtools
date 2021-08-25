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

import base64
import parmed
from openeye import oechem
import pickle


def getData(molecule, tag):
    return molecule.GetData(oechem.OEGetTag(str(tag)))


def encodePyObj(py_obj):
    """Encode Python Objects/ParmEd Structure to BYTES (base64)"""
    pkl_obj = pickle.dumps(py_obj)
    return base64.b64encode(pkl_obj)


def decodePyObj(data):
    """Decode the Base64 encoded Python Object/ParmEd Structure"""
    decoded_obj = base64.b64decode(data)
    return pickle.loads(decoded_obj)