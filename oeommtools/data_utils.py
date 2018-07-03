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