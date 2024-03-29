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

import os

import shutil

from invoke import task, run

import importlib

import oeommtools


PACKAGE_DIR = os.path.dirname(os.path.dirname(oeommtools.__file__))


@task
def setversion(ctx, new_version):
    """
    Set the package version
    """

    clean(ctx)

    fn = os.path.join(PACKAGE_DIR, "oeommtools", "__init__.py")

    with open(fn, "r") as f:
        lines = f.readlines()

    lines = [
        "__version__ = '{}'\n".format(new_version) if "__version__" in line else line
        for line in lines
    ]

    with open(fn, "w") as f:
        f.writelines(lines)

    fn = os.path.join(PACKAGE_DIR, "devtools/conda-recipe", "meta.yaml")

    with open(fn, "r") as f:
        lines = f.readlines()

    lines = [
        "  version: '{}'\n".format(new_version) if "version" in line else line
        for line in lines
    ]

    with open(fn, "w") as f:
        f.writelines(lines)

    importlib.reload(oeommtools)


@task
def version(ctx):
    print(oeommtools.__version__)


@task
def anaconda_upload(ctx, username, label="OrionDev"):
    conda_build_path = os.environ['CONDA_PREFIX'].split('envs')[0] + "conda-bld/noarch/"
    fn = oeommtools.__name__ + "-" + oeommtools.__version__ + "-py_0.tar.bz2"

    full_path = os.path.join(conda_build_path, fn)

    if not os.path.isfile(full_path):
        try:
            print("\n>>>>>>>>>>> Building the conda recipe <<<<<<<<<<<<<\n")
            pkg_meta_file = os.path.join(os.path.dirname(oeommtools.__path__[0]), "devtools/conda-recipe/")
            run("conda build {} -c conda-forge -c omnia -c OpenEye/label/Orion".format(pkg_meta_file))
            print("\n>>>>>>>>>>> End building the conda recipe <<<<<<<<<<<<<\n")
        except Exception as e:
            if not os.path.isfile(full_path):
                raise ValueError("Something went wrong during the conda Building {}".format(str(e)))
    try:
        print("\nAnaconda Uploading....\n")
        run("anaconda upload --user {} -l {} {}".format(username, label, full_path))
        print("\nDone\n")
    except Exception as e:
        raise IOError("It was not possible to upload the file: {}".format(str(e)))


@task
def clean(ctx):
    """
    Clean up doc and package builds
    """
    clean_pyc(ctx)
    clean_docs(ctx)
    clean_pycache(ctx)
    shutil.rmtree("dist", ignore_errors=True)
    shutil.rmtree("build", ignore_errors=True)
    egg_path = "{}.egg-info".format("oeommtools".replace("-", "_"))
    if os.path.isfile(egg_path):
        os.remove(egg_path)
    elif os.path.isdir(egg_path):
        shutil.rmtree(egg_path)
    shutil.rmtree(".pytest_cache", ignore_errors=True)


@task
def clean_pyc(ctx):
    """
    cleans out .pyc files
    """
    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith(".pyc"):
                filename = os.path.join(root, file)
                if os.path.exists(filename):
                    os.unlink(filename)


@task
def clean_pycache(ctx):
    """
    cleans out __pycache__ dirs
    """
    for dirpath, dirs, files in os.walk(os.getcwd()):
        for dir in dirs:
            if dir == "__pycache__":
                del_dir = os.path.join(dirpath, dir)
                shutil.rmtree(del_dir)


@task
def clean_docs(ctx):
    doc_dir = "docs/build/html"
    _clean_out_dir(doc_dir)

    if os.path.isdir("docs/build/doctrees"):
        shutil.rmtree("docs/build/doctrees")


def _clean_out_dir(dir_path):
    if os.path.isdir(dir_path):
        for the_file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, the_file)
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(e)
