import sys
from setuptools import setup, find_packages, Extension

# With gratitude to the SqlAlchemy setup.py authors

from distutils.command.build_ext import build_ext
from distutils.errors import (CCompilerError, DistutilsExecError,
                              DistutilsPlatformError)

ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)
if sys.platform == 'win32':
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    ext_errors += (IOError,)

try:
    from Cython.Build import cythonize
    extensions = cythonize(
        [Extension("glypy.composition.ccomposition", ["glypy/composition/ccomposition.pyx"]),
         Extension("glypy.utils.cenum", ['glypy/utils/cenum.pyx']),
         Extension("glypy._c.utils", ["glypy/_c/utils.pyx"]),
         Extension("glypy._c.structure.glycan_composition", [
                   "glypy/_c/structure/glycan_composition.pyx"]),
         ],
        annotate=True)
except (ImportError, AttributeError):
    print("No Cython")
    extensions = [
        Extension('glypy.composition.ccomposition', sources=[
                  'glypy/composition/ccomposition.c']),
        Extension("glypy.utils.cenum", ['glypy/utils/cenum.c']),
        Extension("glypy._c.utils", ["glypy/_c/utils.c"]),
        Extension("glypy._c.structure.glycan_composition", [
                  "glypy/_c/structure/glycan_composition.c"]),
    ]

cmdclass = {}


with open("glypy/version.py") as version_file:
    version = None
    for line in version_file.readlines():
        if "version = " in line:
            version = line.split(" = ")[1].replace("\"", "").strip()
            break
    else:
        print("Cannot determine version")

long_description = ''
try:
    with open("README.rst") as readme_file:
        long_description = readme_file.read()
except Exception:
    pass


class BuildFailed(Exception):

    def __init__(self):
        self.cause = sys.exc_info()[1]  # work around py 2/3 different syntax

    def __str__(self):
        return str(self.cause)


class ve_build_ext(build_ext):
    # This class allows C extension building to fail.

    def run(self):
        try:
            build_ext.run(self)
        except DistutilsPlatformError:
            raise BuildFailed()

    def build_extension(self, ext):
        try:
            build_ext.build_extension(self, ext)
        except ext_errors:
            raise BuildFailed()
        except ValueError:
            # this can happen on Windows 64 bit, see Python issue 7511
            if "'path'" in str(sys.exc_info()[1]):  # works with both py 2/3
                raise BuildFailed()
            raise


cmdclass['build_ext'] = ve_build_ext


def status_msgs(*msgs):
    print('*' * 75)
    for msg in msgs:
        print(msg)
    print('*' * 75)


required = ["hjson", "six"]


extras = {
    'plot': ["matplotlib>=2.2.0"],
    'glyspace': ['requests', 'rdflib', "SPARQLWrapper"]
}

extras['all'] = list({d for extra in extras.values() for d in extra})


def run_setup(include_cext=True):
    setup(
        name='glypy',
        version=version,
        packages=find_packages(),
        include_package_data=True,
        package_data={
            "glypy.structure": ["glypy/structure/data/*"],
            "glypy.io": ["glypy/io/data/*"],
            "glypy.io.nomenclature": ["glypy/io/nomenclature/data/*"],
        },
        install_requires=required,
        extras_require=extras,
        cmdclass=cmdclass,
        zip_safe=False,
        keywords="glycomics glycan carbohydrate glycoinformatics glypy n-linked o-linked glycosaminoglycan",
        description="A Glycoinformatics Toolkit",
        long_description=long_description,
        ext_modules=extensions if include_cext else None,
        url="https://github.com/mobiusklein/glypy",
        maintainer='Joshua Klein',
        maintainer_email="jaklein@bu.edu",
        classifiers=[
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: Apache Software License',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3'],
        project_urls={
            "Documentation": "https://glypy.readthedocs.io/",
            "Repository": "https://github.com/mobiusklein/glypy",
        }
    )


try:
    run_setup(True)
except Exception as exc:
    status_msgs(
        str(exc),
        "WARNING: The C extension could not be compiled, " +
        "speedups are not enabled.",
        "Failure information, if any, is above.",
        "Retrying the build without the C extension now."
    )

    run_setup(False)

    status_msgs(
        "WARNING: The C extension could not be compiled, " +
        "speedups are not enabled.",
        "Plain-Python build succeeded."
    )
