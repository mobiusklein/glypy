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
       ],
      annotate=True)
except ImportError, AttributeError:
    print("No Cython")
    extensions = [
      Extension('glypy.composition.ccomposition', sources=['glypy/composition/ccomposition.c']),
    ]

cmdclass = {}


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
    'plot': ["matplotlib>=1.4.3"],
    'glycomedb': ['lxml', 'requests'],
    'glyspace': ['requests', 'rdflib', "SPARQLWrapper"]
}

extras['all'] = list({d for extra in extras.values() for d in extra})


def run_setup(include_cext=True):
    setup(
          name='glypy',
          version='0.0.8',
          packages=find_packages(),
          include_package_data=True,
          package_data={
              "glypy.structure": ["glypy/structure/data/*"],
              "glypy.io.nomenclature": ["glypy/io/nomenclature/data/*"]
          },
          install_requires=required,
          extras_require=extras,
          # namespace_packages=[
          #   "glypy",
          #   "glypy.algorithms",
          #   "glypy.io",
          #   "glypy.io.nomenclature",
          #   "glypy.composition",
          #   "glypy.tests"
          # ],
          cmdclass=cmdclass,
          zip_safe=False,
          keywords="glycomics glycan carbohydrate glycoinformatics",
          description="A Glycoinformatics Toolkit",
          ext_modules=extensions if include_cext else None,
          url="https://github.com/mobiusklein/glypy",
          maintainer='Joshua Klein',
          maintainer_email="jaklein@bu.edu"
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
