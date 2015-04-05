import sys
from setuptools import setup, find_packages, Extension
from Cython.Build import cythonize

# With gratitude to the SqlAlchemy setup.py authors

from distutils.command.build_ext import build_ext
from distutils.errors import (CCompilerError, DistutilsExecError,
                              DistutilsPlatformError)

ext_errors = (CCompilerError, DistutilsExecError, DistutilsPlatformError)
if sys.platform == 'win32':
    # 2.6's distutils.msvc9compiler can raise an IOError when failing to
    # find the compiler
    ext_errors += (IOError,)

extensions = cythonize([Extension("pygly2.composition.ccomposition", ["pygly2/composition/ccomposition.pyx"])],
                       annotate=True,
                       profile=True)

cmdclass = {}


class BuildFailed(Exception):

    def __init__(self):
        self.cause = sys.exc_info()[1]  # work around py 2/3 different syntax


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


def run_setup(include_cext=True):
    setup(
          name='pygly2',
          version='0.0.5',
          packages=find_packages(),
          include_package_data=True,
          package_data={
              "pygly2.structure": ["pygly2/structure/data/*"],
              "pygly2.io.nomenclature": ["pygly2/io/nomenclature/data/*"]
          },
          namespace_packages=["pygly2", "pygly2.io"],
          cmdclass=cmdclass,
          zip_safe=False,
          ext_modules=extensions if include_cext else None,
    )

try:
    run_setup(True)
except BuildFailed as exc:
    status_msgs(
        exc.cause,
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
