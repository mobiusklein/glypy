from setuptools import setup, find_packages

setup(
    name='pygly2',
    version='0.0.5',
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "pygly2.structure": ["pygly2/structure/data/*"],
        "pygly2.io.nomenclature": ["pygly2/io/nomenclature/data/*"]
    },
    namespace_packages=["pygly2", "pygly2.io"]
)