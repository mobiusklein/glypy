from setuptools import setup, find_packages

setup(
    name='pygly2-search',
    version='0.0.5',
    packages=find_packages(),
    namespace_packages=["pygly2", "pygly2.search", "pygly2.tests"]
)
