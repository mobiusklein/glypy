from setuptools import setup, find_packages

setup(
    name='pygly2-search',
    version='0.0.5',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
                "pygly-ms2 = pygly2.search.app:taskmain",
        ],
    },
    namespace_packages=["pygly2", "pygly2.search", "pygly2.tests"]
)
