from setuptools import setup, find_packages

setup(
    name='pygly2-search',
    version='0.0.5',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
                "glycan-ms1 = pygly2.search.ms1app:taskmain",
                "glycan-ms2 = pygly2.search.ms2app:taskmain",
                "glycan-render = pygly2.search.ms2app:rerendermain",
                "glycan-serve = pygly2.search.flask_view:main"
        ],
    },
    namespace_packages=["pygly2", "pygly2.search", "pygly2.tests"]
)
