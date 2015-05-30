from setuptools import setup, find_packages

setup(
    name='glypy-search',
    version='0.0.5',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
                "glycan-ms1 = glypy.search.ms1app:taskmain",
                "glycan-ms2 = glypy.search.ms2app:taskmain",
                "glycan-render = glypy.search.ms2app:rerendermain",
                "glycan-serve = glypy.search.flask_view:main"
        ],
    },
    namespace_packages=["glypy", "glypy.search", "glypy.tests"]
)
