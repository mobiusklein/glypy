import os
import sys

command_string = (r"nosetests --with-coverage --with-timer --cover-package=glypy "
                  r"--cover-html --cover-html-dir=test_reports --logging-level=DEBUG -v --with-id %s")


def run(*tests):
    if len(tests) == 0:
        tests = [default_tests]
    cmd = (command_string % ' '.join(tests))
    print cmd
    os.system(cmd)

default_tests = "glypy/tests/"

if __name__ == '__main__':
    run(*sys.argv[1:])
