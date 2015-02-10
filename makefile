all: test sphinx


test:
	python tests.py

sphinx:
	make -f Makefile clean -C doc html

develop:
	python setup.py develop