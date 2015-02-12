all: test sphinx


test:
	coverage run tests.py
	coverage html -d test_reports

sphinx:
	make -f Makefile clean -C doc html

develop:
	python setup.py develop