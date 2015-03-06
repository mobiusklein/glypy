all: test docs


test:
	coverage run tests.py
	coverage html -d test_reports

docs:
	make -f Makefile clean -C doc html

develop:
	python setup.py develop

serve-docs:
	python -m webbrowser -n doc\build\html\index.html