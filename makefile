all: test docs


test:
	nosetests --with-coverage --cover-package=pygly2 --cover-html --cover-html-dir=test_reports --logging-level=CRITICAL

docs:
	make -f Makefile clean -C doc html

develop:
	python setup.py develop

serve-docs:
	python -m webbrowser -n doc\build\html\index.html
