all: develop test docs


test:
	nosetests --traverse-namespace --with-coverage --cover-package=pygly2,pygly2.plot --cover-html --cover-html-dir=test_reports --logging-level=DEBUG -v --with-id

retest:
	nosetests --traverse-namespace --cover-package=pygly2,pygly2.plot --logging-level=CRITICAL -v --with-id --failed
clean:
	@rm  pygly2/*/*.pyd
	@rm  pygly2/*/*.so

docs:
	make -f Makefile clean -C doc html

develop: clean
	python setup.py develop

serve-docs:
	python -m webbrowser -n doc\build\html\index.html
