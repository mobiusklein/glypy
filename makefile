all: develop test docs


test:
	nosetests --with-coverage --with-timer --cover-package=glypy --cover-html --cover-html-dir=test_reports --logging-level=DEBUG -v --with-id glypy/tests/

retest:
	nosetests --cover-package=glypy --logging-level=DEBUG -v --with-id --failed glypy/tests/
clean:
	@rm  glypy/*/*.pyd
	@rm  glypy/*/*.so

sphinx:
	python setup.py build_sphinx

develop:
	python setup.py develop

serve-docs:
	python -m webbrowser -n docs\build\html\index.html
