all: develop test docs


test:
	nosetests --with-coverage --cover-package=glypy --cover-html --cover-html-dir=test_reports --logging-level=DEBUG -v --with-id

testall:
	nosetests --traverse-namespace --with-coverage --cover-package=glypy,glypy.plot --cover-html --cover-html-dir=test_reports --logging-level=DEBUG -v --with-id

retest:
	nosetests --traverse-namespace --cover-package=glypy,glypy.plot --logging-level=DEBUG -v --with-id --failed
clean:
	@rm  glypy/*/*.pyd
	@rm  glypy/*/*.so

sphinx:
	python setup.py build_sphinx

develop: clean
	python setup.py develop

serve-docs:
	python -m webbrowser -n docs\build\html\index.html
