all: develop test docs


test:
	py.test -v  glypy --cov=glypy --cov-report html --cov-report term

retest:
	py.test -v glypy --lf

clean:
	@rm  glypy/*/*.pyd
	@rm  glypy/*/*.so

sphinx:
	python setup.py build_sphinx

develop:
	python setup.py develop

serve-docs:
	python -m webbrowser -n docs/build/html/index.html
