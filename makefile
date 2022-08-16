all: develop test docs


test:
	py.test -v  ./tests --cov=glypy --cov-report html --cov-report term

retest:
	py.test -v ./tests --lf

clean:
	@rm  src/glypy/*/*.pyd
	@rm  src/glypy/*/*.so

sphinx:
	python setup.py build_sphinx

develop:
	python setup.py develop
