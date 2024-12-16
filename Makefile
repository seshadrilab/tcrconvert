all: install

.PHONY: install test docs lint

install:
	pip install .
	rm -rf build *.egg-info

test:
	pytest

docs:
	cd docs;\
	make html

lint:
	ruff check tcrconvert