all: install

.PHONY: install test docs lint

install:
	pip install .
	rm -rf build *.egg-info

test:
	pip install .[dev]
	pytest

docs:
	pip install .[docs]
	cd docs;\
	make clean; \
	make html

lint:
	ruff check tcrconvert