all: install

.PHONY: install test docs lint

install:
	pip install .

test:
	pytest

docs:
	cd docs;\
	make html

lint:
	ruff check tcrconvert