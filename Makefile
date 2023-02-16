.PHONY: all lint test test-cov install dev clean distclean

PYTHON ?= python

all: ;

lint:
	q2lint
	flake8

test: all
	py.test

test-cov: all
	py.test --cov=rescript

install: all
	$(PYTHON) setup.py install
	pip install "ncbi-datasets-pylib==12.6.0"

dev: all
	pip install -e .

clean: distclean

distclean: ;
