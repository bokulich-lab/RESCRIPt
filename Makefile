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
	$(PYTHON) -m pip install ncbi-datasets-pylib

dev: all
	pip install -e .

clean: distclean

distclean: ;
