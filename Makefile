.DEFAULT_GOAL := test

setup:
	@sudo apt-get install build-essential git-core python python-nose python-numpy python-scipy python-pyfits

test:
	nosetests tests
