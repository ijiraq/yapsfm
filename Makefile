.DEFAULT_GOAL := test

setup:
	@sudo apt-get install build-essential git-core python python-nose python-numpy python-pip python-matplotlib python-scipy python-pyfits
	@sudo pip install pypng

test:
	nosetests tests
