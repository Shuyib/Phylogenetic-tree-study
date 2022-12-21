install:
	#install requirements
	pip --no-cache-dir install --upgrade pip &&\
		pip --no-cache-dir install -r requirements.txt	
format:
	#format code
	black *.py utils/*.py testing/*.py
docstring:
	# format docstring
	pyment -w -o numpydoc *.py
lint:
	#flake8 or #pylint
	pylint --disable=R,C --errors-only *.py utils/*.py testing/*.py
test:
	#test
	python -m pytest testing/*.py
run_script:
	# run script
	python multifastaloader.py
	python load_multifasta_metrics.py
build:
	# build the container: More important for the CI/CD
	sudo docker build -t phylo-exp .
run_test:
	# linting Dockerfile
	sudo docker run --rm -i hadolint/hadolint < Dockerfile
run:
	# run the container
	sudo docker run -it -p 8888:8888 --rm phylo-exp:latest
remove:
	# remove file
	rm updated_data/sequences.fasta 
	rm updated_data/sequence_metrics.csv
deploy:
	# add step to deploy to cloud provider if any
	echo "todo"

all: install docstring format lint test build run_test run
