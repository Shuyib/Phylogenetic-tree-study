Here is the updated Makefile with `$PIP` and `$PYTHON` variables used in the appropriate sections, including for `black`, `pyment`, and `pylint`:

```makefile
# .ONESHELL tells make to run each recipe line in a single shell
.ONESHELL:

# .DEFAULT_GOAL tells make which target to run when no target is specified
.DEFAULT_GOAL := all

# Specify python location in virtual environment
# Specify pip location in virtual environment
PYTHON := .venv/bin/python3
PIP := .venv/bin/pip3

venv/bin/activate: requirements.txt
	# create virtual environment
	python3 -m venv .venv
	# make command executable
	chmod +x .venv/bin/activate
	# activate virtual environment
	. .venv/bin/activate
  
activate:
	# activate virtual environment
	. .venv/bin/activate

install: requirements.txt # prerequisite
	# install commands using pip from the virtual environment
	$(PIP) install --no-cache-dir --upgrade pip &&\
	$(PIP) install --no-cache-dir -r requirements.txt

format: install
	# format code using black
	$(PYTHON) -m black *.py utils/*.py testing/*.py

clean: 
	# clean directory of cache files, virtual environment, and data
	rm -rf __pycache__ &&\
	rm -rf utils/__pycache__ &&\
	rm -rf testing/__pycache__ &&\
	rm -rf .pytest_cache &&\
	rm -rf .venv &&\
	rm updated_data/sequences.fasta &&\
	rm updated_data/sequence_metrics.csv

docstring: activate install
	# format docstring using pyment
	$(PYTHON) -m pyment -w -o numpydoc *.py

lint: activate install format
	# lint code using pylint
	$(PYTHON) -m pylint --disable=R,C --errors-only *.py utils/*.py testing/*.py

test: activate install format lint
	# run tests using pytest
	$(PYTHON) -m pytest testing/*.py

run_script:
	# run script: load fasta, create sequence_metrics
	# calculate cosine similarity averages
	$(PYTHON) multifastaloader.py
	$(PYTHON) load_multifasta_metrics.py
	$(PYTHON) tokenize_compare_sequences.py

docker_build: Dockerfile requirements.txt 
	# build the container: More important for the CI/CD
	sudo docker build -t phylo-exp .

docker_run_test: Dockerfile
	# lint Dockerfile
	sudo docker run --rm -i hadolint/hadolint < Dockerfile

docker_run: docker_build 
	# run the container
	sudo docker run -it -p 8888:8888 --rm phylo-exp:latest

deploy:
	# add step to deploy to cloud provider if any
	echo "todo"

# .PHONY tells make that these targets do not represent actual files
.PHONY: activate format clean lint test build run docker_build docker_run docker_push docker_clean docker_run_test

# what all steps to run
all: install docstring format lint test docker_build docker_run
```

This ensures that `black`, `pyment`, and `pylint` use the appropriate Python and pip from the virtual environment.
