SHELL := /usr/bin/env bash

#######
# Help
#######

.DEFAULT_GOAL := help
.PHONY: help

help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-30s\033[0m %s\n", $$1, $$2}'

###################
# Conda Enviroment
###################

PY_VERSION := 3.8
CONDA_ENV_NAME ?= phylo-env
ACTIVATE_ENV = source activate ./$(CONDA_ENV_NAME)

.PHONY: build-conda-env
build-conda-env: $(CONDA_ENV_NAME)  ## Build the conda environment
$(CONDA_ENV_NAME):
	conda create -p $(CONDA_ENV_NAME)  --copy -y  python=$(PY_VERSION)
	$(ACTIVATE_ENV) && pip --no-cache-dir install --upgrade pip &&
	python -s -m pip --no-cache-dir install -r requirements.txt

.PHONY: clean-conda-env
clean-conda-env:  ## Remove the conda environment and the relevant file
	rm -rf $(CONDA_ENV_NAME)
	rm -rf $(CONDA_ENV_NAME).zip

.PHONY: add-to-jupyter
add-to-jupyter: ## Register the conda environment to Jupyter
	$(ACTIVATE_ENV) && python -s -m ipykernel install --user --name $(CONDA_ENV_NAME)

.PHONY: remove-from-jupyter
remove-from-jupyter: ## Remove the conda environment from Jupyter
	jupyter kernelspec uninstall $(CONDA_ENV_NAME)
	
build:
	# build the container: More important for the CI/CD
	docker build -t phylo-exp .
	
run:
	# run the container
	docker run -it -p 9999:9999 --rm phylo-exp:latest

all: build-conda-env clean-conda-env add-to-jupyter remove-from-jupyter build run
