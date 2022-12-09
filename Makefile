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
build:
	# build the container: More important for the CI/CD
	docker build -t phylo-exp .
run_test:
	# linting Dockerfile
	docker run --rm -i hadolint/hadolint < Dockerfile
run:
	# run the container
	docker run -it -p 9999:9999 --rm phylo-exp:latest
deploy:
	# add step to deploy to cloud provider if any
	echo "todo"

all: install docstring format lint test build run_test run
