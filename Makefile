install:
	# install commands
	pip --no-cache-dir install --upgrade pip &&\
		pip --no-cache-dir install -r requirements.txt


build:
	# build the container: More important for the CI/CD
	docker build -t phylo-exp .
	
run:
	# run the container
	docker run -it -p 9999:9999 --rm phylo-exp:latest

all: install build run
