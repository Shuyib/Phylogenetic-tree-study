# Use latest Python runtime as a parent image
FROM python:3.5-slim-buster

# Meta-data
LABEL maintainer="Shuyib" \
      description="Comparing different ways of constructing phylogenetic trees"
      
# Set the working directory to /app
WORKDIR /app

# ensures that the python output is sent to the terminal without buffering
ENV PYTHONUNBUFFERED=TRUE

# Copy the current directory contents into the container at /app
COPY . /app

# create a virtual environment and activate it
RUN python3 -m venv phylo-env

# activate virtual environment
CMD source phylo-env/bin/activate

# Install the required libraries
RUN pip --no-cache-dir install --upgrade pip &&\
		pip --no-cache-dir install -r requirements.txt

# Make port 1111 available to the world outside this container
EXPOSE 1111

# Create mountpoint
VOLUME /app

# Run jupyter when container launches
CMD ["jupyter", "notebook", "--ip='*'", "--port=1111", "--no-browser", "--allow-root"]

