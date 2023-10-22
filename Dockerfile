# Use latest Python runtime as a parent image
FROM python:3.8.16-slim-buster

# Meta-data
LABEL maintainer="Shuyib" \
      description="Comparing different ways of constructing phylogenetic trees"
      
# Set the working directory to /app
WORKDIR /app

# ensures that the python output is sent to the terminal without buffering
ENV PYTHONUNBUFFERED=TRUE

# Copy the current directory contents into the container at /app
COPY . /app

# create a virtual environment, activate it and install the required libraries
RUN python3 -m venv phylo-env &&\
            . phylo-env/bin/activate &&\
            pip --no-cache-dir install --upgrade pip &&\
            pip --no-cache-dir install --requirement requirements.txt

# Make port 8888 available to the world outside this container
EXPOSE 8888

# Create mountpoint
VOLUME /app

# Run jupyter when container launches
CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--port=8888", "--no-browser", "--allow-root"]



