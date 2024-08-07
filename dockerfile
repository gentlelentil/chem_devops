# Use Miniconda as the base image
FROM continuumio/miniconda3


# environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# set work directory
WORKDIR /app

COPY environment.yml /app/


#create conda environment
RUN conda env create -f environment.yml

# Activate the Conda environment
SHELL ["conda", "run", "-n", "smilesapp_dev", "/bin/bash", "-c"]

COPY . /app/

RUN mkdir -p /app/instance

#copy database files
COPY instance/users.db instance/users.db

#expose port
EXPOSE 5000

# run server
CMD ["conda", "run", "--no-capture-output", "-n", "smilesapp_dev", "flask", "run", "--host=0.0.0.0", "-p", "5000"]
