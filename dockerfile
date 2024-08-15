# Use Miniconda as the base image
FROM continuumio/miniconda3


# environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# set work directory
WORKDIR /app

COPY environment.yml /app/


#create conda environment and remove cache
RUN conda env create -f environment.yml && \
    conda clean --all --yes

# Activate the Conda environment
SHELL ["conda", "run", "-n", "smilesapp_dev", "/bin/bash", "-c"]

COPY . /app/

RUN mkdir -p /app/instance

#copy database files not for initial instance
# COPY instance/users.db instance/users.db

#expose port
EXPOSE 5000

RUN pip install gunicorn

# run server with gunicorn for production
# CMD ["conda", "run", "--no-capture-output", "-n", "smilesapp_dev", "flask", "run", "--host=0.0.0.0", "-p", "5000"]
CMD ["conda", "run", "--no-capture-output", "-n", "smilesapp_dev", "gunicorn", "run", "-w", "4", "-b", "0.0.0.0:5000", "run:app"]
