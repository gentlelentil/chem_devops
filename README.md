SMILES2Mol

Chemistry - DevOps project build initial

Setup simple chemistry web application allowing user to input SMILES string to see an image of molecule.

User must login to use the application.

Users stored in simple SQLite database.

Dockerfile built with:

$ docker build -t my-app .

Docker built on linux/amd64 exposed on localhost:5000

Either run through Docker Desktop or:
$ docker run --platform linux/amd64 -d ghcr.io/gentlelentil/chem_devops_docker
