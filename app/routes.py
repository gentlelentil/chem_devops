from flask import Blueprint, render_template, redirect, url_for, flash, request, jsonify, send_file
from flask import after_this_request
from flask_login import login_user, login_required, logout_user, current_user
import requests
from werkzeug.utils import secure_filename
import zipfile
import uuid
import shutil

from . import db

from rdkit import Chem
from rdkit.Chem import Draw
import os
import time
import json

from .models import User, RegistrationForm, LoginForm

auth = Blueprint('auth', __name__)

IMAGE_FOLDER = os.path.join('app', 'static', 'images')


def extract_smiles_from_sdf(file_path):
    smiles_list = []
    supplier = Chem.SDMolSupplier(file_path)
    for mol in supplier:
        if mol is not None:
            smiles_list.append(Chem.MolToSmiles(mol))
    return smiles_list

#homepage
@auth.route('/')
def home():
    if current_user.is_authenticated:
        return redirect(url_for('auth.smiles_generator'))
    else:
        return redirect(url_for('auth.login'))


# @app.route('/register', methods=['GET', 'POST'])
@auth.route('/register', methods=['GET', 'POST'])
def register():
    form = RegistrationForm()
    error_message = None

    if form.validate_on_submit():
        existing_user = User.query.filter_by(username=form.username.data).first()
        if existing_user:
            error_message = 'Username already exists! Please choose a different username'
            # flash('Username already exists! Please choose a different username', 'danger')
            # return redirect(url_for('register.html', form=form))
        else:
            hashed_password = form.password.data
            user = User(username=form.username.data, password=hashed_password)
            db.session.add(user)
            db.session.commit()
            flash('Account created successfully!', 'success')
            return redirect(url_for('auth.login'))
    return render_template('register.html', form=form, error_message=error_message, title='Registration')

# @app.route('/login', methods=['GET', 'POST'])
@auth.route('/login', methods=['GET', 'POST'])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user = User.query.filter_by(username=form.username.data).first()
        if user and user.password == form.password.data:
            login_user(user)
            flash('Login successful!', 'success')
            # return redirect(url_for('auth.index'))
            return redirect(url_for('auth.smiles_generator'))
        else:
            flash('Login unsuccessful. Please check username and password', 'danger')
    return render_template('login.html', form=form, title='Login')

# @app.route('/logout')
@auth.route('/logout')
@login_required
def logout():
    logout_user()
    return redirect(url_for('auth.login'))

# @app.route('/', methods=['GET', 'POST'])
# @auth.route('/', methods=['GET', 'POST'])
# @login_required
# def index():
#     if request.method == 'POST':
#         smiles = request.form['smiles']
#         if smiles:
#             try:
#                 #smiles to mol
#                 mol = Chem.MolFromSmiles(smiles)
#                 if mol is None:
#                     raise ValueError('Invalid SMILES string')
                
#                 img_filename = f'molecule_{int(time.time())}.png'
#                 img_path = os.path.join(IMAGE_FOLDER, img_filename)

#                 Draw.MolToFile(mol, img_path)
#                 print(img_path)

#                 return render_template('index.html', img_filename=img_filename, smiles=smiles)
#             except Exception as e:
#                 return render_template('index.html', error=str(e))
#     return render_template('index.html')

#smiles route
@auth.route('/smiles-generator', methods=['GET', 'POST'])
@login_required
def smiles_generator():
    if request.method == 'POST':
        smiles = request.form['smiles']
        if smiles:
            try:
                #smiles to mol
                mol = Chem.MolFromSmiles(smiles)
                if mol is None:
                    raise ValueError('Invalid SMILES string')
                
                img_filename = f'molecule_{int(time.time())}.png'
                img_path = os.path.join(IMAGE_FOLDER, img_filename)

                Draw.MolToFile(mol, img_path)
                print(img_path)

                return render_template('smiles-generator.html', img_filename=img_filename, smiles=smiles, title='SMILES Generator')
            except Exception as e:
                return render_template('smiles-generator.html', error=str(e))
    return render_template('smiles-generator.html', title='SMILES Generator')

#model url
ADMET_MODEL_URL = 'http://ml_model:5001/predict'
dev_ADMET_MODEL_URL = 'http://0.0.0.0:5001/predict'

model_categories = {
    'Absorption': {
        'oral_abs_class': 'Oral Absorption',
        'hia_class': 'Human Intestinal Absorption'
    },
    'Metabolism': {
        'ML_input_p450-cyp3a4': 'CYP3A4 Inhibition',
        'ML_input_p450-cyp2c19': 'CYP2C19 Inhibition',
        'ML_input_p450-cyp1a2': 'CYP1A2 Inhibition',
        'ML_input_p450-cyp2d6': 'CYP2D6 Inhibition',
        'ML_input_p450-cyp2c9': 'CYP2C9 Inhibition'
    },

    'Distribution': {
        'bbb_class': 'Blood Brain Barrier Permeability'
    },
    'Toxicity': {
        'crl_toxicity_class': 'CRL Toxicity',
        'hep_g2_toxicity_class': 'HepG2 Toxicity',
        'hek_toxicity_class': 'HEK Toxicity',
        'hacat_toxicity_class': 'HaCaT Toxicity',
        'nih_toxicity_class': 'NIH Toxicity',
        'ames_mutagenicity_class': 'Ames Mutagenicity'
    },
    'Cardiotoxicity': {
        'herg_blockers_class': 'hERG Blockers'
    },
}


#prediction route
@auth.route('/admetpredictor', methods=['GET', 'POST'])
@login_required
def admetpredictor():

    if request.method == 'POST':
        smiles_string = request.form.get('smiles_string')
        
        
        if smiles_string:
            try:
                response = requests.post(ADMET_MODEL_URL, json={'smiles': smiles_string})
            except requests.exceptions.RequestException as e:
                response = requests.post(dev_ADMET_MODEL_URL, json={'smiles': smiles_string})

            if response.status_code == 200:
                result = response.json()

                if 'prediction' in result:
                    categorized_results = {category: [] for category in model_categories}

                    for model_name, prediction in result['prediction'].items():
                        for category, subcategories in model_categories.items():
                            if model_name in subcategories:
                                categorized_results[category].append((subcategories[model_name], prediction))
                   
                    # format_results = [(key, value) for key, value in result['prediction'].items()]



                else:
                    format_results = 'No prediction available'

                #admet_result because it is the FUNCTION NAME NOT THE HTML FILENAME

                # return render_template('admet-predictor.html', prediction=format_results, title='ADMET Predictor')
                return render_template('admet-predictor.html', prediction=categorized_results, title='ADMET Predictor')

                # return redirect (url_for('auth.admet_result', result=result))
            else:
                return 'Failed to get prediction', response.status_code
            


    return render_template('admet-predictor.html', title='ADMET Predictor', prediction=None)

@auth.route('/alphafold3', methods=['GET', 'POST'])
@login_required
def alphafold3():
    if request.method == 'POST':
        # Collect input data from the form
        name = request.form.get('name')
        if not name:
            return jsonify({"error": "Name is required."}), 400

        # Dynamically determine the number of sequences from the form keys
        # Collect sequences from the form
        sequence_keys = [key for key in request.form.keys() if key.startswith('sequence_')]
        sequences = [request.form.get(key) for key in sequence_keys if len(request.form.get(key).strip()) > 0]  # Filter empty strings

        if not sequences or any(not seq for seq in sequences):
            return jsonify({"error": "At least one valid sequence is required."}), 400
        if len(sequences) > 3:
            return jsonify({"error": "You can only submit up to 3 sequences."}), 400

        # Handle uploaded SDF file
        sdf_file = request.files.get('sdf_file')

        # Collect ligand SMILES strings from the form
        ligand_keys = [key for key in request.form.keys() if key.startswith('ligand_smiles_')]
        lig_check = [request.form.get(key) for key in ligand_keys if len(request.form.get(key).strip()) > 0]  # Filter empty strings

        # If ligands are provided, check the number of ligands
        if ligand_keys and len(ligand_keys) > 5:
            return jsonify({"error": "You can only submit up to 10 ligands."}), 400

        if sdf_file and lig_check:
            return jsonify({"error": "Please provide either an SDF file or manual ligands, not both."}), 400

        # Prepare files lists
        sequence_files = []
        ligand_files = []

        # Handling sequence files
        for i, seq in enumerate(sequences, start=1):
            sequence_filename = f"sequence_{i}.txt"
            sequence_files.append((sequence_filename, seq))

        # Handle SDF upload
        if sdf_file:
            ligands = []
            if not sdf_file.filename.endswith('.sdf'):
                return jsonify({"error": "Invalid file type. Please upload an SDF file."}), 400
            
            filename = secure_filename(sdf_file.filename)
            sdf_path = os.path.join("/tmp", filename)
            sdf_file.save(sdf_path)

            try:
                supplier = Chem.SDMolSupplier(sdf_path)
                ligands = [Chem.MolToSmiles(mol) for mol in supplier if mol is not None]
            except Exception as e:
                return jsonify({"error": f"Failed to parse SDF file: {str(e)}"}), 400
            finally:
                os.remove(sdf_path)  # Clean up the temporary file

            if not ligands:
                return jsonify({"error": "No valid ligands found in the SDF file."}), 400

        # Handle manual ligand input
        if lig_check:
            ligands = [request.form.get(key) for key in ligand_keys if request.form.get(key)]

        if not ligands:
            return jsonify({"error": "At least one ligand is required (SDF or manual input)."}), 400

        # Handling ligand files
        for i, ligand in enumerate(ligands, start=1):
            ligand_filename = f"ligand_{i}.txt"
            ligand_files.append((ligand_filename, ligand))

        # Prepare the JSON data
        json_data = []
        for seq_idx, sequence in enumerate(sequences, start=1):
            for lig_idx, ligand in enumerate(ligands, start=1):
                data = {
                    "name": name,
                    "sequences": [
                        {
                            "protein": {
                                "id": f"sequence_{seq_idx}",
                                "sequence": sequence
                            }
                        },
                        {
                            "ligand": {
                                "id": f"ligand_{lig_idx}",
                                "smiles": ligand
                            }
                        }
                    ],
                    "modelSeeds": [1],
                    "dialect": "alphafold3",
                    "version": 1
                }
                json_data.append(data)


        # Generate a unique folder for the job based on the job name and a UUID
        unique_id = uuid.uuid4().hex
        job_folder = os.path.join("/home/nathaniel/Desktop/flask/app/static/af3_generated_inputs", f"{name}_{unique_id}")
        os.makedirs(job_folder, exist_ok=True)


        # Save each JSON to a file, creating a unique filename for each combination of sequence and ligand
        for idx, data in enumerate(json_data, start=1):
            json_filename = f"{name}_seq_{data['sequences'][0]['protein']['id']}_lig_{data['sequences'][1]['ligand']['id']}.json"
            json_filepath = os.path.join(job_folder, json_filename)
            with open(json_filepath, 'w') as json_file:
                json.dump(data, json_file, indent=2)

        # # Save generated files into the unique job folder
        # for i, sequence in enumerate(sequences, start=1):
        #     sequence_filename = os.path.join(job_folder, f"sequence_{i}.txt")
        #     with open(sequence_filename, 'w') as seq_file:
        #         seq_file.write(sequence)

        # for i, ligand in enumerate(ligands, start=1):
        #     ligand_filename = os.path.join(job_folder, f"ligand_{i}.json")
        #     with open(ligand_filename, 'w') as lig_file:
        #         json.dump({"ligand_smiles": ligand}, lig_file, indent=2)

        # Create a subdirectory for the ZIP contents
        job_subdir = os.path.join(job_folder, f"{name}_inputs")
        os.makedirs(job_subdir, exist_ok=True)

        # Move generated files into the job-specific subdirectory
        for file in os.listdir(job_folder):
            file_path = os.path.join(job_folder, file)
            if os.path.isfile(file_path):  # Skip directories
                shutil.move(file_path, os.path.join(job_subdir, file))

        # Create a ZIP file containing the job folder
        zip_filename = f"{job_folder}.zip"
        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(job_folder):
                for file in files:
                    zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), job_folder))

        # Delete the job folder and its files (but keep the ZIP)
        for root, _, files in os.walk(job_folder):
            for file in files:
                os.remove(os.path.join(root, file))
        # Now delete everything in the output directory recursively (including the job_subdir)
        shutil.rmtree(job_folder)  # Recursively removes the output_dir and everything inside it

        # try:
        #     os.rmdir(job_folder)  # Remove the empty job folder
        # except OSError as e:
        #     print(f"Error deleting directory {job_folder}: {e.strerror}")

        # Redirect to the success page with the ZIP filename
        return redirect(url_for('auth.success', zip_filename=os.path.basename(zip_filename)))


    return render_template('alphafold3.html')

@auth.route('/success')
@login_required
def success():
    zip_filename = request.args.get('zip_filename')
    return render_template('success.html', zip_filename=zip_filename)

@auth.route('/download_zip/<filename>')
@login_required
def download_zip(filename):
    zip_path = os.path.join('/home/nathaniel/Desktop/flask/app/static/af3_generated_inputs', filename)

    @after_this_request
    def remove_file(response):
        try:
            os.remove(zip_path)  # Delete the ZIP file after serving it
        except Exception as e:
            print(f"Error deleting ZIP file {zip_path}: {e}")
        return response

    return send_file(zip_path, as_attachment=True)
