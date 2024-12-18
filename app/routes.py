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

# Function to generate the list of chain letters based on the numeric value
def generate_chain_letters(chain_count):
    # Chain letters: ['A', 'B', 'C', 'D', ...]
    return [chr(65 + i) for i in range(chain_count)]  # 65 is the ASCII value for 'A'

# @auth.route('/alphafold3', methods=['GET', 'POST'])
# @login_required
# def alphafold3():
#     if request.method == 'POST':
#         # Collect input data from the form
#         name = request.form.get('name')
#         if not name:
#             return jsonify({"error": "Name is required."}), 400

#         # Dynamically collect sequences and their respective chain counts
#         sequence_keys = [key for key in request.form.keys() if key.startswith('sequence_')]
#         sequences = [request.form.get(key).strip() for key in sequence_keys if len(request.form.get(key).strip()) > 0]
        
#         chain_keys = [key for key in request.form.keys() if key.startswith('chain_')]
#         chains = [int(request.form.get(key).strip()) for key in chain_keys if len(request.form.get(key).strip()) > 0]

#         # Validate input
#         if not sequences or any(not seq for seq in sequences):
#             return jsonify({"error": "At least one valid sequence is required."}), 400
#         if len(sequences) != len(chains):
#             return jsonify({"error": "Each sequence must have an associated chain count."}), 400

#         # Handle ligand inputs
#         sdf_file = request.files.get('sdf_file')
#         ligand_keys = [key for key in request.form.keys() if key.startswith('ligand_smiles_')]
#         ligands = [request.form.get(key).strip() for key in ligand_keys if len(request.form.get(key).strip()) > 0]

#         # Ensure maximum ligand count is not exceeded but only for input
#         if ligands:
#             if len(ligands) > 10:
#                 return jsonify({"error": "You can only submit up to 10 ligands."}), 400

#         if sdf_file and ligands:
#             return jsonify({"error": "Please provide either an SDF file or manual ligands, not both."}), 400

#         if sdf_file:
#             if not sdf_file.filename.endswith('.sdf'):
#                 return jsonify({"error": "Invalid file type. Please upload an SDF file."}), 400

#             filename = secure_filename(sdf_file.filename)
#             sdf_path = os.path.join("/tmp", filename)
#             sdf_file.save(sdf_path)

#             try:
#                 supplier = Chem.SDMolSupplier(sdf_path)
#                 ligands = [Chem.MolToSmiles(mol) for mol in supplier if mol is not None]
#             except Exception as e:
#                 return jsonify({"error": f"Failed to parse SDF file: {str(e)}"}), 400
#             finally:
#                 os.remove(sdf_path)  # Clean up the temporary file

#             if not ligands:
#                 return jsonify({"error": "No valid ligands found in the SDF file."}), 400

#         if not ligands:
#             return jsonify({"error": "At least one ligand is required (SDF or manual input)."}), 400



#         # Prepare the output folder
#         job_folder = os.path.join("/home/nathaniel/Desktop/flask/app/static/af3_generated_inputs", f"{name}_{uuid.uuid4().hex}")
#         os.makedirs(job_folder, exist_ok=True)

#         # Generate JSON files for all combinations of sequences and ligands
#         for seq_idx, (sequence, chain_count) in enumerate(zip(sequences, chains)):
#             # Dynamically generate chain IDs (A, B, C, ...)
#             chain_ids = [chr(65 + i) for i in range(chain_count)]

#             for lig_idx, ligand in enumerate(ligands):
#                 # Dynamically assign a ligand ID (after the chains)
#                 ligand_id = chr(65 + chain_count)  # Start after the chain IDs

#                 # Create the JSON object
#                 json_data = {
#                     "name": name,
#                     "sequences": [{
#                         "protein": {
#                             "id": chain_ids,
#                             "sequence": sequence
#                         }
#                     }],
#                     "ligands": [{
#                         "ligand": {
#                             "id": ligand_id,
#                             "smiles": ligand
#                         }
#                     }],
#                     "modelSeeds": [1],
#                     "dialect": "alphafold3",
#                     "version": 1
#                 }

#                 # Save the JSON file
#                 json_filename = f"{name}_seq{seq_idx + 1}_ligand{lig_idx + 1}_input.json"
#                 json_filepath = os.path.join(job_folder, json_filename)
#                 with open(json_filepath, 'w') as json_file:
#                     json.dump(json_data, json_file, indent=2)

#         # Create a ZIP file containing all the generated JSON files
#         zip_filename = f"{job_folder}.zip"
#         with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
#             for root, _, files in os.walk(job_folder):
#                 for file in files:
#                     zipf.write(os.path.join(root, file), os.path.relpath(os.path.join(root, file), job_folder))

#         # Clean up the temporary job folder
#         shutil.rmtree(job_folder)

#         # Redirect to the success page
#         return redirect(url_for('auth.success', zip_filename=os.path.basename(zip_filename)))

#     return render_template('alphafold3.html')


@auth.route('/alphafold3', methods=['GET', 'POST'])
@login_required
def alphafold3():
    if request.method == 'POST':
        # Collect input data from the form
        name = request.form.get('name')
        if not name:
            return jsonify({"error": "Name is required."}), 400

        # Dynamically collect sequences and their respective chain counts
        sequence_keys = [key for key in request.form.keys() if key.startswith('sequence_')]
        sequences = [request.form.get(key).strip() for key in sequence_keys if len(request.form.get(key).strip()) > 0]
        
        chain_keys = [key for key in request.form.keys() if key.startswith('chain_')]
        chains = [int(request.form.get(key).strip()) for key in chain_keys if len(request.form.get(key).strip()) > 0]

        # Validate input
        if not sequences or any(not seq for seq in sequences):
            return jsonify({"error": "At least one valid sequence is required."}), 400
        if len(sequences) != len(chains):
            return jsonify({"error": "Each sequence must have an associated chain count."}), 400

        # Handle ligand inputs
        sdf_file = request.files.get('sdf_file')
        ligand_keys = [key for key in request.form.keys() if key.startswith('ligand_smiles_')]
        ligands = [request.form.get(key).strip() for key in ligand_keys if len(request.form.get(key).strip()) > 0]

        # Handle the new dropdown for duplicate ligands per sequence
        duplicate_ligands_per_sequence = int(request.form.get('ligand_duplicates', 1))  # Default to 1 if not provided
        if duplicate_ligands_per_sequence < 1 or duplicate_ligands_per_sequence > 4:
            return jsonify({"error": "Duplicate ligands per sequence must be between 1 and 4."}), 400

        if ligands:
            if len(ligands) > 10:
                return jsonify({"error": "You can only submit up to 10 ligands."}), 400

        if sdf_file and ligands:
            return jsonify({"error": "Please provide either an SDF file or manual ligands, not both."}), 400

        if sdf_file:

            # Handle the new dropdown for duplicate ligands per sequence
            duplicate_ligands_per_sequence = int(request.form.get('sdf_duplicates', 1))  # Default to 1 if not provided
            if duplicate_ligands_per_sequence < 1 or duplicate_ligands_per_sequence > 4:
                return jsonify({"error": "Duplicate ligands per sequence must be between 1 and 4."}), 400

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
                os.remove(sdf_path)

            if not ligands:
                return jsonify({"error": "No valid ligands found in the SDF file."}), 400

        if not ligands:
            return jsonify({"error": "At least one ligand is required (SDF or manual input)."}), 400

        # Prepare the output folder
        job_folder = os.path.join("/home/nathaniel/Desktop/flask/app/static/af3_generated_inputs", f"{name}_{uuid.uuid4().hex}")
        os.makedirs(job_folder, exist_ok=True)

        # Create lists for sequence IDs and ligand SMILES to be written into text files
        sequence_ids = [f"SEQ{idx + 1}: {seq}" for idx, seq in enumerate(sequences)]
        ligand_ids_and_smiles = [f"LIG{idx + 1}: {ligand}" for idx, ligand in enumerate(ligands)]

        print(f"DEBUG number of duplicate ligands: {duplicate_ligands_per_sequence}")

        # Generate JSON files for all combinations of sequences and ligands
        for seq_idx, (sequence, chain_count) in enumerate(zip(sequences, chains)):
            chain_ids = [chr(65 + i) for i in range(chain_count)]


            for lig_idx, ligand in enumerate(ligands):

                jsonname = f"{name}_seq{seq_idx + 1}_ligand{lig_idx + 1}"
                json_filename = f"{name}_seq{seq_idx + 1}_ligand{lig_idx + 1}_input.json"

                # Initialize the sequence data
                sequence_data = {
                    "protein": {
                        "id": chain_ids,
                        "sequence": sequence
                    }
                }

                # Create a list for ligands with the proper structure
                ligands_list = []

                # Add the requested number of duplicates for the current ligand
                for duplicate_idx in range(duplicate_ligands_per_sequence):
                    ligand_id = chr(65 + chain_count + duplicate_idx)
                    ligand_entry = {
                        "ligand": {
                            "id": ligand_id,
                            "smiles": ligand
                        }
                    }
                    ligands_list.append(ligand_entry)  # Append ligand to list

                # Create the full JSON data structure
                json_data = {
                    "name": jsonname,
                    "sequences": [
                        sequence_data,
                        *ligands_list  # Include protein and ligands as separate entries
                    ],
                    "modelSeeds": [1],
                    "dialect": "alphafold3",
                    "version": 1
                }
                        # Save the JSON file
                json_filepath = os.path.join(job_folder, json_filename)
                with open(json_filepath, 'w') as json_file:
                    json.dump(json_data, json_file, indent=2)
 

        # Create a text file for sequence IDs and their associated sequences
        sequence_ids_filename = os.path.join(job_folder, "sequences.txt")
        with open(sequence_ids_filename, 'w') as seq_file:
            for seq_id in sequence_ids:
                seq_file.write(f"{seq_id}\n")

        # Create a text file for ligand IDs and their associated SMILES
        ligand_ids_filename = os.path.join(job_folder, "ligands.txt")
        with open(ligand_ids_filename, 'w') as lig_file:
            for lig_id_and_smiles in ligand_ids_and_smiles:
                lig_file.write(f"{lig_id_and_smiles}\n")

        # Create a ZIP file containing all the generated JSON files and text files
        zip_filename = f"{job_folder}.zip"
        with zipfile.ZipFile(zip_filename, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, _, files in os.walk(job_folder):
                for file in files:
                    zipf.write(os.path.join(root, file), 
                       os.path.relpath(os.path.join(root, file), os.path.dirname(job_folder)))

        # Clean up the temporary job folder
        shutil.rmtree(job_folder)

        # Redirect to the success page
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
