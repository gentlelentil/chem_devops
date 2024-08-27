from flask import Blueprint, render_template, redirect, url_for, flash, request, jsonify
from flask_login import login_user, login_required, logout_user, current_user
import requests

from . import db

from rdkit import Chem
from rdkit.Chem import Draw
import os
import time

from .models import User, RegistrationForm, LoginForm

auth = Blueprint('auth', __name__)

IMAGE_FOLDER = os.path.join('app', 'static', 'images')

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
# ADMET_MODEL_URL = 'http://0.0.0.0:5001/predict'

#prediction route
@auth.route('/admetpredictor', methods=['GET', 'POST'])
@login_required
def admetpredictor():

    if request.method == 'POST':
        smiles_string = request.form.get('smiles_string')
        
        
        if smiles_string:
            response = requests.post(ADMET_MODEL_URL, json={'smiles': smiles_string})
            if response.status_code == 200:
                result = response.json()

                if 'prediction' in result:
                   
                    format_results = [(key, value) for key, value in result['prediction'].items()]
                else:
                    format_results = 'No prediction available'

                #admet_result because it is the FUNCTION NAME NOT THE HTML FILENAME

                return render_template('admet-predictor.html', prediction=format_results, title='ADMET Predictor')

                # return redirect (url_for('auth.admet_result', result=result))
            else:
                return 'Failed to get prediction', response.status_code
            


    return render_template('admet-predictor.html', title='ADMET Predictor', prediction=None)

# @auth.route('/admetresult')
# @login_required
# def admet_result():
#     prediction = request.args.get('result')




#     return render_template('admet-result.html', prediction=prediction, title='ADMET Result')