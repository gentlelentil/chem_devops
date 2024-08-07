# import os
# import time

# from flask import Flask, render_template, request, redirect, url_for, flash
# from flask_sqlalchemy import SQLAlchemy
# from flask_login import LoginManager, login_user, logout_user, login_required, UserMixin

# from apscheduler.schedulers.background import BackgroundScheduler

# from models import db, User  # Import the model from models.py
# #TODO add routes file to refactor code


# from wtforms import StringField, PasswordField, SubmitField
# from wtforms.validators import DataRequired, Length, EqualTo
# from flask_wtf import FlaskForm


# from rdkit import Chem
# from rdkit.Chem import Draw

# # db = SQLAlchemy()


# def create_app():

#     app = Flask(__name__)
#     # login code
#     app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY') or 'test_key'
#     app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db.sqlite'

#     db.init_app(app)
#     scheduler = BackgroundScheduler()
#     # scheduler.init_app(app)

#     login_manager = LoginManager(app)
#     login_manager.login_view = 'login'
    
#     IMAGE_FOLDER = 'static/images'
#     #in seconds
#     CLEANUP_INTERVAL = 3600

#     if not os.path.exists(IMAGE_FOLDER):
#         os.makedirs(IMAGE_FOLDER)

#     def cleanup_images():
#         now = time.time()
#         for filename in os.listdir(IMAGE_FOLDER):
#             file_path = os.path.join(IMAGE_FOLDER, filename)
#             if os.path.isfile(file_path):
#                 #older than 1 hour
#                 if os.stat(file_path).st_mtime < now - CLEANUP_INTERVAL:
#                     os.remove(file_path)
#                     print(f"Deleted {file_path}")


#     scheduler.add_job(func=cleanup_images, trigger='interval', seconds=CLEANUP_INTERVAL)
    
#     scheduler.start()

#     # class User(db.Model, UserMixin):
#     #     __tablename__ = 'user'
#     #     id = db.Column(db.Integer, primary_key=True)
#     #     username = db.Column(db.String(150), unique=True, nullable=False)
#     #     password = db.Column(db.String(150), nullable=False)



#     class RegistrationForm(FlaskForm):
#         username = StringField('Username', validators=[DataRequired(), Length(min=4, max=150)])
#         password = PasswordField('Password', validators=[DataRequired(), Length(min=4, max=150)])
#         confirm_password = PasswordField('Confirm Password', validators=[DataRequired(), EqualTo('password')])
#         submit = SubmitField('Register')

#     class LoginForm(FlaskForm):
#         username = StringField('Username', validators=[DataRequired(), Length(min=4, max=150)])
#         password = PasswordField('Password', validators=[DataRequired(), Length(min=4, max=150)])
#         submit = SubmitField('Login')

#     @login_manager.user_loader
#     def load_user(user_id):
#         return User.query.get(int(user_id))

#     @app.route('/register', methods=['GET', 'POST'])
#     def register():
#         form = RegistrationForm()
#         error_message = None

#         if form.validate_on_submit():
#             existing_user = User.query.filter_by(username=form.username.data).first()
#             if existing_user:
#                 error_message = 'Username already exists! Please choose a different username'
#                 # flash('Username already exists! Please choose a different username', 'danger')
#                 # return redirect(url_for('register.html', form=form))
#             else:
#                 hashed_password = form.password.data
#                 user = User(username=form.username.data, password=hashed_password)
#                 db.session.add(user)
#                 db.session.commit()
#                 flash('Account created successfully!', 'success')
#                 return redirect(url_for('login'))
#         return render_template('register.html', form=form, error_message=error_message)

#     @app.route('/login', methods=['GET', 'POST'])
#     def login():
#         form = LoginForm()
#         if form.validate_on_submit():
#             user = User.query.filter_by(username=form.username.data).first()
#             if user and user.password == form.password.data:
#                 login_user(user)
#                 flash('Login successful!', 'success')
#                 return redirect(url_for('index'))
#             else:
#                 flash('Login unsuccessful. Please check username and password', 'danger')
#         return render_template('login.html', form=form)


#     @app.route('/logout')
#     @login_required
#     def logout():
#         logout_user()
#         return redirect(url_for('login'))

#     if not os.path.exists('static/images'):
#         os.makedirs('static/images')

#     @app.route('/', methods=['GET', 'POST'])
#     @login_required
#     def index():
#         if request.method == 'POST':
#             smiles = request.form['smiles']
#             if smiles:
#                 try:
#                     #smiles to mol
#                     mol = Chem.MolFromSmiles(smiles)
#                     if mol is None:
#                         raise ValueError('Invalid SMILES string')
                    
#                     img_filename = f'molecule_{int(time.time())}.png'
#                     img_path = os.path.join(IMAGE_FOLDER, img_filename)

#                     Draw.MolToFile(mol, img_path)

#                     return render_template('index.html', img_path=img_path, smiles=smiles)
#                 except Exception as e:
#                     return render_template('index.html', error=str(e))
#         return render_template('index.html')
    
#     @app.teardown_appcontext
#     def shutdown_session(exception=None):
#         if scheduler.running:
#             scheduler.shutdown(wait=False)

#     return app

# #setup 

# app = create_app()

# with app.app_context():
#     db.create_all()
#     # scheduler.start()

# # if __name__ == '__main__':
# #     # try:
# #     app.run(debug=True)
# #     #     db.create_all()

#     # finally:
#     #     scheduler.shutdown()