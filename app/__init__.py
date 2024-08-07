from flask import Flask
from flask_sqlalchemy import SQLAlchemy
from flask_login import LoginManager
from config import Config
from apscheduler.schedulers.background import BackgroundScheduler
import os

db = SQLAlchemy()
login_manager = LoginManager()
login_manager.login_view = 'auth.login'

def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)
    
    db.init_app(app)
    login_manager.init_app(app)
    # login_manager.login_view = 'auth.login'

    from .models import User  # Import the User model

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    # register the blueprint
    from .routes import auth
    app.register_blueprint(auth)

    #initialise the APScheduler
    scheduler = BackgroundScheduler()
    scheduler.add_job(cleanup_images, 'interval', hours=1)
    scheduler.start()

    with app.app_context():
        db.create_all()  # Create database tables

    return app

def cleanup_images():
    #delete images older than 1hr for space saving
    import time
    current_time = time.time()
    from app.routes import IMAGE_FOLDER
    for filename in os.listdir(IMAGE_FOLDER):
        fpath = os.path.join(IMAGE_FOLDER, filename)
        if os.path.isfile(fpath):
            file_creationtime = os.path.getmtime(fpath)
            if current_time - file_creationtime >= 3600:
                os.remove(fpath)