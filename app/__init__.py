import os

from flask import Flask, g, current_app
from flask_sqlalchemy import SQLAlchemy

def create_app():

    # create and configure the app
    app = Flask(__name__)
    app.config.from_object('claymor.default_settings')
    app.config.from_envvar('CLAYMOR_SETTINGS')

    # ensure the instance folder exists
    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    from claymor.models import db, migrate

    db.init_app(app)
    migrate.init_app(app, db)

    return app

#Get db connection
def get_db():
    if 'db' not in g: #check if the 'db' key is in the flask 'g' global object 
        g.db = SQLAlchemy(current_app) #if it is not then create a SQLAlchemy database connection and point to our app 'current_app'. Connects database to our app
    return g.db