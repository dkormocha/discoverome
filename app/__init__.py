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

    from extensions import db, migrate

    db.init_app(app)
    migrate.init_app(app, db)

    return app

#Get db connection
def get_db():
    if 'db' not in g:
        g.db = SQLAlchemy(current_app) 
    return g.db