import os
from flask import Flask
from .extensions import db, migrate

def create_app():
    app = Flask(__name__)

    app.config.from_object('app.default_settings')
    app.config.from_envvar('DISCOVEROME_SETTINGS')

    try:
        os.makedirs(app.instance_path)
    except OSError:
        pass

    db.init_app(app)
    migrate.init_app(app, db)

    return app