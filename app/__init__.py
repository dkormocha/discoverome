import os
from flask import Flask
from .extensions import db, migrate
from .config import DevelopmentConfig, ProductionConfig, TestingConfig
from prometheus_flask_exporter import PrometheusMetrics


def create_app():

    app = Flask(__name__)

    env = os.getenv("APP_ENV", "development")

    config_map = {
        "development": DevelopmentConfig,
        "production": ProductionConfig,
        "testing": TestingConfig,
    }

    app.config.from_object(config_map[env])

    db.init_app(app)
    migrate.init_app(app, db)

    metrics = PrometheusMetrics(app)
    
    # register blueprints
    from app.commands import databp
    app.register_blueprint(databp)

    return app