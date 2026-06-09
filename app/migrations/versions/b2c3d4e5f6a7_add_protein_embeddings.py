"""add protein_embeddings table for RAG

Revision ID: b2c3d4e5f6a7
Revises: a1b2c3d4e5f6
Create Date: 2026-06-08 00:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql
from pgvector.sqlalchemy import Vector

revision = 'b2c3d4e5f6a7'
down_revision = 'a1b2c3d4e5f6'
branch_labels = None
depends_on = None


def upgrade():
    op.execute('CREATE EXTENSION IF NOT EXISTS vector')
    op.create_table(
        'protein_embeddings',
        sa.Column('id', postgresql.UUID(as_uuid=True), primary_key=True),
        sa.Column('uniprot', sa.String, sa.ForeignKey('core.proteins.uniprot'), nullable=False),
        sa.Column('chunk_text', sa.Text, nullable=False),
        sa.Column('embedding', Vector(768), nullable=False),
        schema='core',
    )
    op.create_index(
        'ix_core_protein_embeddings_uniprot',
        'protein_embeddings',
        ['uniprot'],
        schema='core',
    )


def downgrade():
    op.drop_index('ix_core_protein_embeddings_uniprot', table_name='protein_embeddings', schema='core')
    op.drop_table('protein_embeddings', schema='core')
