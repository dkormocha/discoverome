"""add new columns and list/pathway tables

Revision ID: a1b2c3d4e5f6
Revises: 74a7b5f67fc8
Create Date: 2026-06-03 22:00:00.000000

"""
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import UUID


revision = 'a1b2c3d4e5f6'
down_revision = '74a7b5f67fc8'
branch_labels = None
depends_on = None


def upgrade():
    # --- core.proteins: add organism column ---
    op.add_column('proteins', sa.Column('organism', sa.String(), nullable=True), schema='core')
    op.create_index('ix_core_proteins_organism', 'proteins', ['organism'], unique=False, schema='core')

    # --- core.compoundtreatments: add GP-related columns ---
    op.add_column('compoundtreatments', sa.Column('concentration_two', sa.Float(), nullable=True), schema='core')
    op.add_column('compoundtreatments', sa.Column('concentration_three', sa.Float(), nullable=True), schema='core')
    op.add_column('compoundtreatments', sa.Column('organism', sa.String(), nullable=True), schema='core')
    op.add_column('compoundtreatments', sa.Column('project', sa.String(), nullable=True), schema='core')
    op.add_column('compoundtreatments', sa.Column('data_acq', sa.String(), nullable=True), schema='core')

    # --- proteomics.foldchanges: add new columns ---
    op.add_column('foldchanges', sa.Column('group_id', UUID(as_uuid=True), nullable=True), schema='proteomics')
    op.add_column('foldchanges', sa.Column('no_peptides', sa.Float(), nullable=True), schema='proteomics')
    op.add_column('foldchanges', sa.Column('no_unique_peptides', sa.Float(), nullable=True), schema='proteomics')
    op.create_index('ix_proteomics_foldchanges_group_id', 'foldchanges', ['group_id'], unique=False, schema='proteomics')

    # --- core.targetlists ---
    op.create_table(
        'targetlists',
        sa.Column('id', UUID(as_uuid=True), server_default=sa.text('uuid_generate_v4()'), nullable=False),
        sa.Column('name', sa.String(), nullable=True),
        sa.Column('description', sa.String(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('id'),
        schema='core'
    )

    # --- core.list2protein ---
    op.create_table(
        'list2protein',
        sa.Column('id', UUID(as_uuid=True), server_default=sa.text('uuid_generate_v4()'), nullable=False),
        sa.Column('protein_id', sa.String(), nullable=True),
        sa.Column('targetlist_id', UUID(as_uuid=True), nullable=True),
        sa.ForeignKeyConstraint(['protein_id'], ['core.proteins.uniprot']),
        sa.ForeignKeyConstraint(['targetlist_id'], ['core.targetlists.id']),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('id'),
        schema='core'
    )

    # --- core.pathwaylists ---
    op.create_table(
        'pathwaylists',
        sa.Column('id', UUID(as_uuid=True), server_default=sa.text('uuid_generate_v4()'), nullable=False),
        sa.Column('name', sa.String(), nullable=True),
        sa.Column('description', sa.String(), nullable=True),
        sa.Column('database', sa.String(), nullable=True),
        sa.Column('organism', sa.String(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('id'),
        schema='core'
    )

    # --- core.protein2pathway ---
    op.create_table(
        'protein2pathway',
        sa.Column('id', UUID(as_uuid=True), server_default=sa.text('uuid_generate_v4()'), nullable=False),
        sa.Column('protein_id', sa.String(), nullable=True),
        sa.Column('pathway_id', UUID(as_uuid=True), nullable=True),
        sa.ForeignKeyConstraint(['protein_id'], ['core.proteins.uniprot']),
        sa.ForeignKeyConstraint(['pathway_id'], ['core.pathwaylists.id']),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('id'),
        schema='core'
    )


def downgrade():
    op.drop_table('protein2pathway', schema='core')
    op.drop_table('pathwaylists', schema='core')
    op.drop_table('list2protein', schema='core')
    op.drop_table('targetlists', schema='core')

    op.drop_index('ix_proteomics_foldchanges_group_id', table_name='foldchanges', schema='proteomics')
    op.drop_column('foldchanges', 'no_unique_peptides', schema='proteomics')
    op.drop_column('foldchanges', 'no_peptides', schema='proteomics')
    op.drop_column('foldchanges', 'group_id', schema='proteomics')

    op.drop_column('compoundtreatments', 'data_acq', schema='core')
    op.drop_column('compoundtreatments', 'project', schema='core')
    op.drop_column('compoundtreatments', 'organism', schema='core')
    op.drop_column('compoundtreatments', 'concentration_three', schema='core')
    op.drop_column('compoundtreatments', 'concentration_two', schema='core')

    op.drop_index('ix_core_proteins_organism', table_name='proteins', schema='core')
    op.drop_column('proteins', 'organism', schema='core')
