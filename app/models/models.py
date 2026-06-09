from enum import unique
from operator import index
from flask_sqlalchemy import SQLAlchemy
from flask_migrate import Migrate

from sqlalchemy import ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.orm import relationship

from sqlalchemy.dialects.postgresql import UUID
import uuid
from sqlalchemy_utils import UUIDType
from sqlalchemy.sql import text
from pgvector.sqlalchemy import Vector

from app.extensions import db




# CORE SCHEMA (shared entities)


class Protein(db.Model):
    __tablename__ = "proteins"
    __table_args__ = {"schema": "core"}

    uniprot = db.Column(db.String, primary_key=True)
    description = db.Column(db.String)
    symbol = db.Column(db.String, index=True)
    organism = db.Column(db.String, index=True)

    residues = relationship("Residue", back_populates="protein")
    synonyms = relationship("ProteinSynonym", back_populates="protein")
    chains = relationship("StructureChain", back_populates="protein")
    lists = relationship("TargetList", secondary="core.list2protein", back_populates="proteins", uselist=True, overlaps="proteins")
    pathways = relationship("PathwayList", secondary="core.protein2pathway", back_populates="proteins", uselist=True, overlaps="proteins")

class ProteinSynonym(db.Model):
    __tablename__ = "proteinsynonyms"
    __table_args__ = {"schema": "core"}

    id = db.Column( UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    uniprot = db.Column(db.String, ForeignKey("core.proteins.uniprot"))
    synonym = db.Column(db.String)
    type = db.Column(db.String)

    protein = relationship("Protein", back_populates="synonyms")


class Compound(db.Model):
    __tablename__ = "compounds"
    __table_args__ = {"schema": "core"}

    id = db.Column(db.String, primary_key=True)
    smiles = db.Column(db.String, index=True)
    image = db.Column(db.String, index=True)
    # compoundtreatments = relationship(
    #     "CompoundTreatment",
    #     back_populates="compound"
    # )
    compoundtreatments = relationship(
        "CompoundTreatment", 
        back_populates = "compound", 
        primaryjoin = "and_(Compound.id ==CompoundTreatment.compound_id)" # Needs defining because the CompoundTreatment table also has a key for the probe
    )
class CellType(db.Model):
    __tablename__ = "celltypes"
    __table_args__ = {"schema": "core"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)
    name = db.Column(db.String)
    depmapid = db.Column(db.String)
    description = db.Column(db.String)
    compoundtreatments = relationship("CompoundTreatment", back_populates = "celltype")

class Experiment(db.Model):
    __tablename__ = "experiments"
    __table_args__ = {"schema": "core"}

    id = db.Column(db.String, primary_key=True)
    description = db.Column(db.String)
    plexs = relationship("Plex", back_populates = "experiment")




class Plex(db.Model):
    __tablename__ = "plexs"
    __table_args__ = {"schema": "core"}

    id = db.Column(db.String, primary_key=True)
    experiment_id = db.Column(db.String, ForeignKey("core.experiments.id"))
    description = db.Column(db.String)

    experiment = relationship("Experiment", back_populates = "plexs")
    compoundtreatments = relationship("CompoundTreatment", back_populates = "plex", overlaps = "intensityreadings")
    intensityreadings = relationship("IntensityReading",  back_populates = "plex")
    competitionratios = relationship("CompetitionRatio",  back_populates = "plex")


class CompoundTreatment(db.Model):
    __tablename__ = "compoundtreatments"
    __table_args__ = (
        UniqueConstraint("plex_id", "samplename"),
        {"schema": "core"}
    )

    id = db.Column(UUID(as_uuid=True), primary_key=True, default=uuid.uuid4)

    plex_id = db.Column(db.String, ForeignKey("core.plexs.id"), nullable=False)
    samplename = db.Column(db.String, nullable=False)

    compound_id = db.Column(db.String, ForeignKey("core.compounds.id"), nullable=False, index=True)
    celltype_id = db.Column(UUID(as_uuid=True), ForeignKey("core.celltypes.id"), nullable=False)

    tmtchannel = db.Column(db.String)

    concentration = db.Column(db.Float)
    concentrationunits = db.Column(db.String)

    time = db.Column(db.Float)
    timeunits = db.Column(db.String)

    temperature = db.Column(db.Float)

    probe_id = db.Column(db.String, ForeignKey("core.compounds.id"), index=True)

    probeconcentration = db.Column(db.Float)
    probeconcentrationunits = db.Column(db.String)

    probetime = db.Column(db.Float)
    probetimeunits = db.Column(db.String)

    concentration_two = db.Column(db.Float)
    concentration_three = db.Column(db.Float)
    organism = db.Column(db.String)
    project = db.Column(db.String)
    data_acq = db.Column(db.String)

    comments = db.Column(db.String)
    isreference = db.Column(db.Boolean, index=True)

    plex = relationship("Plex", back_populates="compoundtreatments")
    celltype = relationship("CellType", back_populates="compoundtreatments")

    intensityreadings = relationship("IntensityReading", overlaps = "experiment", back_populates = "compoundtreatment")
    competitionratios = relationship("CompetitionRatio", overlaps = "experiment", back_populates = "compoundtreatment")

    compound = relationship(
        "Compound",
        foreign_keys=[compound_id],
        back_populates="compoundtreatments"
    )

    probe = relationship(
        "Compound",
        foreign_keys=[probe_id]
    )





# CHEMOPROTEOMICS SCHEMA

## RESIDUE LEVEL

class Residue(db.Model):
    __tablename__ = "residues"
    __table_args__ = {"schema": "chemoproteomics"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    uniprot = db.Column(db.String, ForeignKey("core.proteins.uniprot"), index=True)
    position = db.Column(db.Integer)
    type = db.Column(db.String(20), index=True)

    UniqueConstraint("uniprot", "position")

    protein = relationship("Protein", back_populates="residues")
    structureresidues = relationship("StructureResidue", back_populates = "residue")

    intensityreadings = relationship("IntensityReading", back_populates="residue")
    competitionratios = relationship("CompetitionRatio", back_populates="residue")
    lists = relationship(
        "ResidueToList",
        back_populates="residue"
    )

class ResidueList(db.Model):
    __tablename__ = 'residuelists'
    __table_args__ = {"schema": "chemoproteomics"}

    id = db.Column(UUID(as_uuid=True), server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    name = db.Column(db.String())
    description = db.Column(db.String())

    residues = relationship(
    "ResidueToList",
    back_populates="residuelist")

class ResidueToList(db.Model):
    __tablename__ = 'list2residue'
    __table_args__ = {"schema": "chemoproteomics"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True, unique=True)

    residue_id = db.Column(
        UUID(as_uuid=True),
        ForeignKey("chemoproteomics.residues.id"),
        index=True
    )

    residuelist_id = db.Column(
        UUID(as_uuid=True),
        ForeignKey("chemoproteomics.residuelists.id")
    )
    residue = relationship("Residue", back_populates="lists")

    residuelist = relationship("ResidueList", back_populates="residues")

class ResidueFeature(db.Model):
     __tablename__ = 'residuefeatures'
     __table_args__ = {"schema": "chemoproteomics"}

     id = db.Column(UUID(as_uuid=True), server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
     residue_id = db.Column(UUID(as_uuid=True), ForeignKey("chemoproteomics.residues.id"), index = True)
     description = db.Column(db.String(), index = True)
     source = db.Column(db.String(), index = True)

     UniqueConstraint(residue_id, description, source)
     
class Compound_cr_four(db.Model):
    __tablename__ = 'compound_cr_four'
    __table_args__ = {"schema": "chemoproteomics"}

    compound_id = db.Column(db.String, primary_key=True)
    count = db.Column(db.Integer)

class Compound_cr_fifteen(db.Model):
    __tablename__ = 'compound_cr_fifteen'
    __table_args__ = {"schema": "chemoproteomics"}

    compound_id = db.Column(db.String, primary_key=True)
    count = db.Column(db.Integer)

##INTENSITIES AND CR

class IntensityReading(db.Model):
    __tablename__ = "intensityreadings"
    __table_args__ = {"schema": "chemoproteomics"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    plex_id = db.Column(db.String, ForeignKey("core.plexs.id"), index=True)
    compoundtreatment_id = db.Column(UUID(as_uuid=True), ForeignKey("core.compoundtreatments.id"), index=True)
    residue_id = db.Column(UUID(as_uuid=True), ForeignKey("chemoproteomics.residues.id"), index=True)

    scan = db.Column(db.String, index=True)
    value = db.Column(db.Float)
    peptideseq = db.Column(db.String)
    modification = db.Column(db.String)

    multimapper = db.Column(db.Boolean)

    UniqueConstraint("plex_id", "compoundtreatment_id", "scan", "residue_id")

    plex = relationship("Plex", back_populates = "intensityreadings")
    compoundtreatment = relationship("CompoundTreatment", back_populates = "intensityreadings")
    residue = relationship("Residue", back_populates = "intensityreadings")


class CompetitionRatio(db.Model):
    __tablename__ = "competitionratios"
    __table_args__ = {"schema": "chemoproteomics"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    plex_id = db.Column(db.String, ForeignKey("core.plexs.id"), index=True)
    compoundtreatment_id = db.Column(UUID(as_uuid=True), ForeignKey("core.compoundtreatments.id"), index=True)
    residue_id = db.Column(UUID(as_uuid=True), ForeignKey("chemoproteomics.residues.id"), index=True)

    scan = db.Column(db.String, index=True)
    cr = db.Column(db.Float)

    control_rsd = db.Column(db.Float(), index = True) #Relative SD of the control samples (DMSO)
    display_flag = db.Column(db.Boolean(), default=False, nullable=False)
    multimapper = db.Column(db.Boolean(), default=False, nullable=False)
    group_id = db.Column(UUID(as_uuid=True), default=uuid.uuid4)
    group_cr = db.Column(db.Float(), index = True)
    p_value = db.Column(db.Float(), index= True)
    replicate_no = db.Column(db.Float(), index = True)

    UniqueConstraint("plex_id", "compoundtreatment_id", "scan")

    plex = relationship("Plex", back_populates = "competitionratios")
    compoundtreatment = relationship("CompoundTreatment", back_populates = "competitionratios")
    residue = relationship("Residue", back_populates="competitionratios")


#PROTEOMICS SCHEMA

class FoldChange(db.Model):
    __tablename__ = "foldchanges"
    __table_args__ = {"schema": "proteomics"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    plex_id = db.Column(db.String, ForeignKey("core.plexs.id"), index=True)
    compoundtreatment_id = db.Column(UUID(as_uuid=True), ForeignKey("core.compoundtreatments.id"), index=True)
    protein_id = db.Column(db.String, ForeignKey("core.proteins.uniprot"), index=True)

    scan = db.Column(db.String)
    foldchange = db.Column(db.Float)
    p_value = db.Column(db.Float)
    replicate_no = db.Column(db.Float)
    group_id = db.Column(UUID(as_uuid=True), index=True)
    no_peptides = db.Column(db.Float)
    no_unique_peptides = db.Column(db.Float)

    UniqueConstraint("plex_id", "compoundtreatment_id", "scan")

    protein = relationship("Protein")
    compoundtreatment = relationship("CompoundTreatment")
    plex = relationship("Plex")


class ProteinIntensityReading(db.Model):
    __tablename__ = "protein_intensityreadings"
    __table_args__ = {"schema": "proteomics"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    plex_id = db.Column(db.String, ForeignKey("core.plexs.id"), index=True)
    compoundtreatment_id = db.Column(UUID(as_uuid=True), ForeignKey("core.compoundtreatments.id"), index=True)
    protein_id = db.Column(db.String, ForeignKey("core.proteins.uniprot"), index=True)

    scan = db.Column(db.String)
    value = db.Column(db.Float)


# CORE SCHEMA - LISTS AND PATHWAYS

class TargetList(db.Model):
    __tablename__ = 'targetlists'
    __table_args__ = {"schema": "core"}

    id = db.Column(UUID(as_uuid=True), server_default=text("uuid_generate_v4()"), primary_key=True, unique=True)
    name = db.Column(db.String())
    description = db.Column(db.String())

    proteins = relationship("Protein", secondary="core.list2protein", back_populates="lists", uselist=True, overlaps="lists")


class ProteinToList(db.Model):
    __tablename__ = 'list2protein'
    __table_args__ = {"schema": "core"}

    id = db.Column(UUID(as_uuid=True), server_default=text("uuid_generate_v4()"), primary_key=True, unique=True)
    protein_id = db.Column(db.String, ForeignKey("core.proteins.uniprot"))
    targetlist_id = db.Column(UUID(as_uuid=True), ForeignKey("core.targetlists.id"))


class PathwayList(db.Model):
    __tablename__ = 'pathwaylists'
    __table_args__ = {"schema": "core"}

    id = db.Column(UUID(as_uuid=True), server_default=text("uuid_generate_v4()"), primary_key=True, unique=True)
    name = db.Column(db.String())
    description = db.Column(db.String())
    database = db.Column(db.String())
    organism = db.Column(db.String())

    proteins = relationship("Protein", secondary="core.protein2pathway", back_populates="pathways", uselist=True, overlaps="pathways")


class ProteinToPathway(db.Model):
    __tablename__ = 'protein2pathway'
    __table_args__ = {"schema": "core"}

    id = db.Column(UUID(as_uuid=True), server_default=text("uuid_generate_v4()"), primary_key=True, unique=True)
    protein_id = db.Column(db.String, ForeignKey("core.proteins.uniprot"))
    pathway_id = db.Column(UUID(as_uuid=True), ForeignKey("core.pathwaylists.id"))


#STRUCUTRE SCHEMA


class Structure(db.Model):
    __tablename__ = "structures"
    __table_args__ = {"schema": "structure"}

    id = db.Column(db.String(10), primary_key=True)
    type = db.Column(db.String)
    resolution = db.Column(db.String)

    chains = relationship("StructureChain", back_populates = "structure")
    ligands = relationship("Ligand", back_populates = "structure")
    pockets = relationship("Pocket", back_populates = "structure")


class StructureChain(db.Model):
    __tablename__ = "structurechains"
    __table_args__ = {"schema": "structure"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    structure_id = db.Column(db.String, ForeignKey("structure.structures.id"))
    chain = db.Column(db.String(), index = True)
    uniprot_id = db.Column(db.String, ForeignKey("core.proteins.uniprot"))

    UniqueConstraint("structure_id", "chain")

    structure = relationship("Structure", back_populates = "chains")
    residues = relationship("StructureResidue", back_populates = "chain")
    protein = relationship("Protein", back_populates = "chains")

class StructureResidue(db.Model):
    __tablename__ = "structureresidues"
    __table_args__ = {"schema": "structure"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    structure_id = db.Column(db.String, ForeignKey("structure.structures.id"))
    chain_id = db.Column(UUID(as_uuid=True), ForeignKey("structure.structurechains.id"))
    residue_id = db.Column(UUID(as_uuid=True), ForeignKey("chemoproteomics.residues.id"))

    pdb_position = db.Column(db.String(), index = True) # Position from PDB file
    confidence = db.Column(db.Float()) #AF2 confidence
    in_disulfide = db.Column(db.Boolean(), index = True)
    accessibility = db.Column(db.Float(), index = True)
    depth = db.Column(db.Float(), index = True)

    position = association_proxy('residue', 'position')
    type = association_proxy('residue', 'type')
    uniprot = association_proxy('chain', 'uniprot_id')

    UniqueConstraint('structure_id', 'chain_id', 'residue_id')

    chain = relationship("StructureChain", back_populates = "residues")
    residue = relationship("Residue", back_populates = "structureresidues")

    liganddistances = relationship("LigandResidueDistance",back_populates = "residue",uselist = True)
    pockets = relationship("Pocket",secondary = 'structure.pocket2residue',back_populates = "residues",uselist = True)



class Ligand(db.Model):
    __tablename__ = "ligands"
    __table_args__ = {"schema": "structure"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    structure_id = db.Column(db.String, ForeignKey("structure.structures.id"), index=True)

    code = db.Column(db.String, index=True)   # PDB 3-letter code
    chain = db.Column(db.String, index=True)

    smiles = db.Column(db.String, index=True)
    inchi = db.Column(db.String, index=True)
    chembl = db.Column(db.String, index=True)
    name = db.Column(db.String, index=True)

    mw = db.Column(db.Float)
    artefact = db.Column(db.Boolean)

    UniqueConstraint('structure_id', 'code', 'chain')

    structure = relationship("Structure", back_populates = "ligands")
    residuedistances = relationship("LigandResidueDistance", back_populates = "ligand")

class LigandResidueDistance(db.Model):
    __tablename__ = "residue2ligand"
    __table_args__ = {"schema": "structure"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    ligand_id = db.Column(UUID(as_uuid=True), ForeignKey("structure.ligands.id"))
    structureresidue_id = db.Column(UUID(as_uuid=True), ForeignKey("structure.structureresidues.id"))

    distance = db.Column(db.Float, index=True)

    UniqueConstraint(ligand_id, structureresidue_id)

    ligand = relationship("Ligand", back_populates = "residuedistances")
    residue = relationship("StructureResidue", back_populates = "liganddistances")


class PocketResidue(db.Model):
    __tablename__ = "pocket2residue"
    __table_args__ = {"schema": "structure"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    pocket_id = db.Column(UUID(as_uuid=True), ForeignKey("structure.pockets.id"), index=True)
    structureresidue_id = db.Column(UUID(as_uuid=True), ForeignKey("structure.structureresidues.id"), index=True)

    UniqueConstraint("pocket_id", "structureresidue_id")


class Pocket(db.Model):
    __tablename__ = "pockets"
    __table_args__ = {"schema": "structure"}

    id = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)

    structure_id = db.Column(db.String, ForeignKey("structure.structures.id"), index=True)

    pocket_id = db.Column(db.String)  # fpocket / AF2 derived id

    pocket_score = db.Column(db.Float, index=True)
    drug_score = db.Column(db.Float, index=True)

    mean_confidence = db.Column(db.Float)
    median_confidence = db.Column(db.Float)
    min_confidence = db.Column(db.Float)

    pocket_volume_MC = db.Column(db.Float)
    pocket_volume_hull = db.Column(db.Float)

    UniqueConstraint("structure_id", "pocket_id")

    structure = relationship("Structure", back_populates = "pockets")
    residues = relationship("StructureResidue", secondary = 'structure.pocket2residue', back_populates = "pockets",uselist = True, overlaps = "pockets")


class ProteinEmbedding(db.Model):
    __tablename__ = "protein_embeddings"
    __table_args__ = {"schema": "core"}

    id       = db.Column(UUID(as_uuid=True), default=uuid.uuid4, primary_key=True)
    uniprot  = db.Column(db.String, ForeignKey("core.proteins.uniprot"), nullable=False, index=True)
    chunk_text = db.Column(db.Text, nullable=False)
    embedding  = db.Column(Vector(768), nullable=False)

    protein = relationship("Protein", backref="embeddings")
