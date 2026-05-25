from enum import unique
from operator import index
from flask_sqlalchemy import SQLAlchemy

from sqlalchemy import ForeignKey, Integer, String, UniqueConstraint
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.orm import relationship

from sqlalchemy.dialects.postgresql import UUID
import uuid
from sqlalchemy_utils import UUIDType
from sqlalchemy.sql import text

from extensions import db, migrate


class Residue(db.Model):
    __tablename__ = 'residues'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    uniprot = db.Column(db.String(20), ForeignKey("proteins.uniprot"), index = True)
    position = db.Column(db.Integer(), primary_key = True) #UniProt position
    type = db.Column(db.String(20), index = True) #Three character AA code - allows longer in case of edges cases

    UniqueConstraint('uniprot', 'position')

    protein = relationship("Protein", back_populates = "residues") # Link to proteins
    structureresidues = relationship("StructureResidue", back_populates = "residue")
    intensityreadings = relationship("IntensityReading", back_populates = "residue")
    competitionratios = relationship("CompetitionRatio", back_populates = "residue")
    lists = relationship("ResidueList", secondary = "list2residue", back_populates = "residues", uselist = True, overlaps = "residues")

    def __repr__(self):
        return f"{self.uniprot}:{self.type}{self.position}"

class Protein(db.Model):
    __tablename__ = 'proteins'

    uniprot = db.Column(db.String(), primary_key = True)
    description = db.Column(db.String())
    symbol = db.Column(db.String(), index = True)

    residues = relationship("Residue", back_populates = "protein", order_by = "Residue.position")
    chains = relationship("StructureChain", back_populates = "protein")
    synonyms = relationship("ProteinSynonym", back_populates = "protein")
    lists = relationship("TargetList", secondary = "list2protein", back_populates = "proteins", uselist = True, overlaps = "proteins")

    def __repr__(self):
        return f"{self.uniprot}: {self.symbol}"

class TargetList(db.Model):
    __tablename__ = 'targetlists'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    name = db.Column(db.String())
    description = db.Column(db.String())

    proteins = relationship("Protein", secondary = 'list2protein', back_populates = "lists", uselist = True, overlaps = "lists")

class ProteinToList(db.Model):
    __tablename__ = 'list2protein'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    protein_id = db.Column(db.String, ForeignKey("proteins.uniprot"))
    targetlist_id = db.Column(UUIDType, ForeignKey("targetlists.id"))

class ProteinSynonym(db.Model):
    __tablename__ = 'proteinsynonyms'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    uniprot = db.Column(db.String(), ForeignKey("proteins.uniprot"))
    synonym = db.Column(db.String()) #Includes protein 'names' as well as ENSG ids
    type = db.Column(db.String()) #Should be an enum really, but specify whether Ensembl or name
    
    protein = relationship("Protein", back_populates = "synonyms") 

class CompoundTreatment(db.Model):
    __tablename__ = 'compoundtreatments'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)    
    plex_id = db.Column(db.String(),  ForeignKey("plexs.id"))
    samplename = db.Column(db.String())
    compound_id = db.Column(db.String(),  ForeignKey("compounds.id"), index = True)
    celltype_id = db.Column(UUIDType, ForeignKey("celltypes.id"))
    tmtchannel = db.Column(db.String())
    concentration = db.Column(db.Float())
    concentrationunits = db.Column(db.String())
    time = db.Column(db.Float())
    timeunits = db.Column(db.String())
    temperature = db.Column(db.Float()) #To add
    probe_id = db.Column(db.String(), ForeignKey("compounds.id"), index = True)
    probeconcentration = db.Column(db.Float())
    probeconcentrationunits = db.Column(db.String())
    probetime = db.Column(db.Float())
    probetimeunits = db.Column(db.String())
    comments = db.Column(db.String())
    isreference = db.Column(db.Boolean(), index = True)

    UniqueConstraint(plex_id, samplename) #For a given plex the samplename must be unique and the plex/sample combination defines the treatment

    plex = relationship("Plex", back_populates = "compoundtreatments")
    celltype = relationship("CellType",   back_populates = "compoundtreatments")

    intensityreadings = relationship("IntensityReading", overlaps = "experiment", back_populates = "compoundtreatment")
    competitionratios = relationship("CompetitionRatio", overlaps = "experiment", back_populates = "compoundtreatment")

    compound = relationship("Compound", foreign_keys = [compound_id]) #both compound and probe are compounds so need separate foreign keys defined
    probe = relationship("Compound", foreign_keys = [probe_id])

class Plex(db.Model):
    __tablename__ = 'plexs'

    id = db.Column(db.String(), primary_key = True)
    experiment_id = db.Column(db.String(), ForeignKey("experiments.id"))
    description = db.Column(db.String())

    experiment = relationship("Experiment", back_populates = "plexs")
    compoundtreatments = relationship("CompoundTreatment", back_populates = "plex", overlaps = "intensityreadings")
    intensityreadings = relationship("IntensityReading",  back_populates = "plex")
    competitionratios = relationship("CompetitionRatio",  back_populates = "plex")

class Experiment(db.Model):
    __tablename__ = 'experiments'

    id = db.Column(db.String(), primary_key = True)
    description = db.Column(db.String())

    plexs = relationship("Plex", back_populates = "experiment")
       

class CellType(db.Model):
    __tablename__ = 'celltypes'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    name = db.Column(db.String())
    depmapid = db.Column(db.String())
    description = db.Column(db.String())

    compoundtreatments = relationship("CompoundTreatment", back_populates = "celltype")

class Compound(db.Model):
    __tablename__ = 'compounds'

    id = db.Column(db.String(), primary_key = True)
    smiles = db.Column(db.String(), index = True)
    image = db.Column(db.String(), index = True)

    compoundtreatments = relationship(
        "CompoundTreatment", 
        back_populates = "compound", 
        primaryjoin = "and_(Compound.id ==CompoundTreatment.compound_id)" # Needs defining because the CompoundTreatment table also has a key for the probe
    )

class IntensityReading(db.Model):
    __tablename__ = 'intensityreadings'
    
    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    plex_id = db.Column(db.String(),  ForeignKey("plexs.id"), index = True)
    compoundtreatment_id = db.Column(UUIDType, ForeignKey("compoundtreatments.id"), index = True)
    residue_id = db.Column(UUIDType, ForeignKey("residues.id"), index = True)
    scan = db.Column(db.String(), index = True)
    value = db.Column(db.Float())
    peptideseq = db.Column(db.String(), index = True)
    modification = db.Column(db.String(), index = True)
    multimapper = db.Column(db.Boolean())
    
    UniqueConstraint('plex_id', 'compoundtreatment_id', 'scan', 'residue_id') #The combination of plex, compound treatment and scan is unique until we consider multimappers
    #This makes the residue and additional unique constraint. 

    plex = relationship("Plex", back_populates = "intensityreadings")
    compoundtreatment = relationship("CompoundTreatment", back_populates = "intensityreadings")
    residue = relationship("Residue", back_populates = "intensityreadings")

class CompetitionRatio(db.Model):
    __tablename__ = 'competitionratios'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    plex_id = db.Column(db.String(),  ForeignKey("plexs.id"), index = True)
    compoundtreatment_id = db.Column(UUIDType, ForeignKey("compoundtreatments.id"), index = True)
    residue_id = db.Column(UUIDType, ForeignKey("residues.id"), index = True)
    scan = db.Column(db.String(), index = True)
    cr = db.Column(db.Float(), index = True)
    control_rsd = db.Column(db.Float(), index = True) #Relative SD of the control samples (DMSO)
    display_flag = db.Column(db.Boolean(), default=False, nullable=False)
    multimapper = db.Column(db.Boolean(), default=False, nullable=False)
    group_id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), index = True)
    group_cr = db.Column(db.Float(), index = True)
    p_value = db.Column(db.Float(), index= True)
    replicate_no = db.Column(db.Float(), index = True)

    UniqueConstraint('plex_id', 'compoundtreatment_id', 'scan')

    plex = relationship("Plex", back_populates = "competitionratios")
    compoundtreatment = relationship("CompoundTreatment", back_populates = "competitionratios")
    residue = relationship("Residue", back_populates = "competitionratios")

class Structure(db.Model):
    __tablename__ = 'structures'

    id = db.Column(db.String(10), primary_key = True) #PDB code or uniprot ID (for AF2)
    type = db.Column(db.String(10))
    resolution = db.Column(db.String(10))

    chains = relationship("StructureChain", back_populates = "structure")
    ligands = relationship("Ligand", back_populates = "structure")
    pockets = relationship("Pocket", back_populates = "structure")

class StructureChain(db.Model):
    __tablename__ = 'structurechains'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    structure_id = db.Column(db.String(), ForeignKey("structures.id"))
    chain = db.Column(db.String(), index = True) #PDB chain code (may be blank)
    uniprot_id = db.Column(db.String(), ForeignKey("proteins.uniprot"), index = True)

    UniqueConstraint('structure_id', 'chain')

    structure = relationship("Structure", back_populates = "chains")
    residues = relationship("StructureResidue", back_populates = "chain")
    protein = relationship("Protein", back_populates = "chains")

class StructureResidue(db.Model):
    __tablename__ = 'structureresidues'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    structure_id = db.Column(db.String(), ForeignKey("structures.id"), index = True)
    chain_id = db.Column(UUIDType, ForeignKey("structurechains.id"), index = True) #Auto id not the PDB chain ID
    residue_id = db.Column(UUIDType, ForeignKey("residues.id"), index = True)

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
    pockets = relationship("Pocket",secondary = 'pocket2residue',back_populates = "residues",uselist = True)

class Ligand(db.Model):
    __tablename__ = 'ligands'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    structure_id = db.Column(db.String(), ForeignKey("structures.id"), index = True)
    code = db.Column(db.String(), index = True) #Three letter code in PDB
    chain = db.Column(db.String(), index = True)
    smiles = db.Column(db.String(), index = True)
    inchi = db.Column(db.String(), index = True)
    chembl = db.Column(db.String(), index = True)
    name = db.Column(db.String(), index = True)
    mw = db.Column(db.Float(), index = True)
    artefact = db.Column(db.Boolean())

    UniqueConstraint('structure_id', 'code', 'chain')

    structure = relationship("Structure", back_populates = "ligands")
    residuedistances = relationship("LigandResidueDistance", back_populates = "ligand")

class LigandResidueDistance(db.Model):
    __tablename__ = "residue2ligand"

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    ligand_id = db.Column(UUIDType, ForeignKey("ligands.id"), index = True)
    structureresidue_id = db.Column(UUIDType, ForeignKey("structureresidues.id"), index = True)
    distance = db.Column(db.Float(), index = True)

    UniqueConstraint(ligand_id, structureresidue_id)

    ligand = relationship("Ligand", back_populates = "residuedistances")
    residue = relationship("StructureResidue", back_populates = "liganddistances")

class PocketResidue(db.Model):
    __tablename__ = "pocket2residue"

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    pocket_id = db.Column(UUIDType, ForeignKey('pockets.id'), primary_key = True, index = True)
    structureresidue_id = db.Column(UUIDType, ForeignKey('structureresidues.id'), primary_key = True, index = True)

    UniqueConstraint(pocket_id, structureresidue_id)

    #TODO: Perhaps we want this for every pocket nucleophile as a distance as we did for ligands?
    #For now any CYS atom touching a pocket counts (not just SG atom)

class Pocket(db.Model):
    __tablename__ = "pockets"

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    structure_id = db.Column(db.String(), ForeignKey("structures.id"))
    pocket_id = db.Column(db.String()) #Make a string for future proof though fpocket gives numeric ids
    pocket_score = db.Column(db.Float(), index = True)
    drug_score = db.Column(db.Float(), index = True)
    mean_confidence = db.Column(db.Float(), index = True) # Confidende from AF2 data. Set to null or 100 for PDB pockets
    median_confidence = db.Column(db.Float(), index = True)
    min_confidence = db.Column(db.Float(), index = True)
    pocket_volume_MC = db.Column(db.Float(), index = True)
    pocket_volume_hull = db.Column(db.Float(), index = True)

    UniqueConstraint(structure_id, pocket_id)

    structure = relationship("Structure", back_populates = "pockets")
    residues = relationship("StructureResidue", secondary = 'pocket2residue', back_populates = "pockets",uselist = True, overlaps = "pockets")

class ResidueList(db.Model):
    __tablename__ = 'residuelists'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    name = db.Column(db.String())
    description = db.Column(db.String())

    residues = relationship("Residue", secondary = 'list2residue', back_populates = "lists", uselist = True, overlaps = "lists")

class ResidueToList(db.Model):
    __tablename__ = 'list2residue'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    residue_id = db.Column(UUIDType, ForeignKey("residues.id"), index = True)
    residuelist_id = db.Column(UUIDType, ForeignKey("residuelists.id"))
    
class ResidueFeature(db.Model):
     __tablename__ = 'residuefeatures'
     id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
     residue_id = db.Column(UUIDType, ForeignKey("residues.id"), index = True)
     description = db.Column(db.String(), index = True)
     source = db.Column(db.String(), index = True)

     UniqueConstraint(residue_id, description, source)
     
class Compound_cr_four(db.Model):
    __tablename__ = 'compound_cr_four'
    compound_id = db.Column(db.String, primary_key=True)
    count = db.Column(db.Integer)

class Compound_cr_fifteen(db.Model):
    __tablename__ = 'compound_cr_fifteen'
    compound_id = db.Column(db.String, primary_key=True)
    count = db.Column(db.Integer)

