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

from extensions import db, migrate



class GpProtein(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpproteins'
    
    uniprot = db.Column(db.String(), primary_key = True)
    description = db.Column(db.String())
    symbol = db.Column(db.String(), index = True)
    organism = db.Column(db.String())

    gpsynonyms = relationship("GpProteinSynonym", back_populates = "gpprotein")
    foldchanges = relationship("FoldChange",  back_populates = "gpprotein")
    gpintensityreadings = relationship("GpIntensityReading",  back_populates = "gpprotein")
    pathways = relationship("PathwayList", secondary = "protein2pathway", back_populates = "gpprotein", uselist = True, overlaps = "gpprotein")

    def __repr__(self):
        return f"{self.uniprot}: {self.symbol}"
    

class GpProteinSynonym(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpproteinsynonyms'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    uniprot = db.Column(db.String(), ForeignKey("gpproteins.uniprot"))
    synonym = db.Column(db.String()) #Includes protein 'names' as well as ENSG ids
    type = db.Column(db.String()) #Should be an enum really, but specify whether Ensembl or name
    
    gpprotein = relationship("GpProtein", back_populates = "gpsynonyms") 

class GpCompoundTreatment(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpcompoundtreatments'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)    
    plex_id = db.Column(db.String(),  ForeignKey("gpplexs.id"))
    samplename = db.Column(db.String())
    compound_id = db.Column(db.String(),  ForeignKey("gpcompounds.id"), index = True)
    celltype_id = db.Column(UUIDType, ForeignKey("gpcelltypes.id"))
    tmtchannel = db.Column(db.String())
    concentration = db.Column(db.Float())
    concentration_two = db.Column(db.Float())
    concentration_three = db.Column(db.Float())
    concentrationunits = db.Column(db.String())
    data_acq = db.Column(db.String())
    time = db.Column(db.Float())
    timeunits = db.Column(db.String())
    temperature = db.Column(db.String()) #To add
    comments = db.Column(db.String())
    isreference = db.Column(db.Boolean(), index = True)
    project = db.Column(db.String())
    organism = db.Column(db.String())
    UniqueConstraint(plex_id, samplename) #For a given plex the samplename must be unique and the plex/sample combination defines the treatment

    gpplex = relationship("GpPlex", back_populates = "gpcompoundtreatments")
    gpcelltype = relationship("GpCellType",   back_populates = "gpcompoundtreatments")

    gpintensityreadings = relationship("GpIntensityReading", overlaps = "gpexperiment", back_populates = "gpcompoundtreatments")
    foldchanges = relationship("FoldChange", overlaps = "gpexperiment", back_populates = "gpcompoundtreatments")

    gpcompound = relationship("GpCompound", foreign_keys = [compound_id]) #both compound and probe are compounds so need separate foreign keys defined

class GpPlex(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpplexs'

    id = db.Column(db.String(), primary_key = True)
    experiment_id = db.Column(db.String(), ForeignKey("gpexperiments.id"))
    description = db.Column(db.String())

    gpexperiment = relationship("GpExperiment", back_populates = "gpplexs")
    gpcompoundtreatments = relationship("GpCompoundTreatment", back_populates = "gpplex", overlaps = "gpintensityreadings")
    gpintensityreadings = relationship("GpIntensityReading",  back_populates = "gpplex")
    foldchanges = relationship("FoldChange",  back_populates = "gpplex")

class GpExperiment(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpexperiments'

    id = db.Column(db.String(), primary_key = True)
    description = db.Column(db.String())

    gpplexs = relationship("GpPlex", back_populates = "gpexperiment")
       

class GpCellType(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpcelltypes'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    name = db.Column(db.String())
    depmapid = db.Column(db.String())
    description = db.Column(db.String())

    gpcompoundtreatments = relationship("GpCompoundTreatment", back_populates = "gpcelltype")

class GpCompound(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpcompounds'

    id = db.Column(db.String(), primary_key = True)
    smiles = db.Column(db.String(), index = True)
    image = db.Column(db.String(), index = True)

    gpcompoundtreatments = relationship(
        "GpCompoundTreatment", 
        back_populates = "gpcompound", 
        primaryjoin = "and_(GpCompound.id ==GpCompoundTreatment.compound_id)" # Needs defining because the CompoundTreatment table also has a key for the probe
    )

class GpIntensityReading(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'gpintensityreadings'
    
    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    plex_id = db.Column(db.String(),  ForeignKey("gpplexs.id"), index = True)
    compoundtreatment_id = db.Column(UUIDType, ForeignKey("gpcompoundtreatments.id"), index = True)
    protein_id = db.Column(db.String, ForeignKey("gpproteins.uniprot"), index = True)
    scan = db.Column(db.String(), index = True)
    value = db.Column(db.Float())
    
    UniqueConstraint('plex_id', 'compoundtreatment_id', 'scan', 'protein_id') #The combination of plex, compound treatment and scan is unique until we consider multimappers
    #This makes the residue and additional unique constraint. 

    gpplex = relationship("GpPlex", back_populates = "gpintensityreadings")
    gpcompoundtreatments = relationship("GpCompoundTreatment", back_populates = "gpintensityreadings")
    gpprotein = relationship("GpProtein", back_populates = "gpintensityreadings")
    
class FoldChange(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'foldchanges'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    plex_id = db.Column(db.String(),  ForeignKey("gpplexs.id"), index = True)
    compoundtreatment_id = db.Column(UUIDType, ForeignKey("gpcompoundtreatments.id"), index = True)
    protein_id = db.Column(db.String, ForeignKey("gpproteins.uniprot"), index = True)
    scan = db.Column(db.String(), index = True)
    foldchange = db.Column(db.Float(), index = True)
    p_value = db.Column(db.Float(), index= True)
    replicate_no = db.Column(db.Float(), index = True)
    control_rsd = db.Column(db.Float(), index = True) 
    no_peptides = db.Column(db.Float(), index = True)
    no_unique_peptides = db.Column(db.Float(), index = True)
    group_id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), index = True)
    UniqueConstraint('plex_id', 'compoundtreatment_id', 'scan')

    gpplex = relationship("GpPlex", back_populates = "foldchanges")
    gpcompoundtreatments = relationship("GpCompoundTreatment", back_populates = "foldchanges")
    gpprotein = relationship("GpProtein", back_populates = "foldchanges")


class PathwayList(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'pathwaylists'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    name = db.Column(db.String())
    description = db.Column(db.String())
    database = db.Column(db.String())
    organism = db.Column(db.String())

    gpprotein = relationship("GpProtein", secondary = 'protein2pathway', back_populates = "pathways", uselist = True, overlaps = "pathways")



class ProteinToPathway(db.Model):
    __bind_key__ = 'globalproteomics'
    __tablename__ = 'protein2pathway'

    id = db.Column(UUIDType, server_default = text("uuid_generate_v4()"), primary_key = True, unique = True)
    protein_id = db.Column(db.String, ForeignKey("gpproteins.uniprot"))
    pathway_id = db.Column(UUIDType, ForeignKey("pathwaylists.id"))
