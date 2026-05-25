import click
import glob
import os
import re
import pandas as pd
from flask import Blueprint
from os.path import exists
from tqdm import tqdm
import logging

from sqlalchemy import or_, func

import gzip
import xml.etree.ElementTree as ET
import uuid

from rdkit import Chem
from rdkit.Chem import Draw
import statistics as stat
from app.models import CompetitionRatio, Protein, Residue, CompoundTreatment, Compound, CellType, Experiment, IntensityReading

from scipy.stats import ttest_ind


databp = Blueprint("data", __name__)

@databp.cli.command("experiment_import")
@click.option('-d', '--datafile', required=True)
@click.option('-s', '--sampletable', required=True)
@click.option('-l', '--logfile', default='discoverome_data_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def experiment_import(datafile, sampletable, logfile, drop_all):
    """Imports an SLC-ABPP data file and sample sheet"""

    from models.chemopro import db, Protein, Residue, CompoundTreatment, Compound, CellType, Experiment, IntensityReading, Plex

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    logging.info(f"Importing {datafile} and {sampletable}")
    
    if not exists(datafile):
        raise FileNotFoundError

    if not exists(sampletable):
        raise FileNotFoundError

    if drop_all:
        print("***Dropping and creating all tables***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)
   
    #Read the sample table and import each sample as a compound treatment object with associated compounds and probes
    #Now find all the compounds (and probes) and cell types referenced in the file
    sample_table = pd.read_csv(sampletable)
    compounds = set()
    celltypes = set()
    samplenames = set()
    experiment_id = set()
    plex_id = set()
    for _, row in sample_table.iterrows():
        compounds.add(row['Compound'])
        compounds.add(row['Probe'])
        celltypes.add(row['CellType'])
        samplenames.add(row['SampleName'])
        experiment_id.add(row['Experiment'])
        plex_id.add(row['Plex'])

    #There should only be one experiment and one plex per datafile
    if len(experiment_id) > 1 or len(plex_id) > 1:
        logging.error(f"More than one experiment or plex ID in {sampletable}")
        raise Exception(f"More than one experiment or plex ID in {sampletable}")
    else:
        experiment_id = experiment_id.pop()
        plex_id = plex_id.pop()

    #Find the experiment or create it if it is new
    #If the plex already exist then throw an error
    experiment = db.session.query(Experiment).filter(Experiment.id == experiment_id)
    if not db.session.query(experiment.exists()).scalar():
        experiment = Experiment(id = experiment_id, description = "Experiment description")
        db.session.add(experiment)
        db.session.commit()
    else:
        experiment = experiment.first()

    plex = db.session.query(Plex).filter(Plex.id == plex_id)
    if not db.session.query(plex.exists()).scalar():
        plex = Plex(id = plex_id, experiment = experiment, description = "Plex description")
        db.session.add(plex)
        db.session.commit()
    else:
        logging.error(f"Found plex: {plex_id} in db already.")
        raise Exception(f"Found plex: {plex_id} in db already.")

    #Create compound objects for each one as long as it doesn't exist in the db
    for compound_id in compounds:
        logging.info(f"Searching for {compound_id}")
        query = db.session.query(Compound).filter(Compound.id == compound_id)
        if not db.session.query(query.exists()).scalar():
            logging.info(f"Adding compound {compound_id}")
            cpd = Compound(id = compound_id)
            db.session.add(cpd)
    db.session.commit()

    for celltype in celltypes:
        logging.info(f"Searching for cell type {celltype}")
        query = db.session.query(CellType).filter(CellType.name == celltype)
        if not db.session.query(query.exists()).scalar():
            logging.info(f"Adding cell type {celltype}")
            clt = CellType(name = celltype)
            db.session.add(clt)
    db.session.commit()

    #Now create the compound treatment objects (also store these in a dictionary for quick look up later)
    cpdtrts = {}
    cpdtrt_groups = {}
    refs = []
    for _, row, in sample_table.iterrows():
        logging.info(f"Adding compound treatment: {row['SampleName']}")
        compound = db.session.query(Compound).filter(Compound.id == row['Compound']).first()
        probe = db.session.query(Compound).filter(Compound.id == row['Probe']).first()
        celltype = db.session.query(CellType).filter(CellType.name == row['CellType']).first()
        cpdtrt = CompoundTreatment(
            plex = plex,
            samplename = row['SampleName'],
            compound = compound,
            probe = probe,
            celltype = celltype,
            tmtchannel = row['TMTChannel'],
            concentration = row['Concentration'],
            concentrationunits = row['ConcentrationUnits'],
            time = row['Time'],
            timeunits = row['TimeUnits'],
            probeconcentration = row['ProbeConcentration'],
            probeconcentrationunits = row['ProbeConcentrationUnits'],
            probetime = row['ProbeTime'],
            probetimeunits = row['ProbeTimeUnits'],
            comments = row['Comments'],
            isreference = row['IsReference'] == "Y"
        )
        try:
            db.session.add(cpdtrt)
            db.session.flush()
            db.session.refresh(cpdtrt)
            cpdtrts[row['SampleName']] = cpdtrt.id
            cpdtrt_groups[row['SampleName']] = \
                row['Compound']+str(row['Concentration'])+row['ConcentrationUnits']+ \
                str(row['Time'])+row['TimeUnits']+\
                row['Probe']+str(row['ProbeConcentration'])+row['ProbeConcentrationUnits']+ \
                str(row['ProbeTime'])+row['ProbeTimeUnits']+\
                row['CellType'] #Same compound treatment and probe treatment means samples share a 'group'
            if row['IsReference'] == "Y":
                refs.append(row['SampleName'])
        except Exception as err:
            logging.error(f"Error adding compound treatment {row}\n{err}")
            db.session.rollback()
        else:
            db.session.commit()

    #replace the compound ids with uuids that define group membership
    # Create a set of unique values in the dictionary
    unique_cpds = set(cpdtrt_groups.values())

    # Generate a uuid for each unique value and create a mapping between old values and new uuids
    mapping = {val: str(uuid.uuid4()) for val in unique_cpds}

    # Replace each value in the original dictionary with its corresponding uuid from the mapping
    for key, val in cpdtrt_groups.items():
        cpdtrt_groups[key] = mapping[val]

    #Read the data file and import each intensity reading linking it to each compound treatment    
    data = pd.read_csv(datafile, dtype={'Position': str}, keep_default_na=False)

    #We will save the mean reference intensity for each peptide/scan from each residue in this dictionary of dictionaries
    #To use at the end to flag which peptide/scan is to be displayed
    mean_ref_intensities = {}

    for _, row in tqdm(data.iterrows(), total=data.shape[0]):

        #Errors here should back out the plex and compound treatments. Otherwise you cannot reload

        logging.info(f"Searching for residue {row['UniProt']} - {row['Position']}")

        #Find the proteins and positions
        uniprots = row['UniProt'].split(";")
        positions = str(row['Position']).split(" ; ")

        if len(uniprots) > len(positions):
            logging.error(f"{row}: Number of UniProt IDs ({len(uniprots)}) not equal to positions ({len(positions)})")
            continue
        elif len(uniprots) < len(positions): #This scenario could be where there are multiple modified positions in the same peptide - not dealing with these currently TODO
            logging.error(f"{row}: Number of UniProt IDs ({len(uniprots)}) not equal to positions ({len(positions)})")
            continue

        for uniprot, position in zip(uniprots, positions):

            protein = db.session.query(Protein).filter(Protein.uniprot == uniprot).first()
            if protein == None:
                logging.warning(f"Protein not found in db: {uniprot}")
                continue

            if position == "NA" or position == "NaN": #We occasionally see unknown residues - this is an error in upstream data and needs fixing
                logging.info(f"Position of residue is NA: {uniprot}:{position}")
                continue
            else:
                try:
                    position=int(position)
                except ValueError:
                    print(f"Error parsing position in {row}")
                    raise

            residue_type = row['ModifiedResidue'].upper()[0:3] #Should always be 'CYS' or other upcase three letter code

            #Find the residue
            residue = db.session.query(Residue).filter(
                Residue.uniprot == uniprot,
                Residue.position == position,
                Residue.type == residue_type
            ).first()
            if residue == None:
                logging.warning(f"Residue - {uniprot}:{position}:{residue_type} not found")
                continue # move to next position if the residue is not found

            db.session.commit()

            #Add an entry to the dictionary of residues if we've not before
            residue_id = f"{uniprot}:{position}"
            if not residue_id in mean_ref_intensities:
                mean_ref_intensities[residue_id] = {}

            #Now go through each sample name from the sample table and add that intensity reading
            #Find reference samples
            if len(refs) == 0:
                logging.error("No references found - check if compound treatments added?")

            ref_values = []
            for ref_sample in refs:
                ref_values.append(row[ref_sample])

            ref_value_mean = stat.mean(ref_values)

            if ref_value_mean == 0:
                logging.info(f"All control channels are 0 - skipping")
                continue

            ref_value_rsd = (stat.stdev(ref_values) / ref_value_mean) * 100
            
            mean_ref_intensities[residue_id][row['Scan']] = ref_value_mean

            logging.info(f"Reference values {ref_values}")
            logging.info(f"Adding itensity readings and CR - mean reference {ref_value_mean} / rsd {ref_value_rsd}")
            
            group_mean_values = {}
            p_values = {}

            for sample in samplenames:

                value = row[sample]
                if value == "NA":
                    value = float('nan')
                else:
                    value = float(value)
                if value == 0:
                    value = 1 #To avoid infinity

                group_id = cpdtrt_groups[sample] # Adds the uuid(value) for the given sample name
                residue_id = str(residue.id)
                scan = row['Scan']
                # Stores group_id and residue_id with the intensity values for each instances in the data table rows
                if group_id not in group_mean_values:
                    group_mean_values[group_id] = {}
                    group_mean_values[group_id][residue_id] = {}
                    group_mean_values[group_id][residue_id][scan] = [value]
                else:
                    if residue_id not in group_mean_values[group_id]:
                        group_mean_values[group_id][residue_id] = {}
                        group_mean_values[group_id][residue_id][scan] = [value]
                    else:
                        if scan not in group_mean_values[group_id][residue_id]:
                            group_mean_values[group_id][residue_id][scan] = [value]
                        else:
                            group_mean_values[group_id][residue_id][scan].append(value)



                intensityreading = IntensityReading(
                    plex = plex,
                    compoundtreatment_id = cpdtrts[sample],
                    residue = residue,
                    scan = str(row['Scan']),
                    value = value,
                    peptideseq = row['Peptide'],
                    modification = row['Modifications'],
                    multimapper = len(uniprots)>1
                )
                competitionratio = CompetitionRatio(
                    plex = plex,
                    compoundtreatment_id = cpdtrts[sample],
                    residue = residue,
                    scan = str(row['Scan']),
                    cr = ref_value_mean / value, 
                    control_rsd = ref_value_rsd,
                    display_flag = True, #Label to say if this CR is the one to display to the end user we default to true and then set to false others later,
                    multimapper = len(uniprots)>1,
                    group_id = group_id,

                )
                db.session.add(intensityreading)
                db.session.add(competitionratio)
            try:
                db.session.commit()
            except:
                logging.info(f"Error adding {row} rolling back")
                db.session.rollback()
                raise Exception

            # Calculates the average for all the intensity values in the list in the dictionary
            for group_id, group_data in group_mean_values.items():
                for residue_id, scans_data in group_data.items():
                    for scan, values in scans_data.items():
                        if values:
                            # Calculate the mean value
                            mean_value = sum(values) / len(values)

                            # Initialize p_value as None
                            p_value = None

                            # Check if there are more than 1 value for the t-test
                            if len(values) > 1:
                                # Perform 2 sample T-test
                                t_stat, p_value = ttest_ind(ref_values, values)

                            # Update the dictionary with the mean value and p-value
                            group_mean_values[group_id][residue_id][scan] = {
                                'mean_value': mean_value,
                                'p_value': p_value
                            }

                            # Query the competition ratio table and filter using residue_id, group_id, and scan
                            comp_group_cr = db.session.query(CompetitionRatio).filter(
                                CompetitionRatio.residue_id == residue_id,
                                CompetitionRatio.group_id == group_id,
                                CompetitionRatio.scan == str(scan)
                            ).all()

                            # Calculate the CR and update group_cr and p_value columns
                            for comp in comp_group_cr:
                                comp.group_cr = ref_value_mean / mean_value
                                comp.p_value = p_value
                                comp.replicate_no = len(values)

                            db.session.commit()


    #Go through each residue in mean_ref_intensities
    #If there's only one Scan id (key) then all is fine. Otherwise find the scans that are not max and update them to display = False
    for residue in tqdm(mean_ref_intensities):
        max_intens = 0
        max_scan = None
        if len(mean_ref_intensities[residue]) > 1:
            logging.info(f"Residue {residue} has multiple scans/peptides")
            for scan in mean_ref_intensities[residue]:
                if mean_ref_intensities[residue][scan] > max_intens:
                    max_intens = mean_ref_intensities[residue][scan]
                    max_scan = scan
            #Find residue with this ID
            residue = db.session.query(Residue).filter(
                Residue.uniprot == residue.split(':')[0],
                Residue.position == residue.split(':')[1],
            ).first()
            logging.info(f"Scan {max_scan} has highest mean intensity in reference channels - setting all other scans to not disply")
            #find all the Crs for this residue in this plex and mark as display false unless they are the scan with the max intensity
            for cr in db.session.query(CompetitionRatio).filter(
                CompetitionRatio.plex == plex,
                CompetitionRatio.residue == residue,
            ):
                if not cr.scan == max_scan:
                    logging.info(f"Setting scan {cr.scan}-{cr.compoundtreatment_id} to false")
                    cr.display_flag = False
                    db.session.commit()

@databp.cli.command("sbio_import_ss")
@click.option('-d', '--dataroot', required=True)
@click.option('-l', '--logfile', default='discoverome_sbio_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def sbio_import(dataroot, logfile, drop_all):
    """Imports discoverome structural biology information"""

    from models.chemopro import db, Protein, Residue, Structure, StructureChain, StructureResidue, Ligand, LigandResidueDistance

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)

    #read through the three types of files in the root directory - disulfides, accessbility and ligand distance
    #disulfide files are only there for CYS but form a root of the analysis
    disulfide_files = glob.glob(f"{dataroot}/*/*_SS_CYS.csv")
    residuetype = "CYS"

    for f in tqdm(disulfide_files, smoothing=0):

        structure_code = os.path.basename(f).split('_')[0]

        #Make structure
        logging.info(f"Making structure - {structure_code}")
        structure = db.session.query(Structure).filter(Structure.id == structure_code)
        if not db.session.query(structure.exists()).scalar():
            structure = Structure(id = structure_code)
            db.session.add(structure)
            db.session.commit()
        else:
            structure = structure.first()

        logging.info(f"Processing {f} as {structure_code} SS info")

        data = pd.read_csv(f, dtype={'PDBChain1':str, 'PDBChain2':str})
        for index, row in data.iterrows():

            logging.info(f"Processing row: {','.join((str(x) for x in row.to_list()))}")

            #Skip any rows that do not map to UniProt
            if row['UPAccession1'] == "-" or row['UPAccession2'] == "-":
                logging.info(f"Skipping row {index} due to missing UniProt")
                continue

            #Make protein
            protein1 = db.session.query(Protein).filter(Protein.uniprot == row["UPAccession1"])
            if not db.session.query(protein1.exists()).scalar():
                logging.warning(f"Protein not found in db: {row['UPAccession1']}")
                db.session.rollback()
                continue
            else:
                protein1 = protein1.first()

            if row["UPAccession1"] != row["UPAccession2"]:
                protein2 = db.session.query(Protein).filter(Protein.uniprot == row["UPAccession2"])
                if not db.session.query(protein2.exists()).scalar():
                    logging.warning(f"Protein not found in db: {row['UPAccession2']}")
                    db.session.rollback()
                    continue
                else:
                    protein2 = protein2.first()
            else:
                protein2 = protein1

            #Make chains
            chain1 = db.session.query(StructureChain).filter(StructureChain.structure == structure, StructureChain.chain == row['PDBChain1'])
            if not db.session.query(chain1.exists()).scalar():
                chain1 = StructureChain(structure = structure, chain = row["PDBChain1"], uniprot_id = row["UPAccession1"])
                db.session.add(chain1)
            else:
                chain1 = chain1.first()

            if row["PDBChain1"] != row["PDBChain2"]:
                chain2 = db.session.query(StructureChain).filter(StructureChain.structure == structure, StructureChain.chain == row['PDBChain2'])
                if not db.session.query(chain2.exists()).scalar():
                    chain2 = StructureChain(structure = structure, chain = row["PDBChain2"], uniprot_id = row["UPAccession2"])
                    db.session.add(chain2)
                else:
                    chain2 = chain2.first()
            else:
                chain2 = chain1

            #Make residue1
            res1 = db.session.query(Residue).filter(Residue.uniprot == row["UPAccession1"], Residue.position == row['UPResidue1'])
            if not db.session.query(res1.exists()).scalar():
                logging.info(f"Residue not found in database: {row['UPAccession1']}:{row['UPResidue1']}")
                db.session.rollback()
                continue
            else:
                res1 = res1.first()

            #Make residue 2
            res2 = db.session.query(Residue).filter(Residue.uniprot == row["UPAccession2"], Residue.position == row['UPResidue2'])
            if not db.session.query(res2.exists()).scalar():
                logging.info(f"Residue not found in database: {row['UPAccession2']}:{row['UPResidue2']}")
                db.session.rollback()
                continue
            else:
                res2 = res2.first()

            #Make structure residue1
            sres1 = db.session.query(StructureResidue).filter(
                StructureResidue.structure_id == structure_code, 
                StructureResidue.chain == chain1, StructureResidue.position == row['UPResidue1']
            )
            if not db.session.query(sres1.exists()).scalar():
                sres1 = StructureResidue(
                    structure_id = structure_code, chain = chain1, residue = res1, 
                    pdb_position = row['PDBResidue1'], type = residuetype, 
                    confidence = 999, in_disulfide = True
                )
                db.session.add(sres1)
            else:
                sres1 = sres1.first()

            sres2 = db.session.query(StructureResidue).filter(
                StructureResidue.structure_id == structure_code, 
                StructureResidue.chain == chain2, StructureResidue.position == row['UPResidue2'])
            if not db.session.query(sres2.exists()).scalar():
                sres2 = StructureResidue(
                    structure_id = structure_code, chain = chain2, residue = res2, 
                    pdb_position = row['PDBResidue2'], type = residuetype, 
                    confidence = 999, in_disulfide = True
                )
                db.session.add(sres2)
            else:
                sres2 = sres2.first()

            db.session.commit()

@databp.cli.command("sbio_import_sasa")
@click.option('-d', '--dataroot', required=True)
@click.option('-r', '--residuetype', required=True)
@click.option('-l', '--logfile', default='discoverome_sbio_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def sbio_import(dataroot, logfile, residuetype, drop_all):
    """Imports discoverome structural biology information"""

    from models.chemopro import db, Protein, Residue, Structure, StructureChain, StructureResidue, Ligand, LigandResidueDistance

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)

    #read through the three types of files in the root directory - disulfides, accessbility and ligand distance
    #disulfide files are only there for CYS but form a root of the analysis
    files = glob.glob(f"{dataroot}/*/*_SASA_{residuetype}.csv")

    for f in tqdm(files, smoothing=0):

        structure_code = os.path.basename(f).split('_')[0]
        logging.info(f"Processing {f} as {structure_code} SASA info")

        if not os.path.exists(f):
            logging.warn(f"Missing SASA file {f} for {structure_code}")
            continue

        #Make structure
        structure = db.session.query(Structure).filter(Structure.id == structure_code)
        if not db.session.query(structure.exists()).scalar():
            logging.info(f"Making structure - {structure_code}")
            structure = Structure(id = structure_code)
            db.session.add(structure)
            db.session.commit()
        else:
            structure = structure.first()
            
        data = pd.read_csv(f, dtype={'PDBChain':str})
        for index, row in data.iterrows():

            #Hack TO AVOID NAs in chain IDs (that are called 'NA' like in 6WAT)
            if not isinstance(row['PDBChain'], str):
                row['PDBChain'] = "NA"

            logging.info(f"Processing row: {','.join((str(x) for x in row.to_list()))}")

            #Skip any rows that do not map to UniProt
            if row['UPAccession'] == "-":
                logging.info(f"Skipping row {index} due to missing UniProt")
                continue

            #Make protein
            protein1 = db.session.query(Protein).filter(Protein.uniprot == row["UPAccession"])
            if not db.session.query(protein1.exists()).scalar():
                logging.warning(f"Protein not found in db: {row['UPAccession']}")
                db.session.rollback()
                continue
            else:
                protein1 = protein1.first()

            #Make chain
            chain1 = db.session.query(StructureChain).filter(StructureChain.structure == structure, StructureChain.chain == row['PDBChain'])
            if not db.session.query(chain1.exists()).scalar():
                chain1 = StructureChain(structure = structure, chain = row["PDBChain"], uniprot_id = row["UPAccession"])
                db.session.add(chain1)
            else:
                chain1 = chain1.first()

            #Make residue
            res1 = db.session.query(Residue).filter(Residue.uniprot == row["UPAccession"], Residue.position == row['UPResidue'])
            if not db.session.query(res1.exists()).scalar():
                logging.info(f"Residue not found in database: {row['UPAccession']}:{row['UPResidue']}")
                db.session.rollback()
                continue
            else:
                res1 = res1.first()

            #Make structure residue1
            sres1 = db.session.query(StructureResidue).filter(
                StructureResidue.structure_id == structure_code, 
                StructureResidue.chain == chain1, StructureResidue.pdb_position == str(row['PDBResidue'])
            )
            if not db.session.query(sres1.exists()).scalar():
                logging.info(f"StructureResidue not found so creating new")
                sres1 = StructureResidue(
                    structure_id = structure_code, chain = chain1, residue = res1, 
                    pdb_position = row['PDBResidue'], #type = residuetype, 
                    confidence = 999, in_disulfide = False, accessibility = row["SASA"], depth = row["Depth"]
                )
                db.session.add(sres1)
            else:
                sres1 = sres1.first()
                logging.info(f"StructureResidue found so updating")
                sres1.accessibility = row['SASA']
                sres1.depth = row['Depth']

            db.session.commit()

@databp.cli.command("sbio_import_ligand_distance")
@click.option('-d', '--dataroot', required=True)
@click.option('-r', '--residuetype', required=True)
@click.option('-l', '--logfile', default='discoverome_sbio_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def sbio_import(dataroot, logfile, residuetype, drop_all):
    """Imports discoverome structural biology information"""

    from models.chemopro import db, Protein, Residue, Structure, StructureChain, StructureResidue, Ligand, LigandResidueDistance

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)

    #read through the three types of files in the root directory - disulfides, accessbility and ligand distance
    #disulfide files are only there for CYS but form a root of the analysis
    files = glob.glob(f"{dataroot}/*/*_Ligand_{residuetype}.csv")

    for f in tqdm(files, smoothing=0):

        structure_code = os.path.basename(f).split('_')[0]
        logging.info(f"Processing {f} as {structure_code} ligand info")

        if not os.path.exists(f):
            logging.warn(f"Missing ligand info file {f} for {structure_code}")
            continue

        #Make structure
        structure = db.session.query(Structure).filter(Structure.id == structure_code)
        if not db.session.query(structure.exists()).scalar():
            logging.info(f"Making structure - {structure_code}")
            structure = Structure(id = structure_code)
            db.session.add(structure)
            db.session.commit()
        else:
            structure = structure.first()

        data = pd.read_csv(f, dtype={'PDBChain':str, 'LigandID':str, 'LigandChain':str})
        for index, row in data.iterrows():

            #Hack TO AVOID NAs in chain IDs (that are called 'NA' like in 6WAT)
            if not isinstance(row['PDBChain'], str):
                row['PDBChain'] = "NA"
            if not isinstance(row['LigandChain'], str):
                row['LigandChain'] = "NA"
            if not isinstance(row['LigandID'], str):
                row['LigandID'] = "NA"

            logging.info(f"Processing row: {','.join((str(x) for x in row.to_list()))}")

            #Skip any rows with distance > 10A
            distance_limit = 10
            if row['Distance'] > distance_limit:
                logging.info(f"Skipping row {index} due to distance limit ({distance_limit})")
                continue

            #Skip any rows that do not map to UniProt
            if row['UPAccession'] == "-":
                logging.info(f"Skipping row {index} due to missing UniProt")
                continue

            #Make protein
            protein1 = db.session.query(Protein).filter(Protein.uniprot == row["UPAccession"])
            if not db.session.query(protein1.exists()).scalar():
                logging.warning(f"Protein not found in db: {row['UPAccession']}")
                db.session.rollback()
                continue
            else:
                protein1 = protein1.first()

            #Make chain
            chain1 = db.session.query(StructureChain).filter(StructureChain.structure == structure, StructureChain.chain == row['PDBChain'])
            if not db.session.query(chain1.exists()).scalar():
                chain1 = StructureChain(structure = structure, chain = row["PDBChain"], uniprot_id = row["UPAccession"])
                db.session.add(chain1)
            else:
                chain1 = chain1.first()

            #Make residue
            res1 = db.session.query(Residue).filter(Residue.uniprot == row["UPAccession"], Residue.position == row['UPResidue'])
            if not db.session.query(res1.exists()).scalar():
                logging.info(f"Residue not found in database: {row['UPAccession']}:{row['UPResidue']}")
                db.session.rollback()
                continue
            else:
                res1 = res1.first()

            #Make structure residue1
            sres1 = db.session.query(StructureResidue).filter(
                StructureResidue.structure_id == structure_code, 
                StructureResidue.chain == chain1, StructureResidue.pdb_position == str(row['PDBResidue'])
            )
            if not db.session.query(sres1.exists()).scalar():
                logging.info(f"Could not find StructureResidue - should not happen since SASA should cover all.")
                sres1 = StructureResidue(
                    structure_id = structure_code, chain = chain1, residue = res1, 
                    pdb_position = row['PDBResidue'], #type = residuetype, 
                    confidence = 999, in_disulfide = False 
                )
                db.session.add(sres1)
            else:
                sres1 = sres1.first()

            #Make Ligand
            ligand = db.session.query(Ligand).filter(
                Ligand.structure_id == structure_code, 
                Ligand.code == row['LigandID'], Ligand.chain == row['LigandChain']
            )
            if not db.session.query(ligand.exists()).scalar():
                ligand = Ligand(structure_id = structure_code, chain = row['LigandChain'], code = row["LigandID"])
                db.session.add(ligand)
            else:
                ligand = ligand.first()

            #Make ligand-residue distance - sometimes the same ligand seems to appear twice e.g. ACT in 5ezd (two in chain B). In this case report the shortest distance
            ligandresiduedist = db.session.query(LigandResidueDistance).filter(
                LigandResidueDistance.ligand_id == ligand.id,
                LigandResidueDistance.structureresidue_id == sres1.id
            )
            if not db.session.query(ligandresiduedist.exists()).scalar():
                ligandresiduedist = LigandResidueDistance(
                    ligand = ligand, 
                    residue = sres1, 
                    distance = row['Distance'])
                db.session.add(ligandresiduedist)
            else:
                ligandresiduedist = ligandresiduedist.first()
                if row['Distance'] < ligandresiduedist.distance:
                    ligandresiduedist.distance = row['Distance']

            db.session.commit()


#This will only work on structures already imported because as currently speced the fpocket output does not include
#explicit mapping to uniprot IDs
#Note that fpocket uses naive PDB residue numbering so we need to use SIFTS to map them to the correct Structure Residues
@databp.cli.command("fpocket_import")
@click.option('-d', '--dataroot', required=True)
@click.option('-s', '--siftsroot', required=False)
@click.option('-l', '--logfile', default='discoverome_sbio_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def sbio_import(dataroot, siftsroot, logfile, drop_all):
    """Imports discoverome pocket information"""

    from models.chemopro import db, Structure, StructureChain, StructureResidue, Pocket, PocketResidue
 
    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)

    structure_codes = glob.glob(f"{dataroot}/*/*")

    for f in tqdm(structure_codes, smoothing=0):

        structure_code = os.path.split(f)[-1]
        structure_abrev = structure_code[1:3]

        #Find structure
        logging.info(f"Making pockets for structure - {structure_code} - {structure_abrev}")
        structure = db.session.query(Structure).filter(Structure.id == structure_code)
        if not db.session.query(structure.exists()).scalar():
            logging.warn(f"No structure {structure_code} found - need to create first to get uniprot mapping")
            continue
        else:
            structure = structure.first()

        logging.info(f"Processing {f} as {structure_code}")

        if siftsroot:
            sifts_file = f"{siftsroot}/{structure_abrev}/{structure_code}.xml.gz"
            residue_map = {}  # To store map from naive residue numbering to true PDb numbering

            logging.info(f"Parsing SIFTS: {sifts_file}")
            if os.path.exists(sifts_file):
                try:
                    sift = ET.parse(gzip.open(sifts_file, 'r'))
                    for chain in sift.findall(".//entity[@type='protein']", namespaces={"": "http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"}):
                        for residue in chain.findall(".//residue", namespaces={"": "http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"}):
                            pdb_res = residue.findall(".//crossRefDb[@dbSource='PDB']", namespaces={"": "http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd"})[0]
                            residue_map[f"{chain.attrib['entityId']}::{residue.attrib['dbResNum']}"] = [pdb_res.attrib['dbResNum'], pdb_res.attrib['dbChainId']]
                except Exception as e:
                    logging.error(f"Error parsing SIFTS file {sifts_file}: {e}")
            else:
                logging.warning(f"SIFTS file {sifts_file} not found. Skipping SIFTS parsing.")
                # Optionally: continue to the next iteration if SIFTS parsing is critical
                # continue

        pocket_files = glob.glob(f"{f}/pockets/*_atm.cif")
        for pocket_file in pocket_files:
            logging.info(f"Processing: {pocket_file}")

            # Then look for confidence file
            af2_confidence_file = pocket_file.replace("_atm.cif","_conf.csv")

            pocket_confidence_data = None
            if os.path.exists(af2_confidence_file):
                logging.info(f"Reading confidence data from {af2_confidence_file}")
                pocket_confidence_data = pd.read_csv(af2_confidence_file).iloc[0]

            pocket_id = os.path.basename(pocket_file).split("_")[0]

            # CIF parsing code
            pocket_data = {}
            residue_data = set()
            if os.path.exists(pocket_file):
                with open(pocket_file) as pock_f:
                    lines = pock_f.readlines()
                    for line in lines:
                        pocket_match = re.search("^\d+\s+-(.+?)\s+:\s+(\S+)$", line)
                        if pocket_match:
                            pocket_data[pocket_match.group(1).lstrip()] = float(pocket_match.group(2))
                        residue_match = re.search("^ATOM.+?(?:[CL]YS|SER)\s{3}(.).+?(\d+)", line)  # Match the chain ID (group1) and naive PDB position (group2) of every CYS/LYS atom
                        if residue_match:
                            residue_data.add(f"{residue_match.group(1)}::{residue_match.group(2)}")
            else:
                logging.warning(f"Pocket file {pocket_file} not found. Skipping this file.")
                continue

            if len(residue_data) == 0:
                logging.info(f"No cysteines or lysines found in pocket {pocket_id}")
                continue

            # Need to check if pocket already exists? Probably always drop and restart
            pocket = Pocket(
                structure_id=structure_code,
                pocket_id=pocket_id,
                pocket_score=pocket_data.get('Pocket Score'),
                drug_score=pocket_data.get('Drug Score'),
                pocket_volume_MC=pocket_data.get('Pocket volume (Monte Carlo)'),
                pocket_volume_hull=pocket_data.get('Pocket volume (convex hull)')
            )
            if isinstance(pocket_confidence_data, pd.Series):
                pocket.mean_confidence = pocket_confidence_data['MeanConfidence']
                pocket.median_confidence = pocket_confidence_data['MedianConfidence']
                pocket.min_confidence = pocket_confidence_data['MinConfidence']

            db.session.add(pocket)
            db.session.flush()
            db.session.refresh(pocket)

            for res_data in residue_data:
                logging.info(f"Creating link to pocket {pocket.pocket_id} to residue {res_data}")

                # Need to map from naive PDb numbering to correct PDB numbering (not UniProt) here if SIFTS done
                if siftsroot and res_data in residue_map:
                    pdb_position = residue_map[res_data][0]
                else:
                    pdb_position = res_data.split("::")[1]

                res = db.session.query(StructureResidue).filter(
                    StructureResidue.pdb_position == pdb_position,
                    StructureResidue.structure_id == structure_code,
                    StructureResidue.chain_id == StructureChain.id,
                    StructureChain.chain == res_data.split("::")[0]
                )
                if not db.session.query(res.exists()).scalar():
                    logging.warning(f"No residue found for {structure_code} - {res_data.split('::')[0]} - {res_data.split('::')[1]}")
                    continue
                elif res.count() > 1:
                    logging.warning(f"Multiple residues found for {structure_code} - {res_data.split('::')[0]} - {res_data.split('::')[1]}")
                    continue
                else:
                    res = res.first()

                pocketresidue = PocketResidue(structureresidue_id=res.id, pocket_id=pocket.id)
                db.session.add(pocketresidue)

        db.session.commit()


@databp.cli.command("ligand_import")
@click.option('-d', '--datafile', required=True)
@click.option('-l', '--logfile', default='discoverome_pdbligand_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def ligand_import(datafile, logfile, drop_all):
    """Imports PDB ligand information"""

    from models.chemopro import db, Ligand
 
    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)

    df = pd.read_csv(datafile, sep="\t", header=None, dtype={0:str, 1:str, 2:str, 3:float, 4:str, 5:str})

    for _, row in tqdm(df.iterrows(), total=df.shape[0]):
        smiles = row[0]
        code = str(row[1]).upper() #Not sure why this is an issue
        name = row[2]
        mw = row[3]
        inchi = row[4]
        chembl = row[5]

        ligands = db.session.query(Ligand).filter(Ligand.code == code)
        logging.info(f"Updating {code} with data {row}")
        for lig in ligands:
            lig.smiles = smiles
            lig.name = name
            lig.mw = mw
            lig.inchi = inchi
            lig.chembl = chembl

        db.session.commit()

@databp.cli.command("protein_data_import")
@click.option('-pf', '--proteinfile', required = True)
@click.option('-l', '--logfile', default = 'discoverome_protein_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def protein_data_import(proteinfile, logfile, drop_all):
    "Import protein synonyms and the position of all residues"

    from models.chemopro import db, ProteinSynonym, Protein, Residue

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)

    df = pd.read_csv(proteinfile, sep="\t", dtype={'Entry': str, 'Gene Names (primary)': str, 'Gene Names (synonym)': str}, keep_default_na=False)
    #df['Gene Names'] = df['Gene Names'].str.split(' ') # Split the values in the synonym column using the space between each values
    #df = df.explode('Gene Names') # Creates new row for each synonym retaining the same data as the original  
    df['Function [CC]'] = df['Function [CC]'].astype(str).map(lambda x:x.lstrip("FUNCTION:")) # Remove the 'FUNCTION:' string from function column
    df.rename(columns={'Entry':'uniprot', 'Gene Names (primary)': 'symbol', 'Gene Names (synonym)': 'synonyms', 'Function [CC]' : 'description'}, inplace= True)

    for _,row in tqdm(df.iterrows(), total=df.shape[0]):

        #check if protein already exists
        protein = db.session.query(Protein).filter(Protein.uniprot==row['uniprot']).first()

        #add Protein with symbol and description
        if protein:
            protein.symbol = row['symbol']
            protein.description = row['description']
            db.session.commit()
        else:
            protein = Protein(uniprot=row['uniprot'], symbol=row['symbol'], description=row['description'])
            db.session.add(protein)
            db.session.commit()

        #Go through each synonym and add to synonym table
        synonyms = row['synonyms'].split(' ')
        for synonym in synonyms:
            p_synonym = ProteinSynonym(uniprot=row['uniprot'], synonym=synonym, type='uniprot_name')
            logging.info(f"Creating new protein synonym for {row['uniprot']} - {synonym}")
            db.session.add(p_synonym)
        db.session.commit()

        #Go through sequence to find each CYS and add as a Residue
        for i in range(len(row['Sequence'])):
            if row['Sequence'][i] == 'C':
                residue = Residue(uniprot=row['uniprot'], position=i+1, type="CYS")
                db.session.add(residue)
            if row['Sequence'][i] == 'K':
                residue = Residue(uniprot=row['uniprot'], position=i+1, type="LYS")
                db.session.add(residue)
            if row['Sequence'][i] == 'S':
                residue = Residue(uniprot=row['uniprot'], position=i+1, type="SER")
                db.session.add(residue)
        db.session.commit()


@databp.cli.command("smiles_import")
@click.option('-sf', '--smilefile', required = True)
@click.option('-l', '--logfile', default = 'smiles_import.log')
def smiles_import(smilefile, logfile):
    "Import smiles string to compound table"

    from models.chemopro import db, Compound

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    df_smiles = pd.read_csv(smilefile, sep=",") # Load file in as a pandas dataframe
    df_smiles= df_smiles.filter(['ID','SMILES', 'Batch ID', "Formatted ID"], axis = 1) # Filter the table to have only 'ID' and 'Smiles' columns 

    for _,row in tqdm(df_smiles.iterrows(), total=df_smiles.shape[0]): # For the loading function 
        compound_id = row["ID"] # get the compound id's 
        smiles = row["SMILES"] # get the smiles string 
        batch_id = row["Batch ID"]
        formatted_id = row["Formatted ID"]

        compound_smiles = db.session.query(Compound).filter(or_(Compound.id == compound_id, Compound.id.like(f'%{batch_id}%'), Compound.id.like(f'%{formatted_id}%'))) # query the database to see if 'id' in table matched 'id' in dataframe
        if compound_smiles.count() == 0: # if not matches, entry into the log file
            logging.info(f"Compound{compound_id} does not exist in the compound table ")

        else: # If compound 'id' matches, add the smiles string to the database compound table
            logging.info(f"Updating compound {compound_id} with SMILES string {smiles}")
            for c in compound_smiles:
                c.smiles = smiles
                db.session.commit()

@databp.cli.command("cpd_images")
@click.option('-l', '--logfile', default = 'cpd_images.log')
def smiles_import(logfile):
    "Calculate images from SMILES, save into folder and import path to compound table"

    from models.chemopro import db, Compound

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    compounds = db.session.query(Compound).all() # Query compound table 

    for s in compounds: # loop through each compound in the table 

        if s.image is None and s.smiles is not None:
            logging.info(f"For Compound {s.id} creating new SMILES strucutre {s.smiles}")            
            mol = Chem.MolFromSmiles(s.smiles) # create a chemical strucutre from the smiles string
            img = Draw.MolToImage(mol) #draw the strucutre 
            file = f'/resources/smiles-structure/{s.id}.png'
            if os.path.exists(file):
                logging.info(f"File for {s.id} already exists")
            else:
                img.save(file)

            s.image = file
            logging.info(f"Adding image path of {s.id} to image field")
            db.session.commit()

        else:
            if s.image is not None:
                logging.info(f"Image path already exists for {s.id}")
            else:
                logging.info(f"No smiles entry for {s.id}")

@databp.cli.command("target_list_import")
@click.option('-f', '--file', required = True) 
@click.option('-n', '--listname', type=str, required = True) # List name is required unlike description
@click.option('-d', '--description', type =str )
@click.option('-l', '--logfile', default = 'target_list_import.log')
def target_list_import(file, listname, description, logfile): # Function takes in file, name of list and description of list
    "Import target list data"

    from models.chemopro import db, TargetList, ProteinToList, Protein
    
    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    new_list = TargetList(name=listname, description = description) # Creates new list in the target list table taking the list name and description if provided 
    db.session.add(new_list)
    db.session.flush()

    logging.info(f"Created new protein list {listname} - {description} ({new_list.id})")

    with open(file, 'r') as f: # opens the file with uniprot id list
        for line in f:
            uniprot_id = line.strip()

            logging.info(f"Adding {uniprot_id} to list {new_list.id}")
            protein = db.session.query(Protein).filter(Protein.uniprot==uniprot_id).first()
            if protein:
                logging.info(f"{uniprot_id} found in database")
                new_list.proteins.append(protein)
            else:
                logging.info(f"Protein not found in database: {uniprot_id}")
            
    db.session.commit()


@databp.cli.command("residue_list_import")
@click.option('-f', '--file', required = True) 
@click.option('-n', '--listname', type=str, required = True) # List name is required unlike description
@click.option('-d', '--description', type =str )
@click.option('-l', '--logfile', default = 'residue_list_import.log')
def residue_list_import(file, listname, description, logfile): # Function takes in file, name of list and description of list
    "Import residue list data"

    from models.chemopro import db, ResidueList, ResidueToList, Residue
    
    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    list = ResidueList(name=listname, description = description) # Creates new list in the target list table taking the list name and description if provided 
    db.session.add(list)
    db.session.flush()

    with open(file, 'r') as f: # opens the file with uniprot id list
        for line in f:
            line = line.strip()
            uniprot_id,position = line.split(",")
            position = int(position)

            logging.info(f"Adding {uniprot_id}:{position} to list {list.id}")
            residue = db.session.query(Residue).filter(Residue.uniprot==uniprot_id, Residue.position==position).first()
            if residue:
                logging.info(f"{uniprot_id}:{position} found in database")
                list.residues.append(residue)
            else:
                logging.info(f"Residue not found in database: {uniprot_id}:{position}")
            
    db.session.commit()

@databp.cli.command("ensembl_import")
@click.option('-f', '--file', required=True) 
@click.option('-l', '--logfile', default='ensembl_import.log')
def ensembl_import(file, logfile):
    "Import ensembl id for uniprot ids"
    
    from models.chemopro import db, Protein, ProteinSynonym
    
    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )
    
    df = pd.read_csv(file, sep='\t')
    
    for _, row in tqdm(df.iterrows(), total=df.shape[0]):
        uniprot = row['uniprot']
        ensembl = row['ensembl']
        data_type = 'ensembl_ID'
        
        # Check if the ProteinSynonym row already exists
        existing_synonym = db.session.query(ProteinSynonym).filter_by(synonym=ensembl, type=data_type).first()
        if existing_synonym is not None:
            logging.info(f"Synonym '{ensembl}' already exists for uniprot '{existing_synonym.uniprot}'")
            continue
        
        # Add the ProteinSynonym row if it does not exist
        protein = db.session.query(Protein).filter_by(uniprot=uniprot).first()
        if protein is None:
            logging.warning(f"No protein found with uniprot '{uniprot}'")
            continue
        
        synonym = ProteinSynonym(uniprot=uniprot, synonym=ensembl, type=data_type)
        db.session.add(synonym)
        db.session.commit()
        
        logging.info(f"Added synonym '{ensembl}' for uniprot '{uniprot}'")

@databp.cli.command("residue_feature_import")
@click.option('-pf', '--proteinfile', required = True)
@click.option('-l', '--logfile', default = 'discoverome_residue_feature_import.log')
def residue_feature_import(proteinfile, logfile):
    "Import UniProt residue features"

    from models.chemopro import db, Protein, Residue, ResidueFeature

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    df = pd.read_csv(proteinfile, sep="\t", dtype={'Entry': str, 'Gene Names (primary)': str, 'Gene Names (synonym)': str}, keep_default_na=False)
    for _,row in tqdm(df.iterrows(), total=df.shape[0]):

        #Check if protein exists
        protein = db.session.query(Protein).filter(Protein.uniprot==row['Entry']).first()
        if not protein:
            logging.info(f"Protein {row['Entry']} not found")
            continue

        #Add residue features
        #Active site	Binding site	DNA binding	Mutagenesis
        
        active_site = row["Active site"]
        if active_site != "":
            active_sites = active_site.split("ACT_SITE ")
            for site in active_sites[1:]:
                try:
                    position = int(site.split(";")[0])
                except ValueError:
                    print(f"Error: Unable to convert '{site.split(';')[0]}' to integer.")
                    continue  # Skip this iteration and move to the next one

                # Check if residue exists
                residue = db.session.query(Residue).filter(Residue.uniprot == row['Entry'], Residue.position == position).first()

                if residue:
                    feature = ResidueFeature(residue_id = residue.id, description = site, source = "UniProt ACT_SITE")
                    logging.info(f"{row['Entry']}:{position} adding feature: {site}")
                    db.session.add(feature)
                    db.session.commit()
                else:
                    logging.info(f"{row['Entry']}:{position} not in database")

        binding_site = row["Binding site"]
        if binding_site != "":
            binding_sites = binding_site.split("BINDING ")
            try:
                for site in binding_sites[1:]:
                    if ".." in site.split(";")[0]:
                        start = int(site.split(";")[0].split("..")[0])
                        end   = int(site.split(";")[0].split("..")[1])
                        positions = list(range(start, end+1))
                    else:
                        positions = [int(site.split(";")[0])]
                    #Check if residue exists
                    for position in positions:
                        residue = db.session.query(Residue).filter(Residue.uniprot==row['Entry'],Residue.position==position).first()
                        if residue:
                            feature = ResidueFeature(residue_id = residue.id, description = site, source = "UniProt BINDING")
                            logging.info(f"{row['Entry']}:{position} adding feature: {site}")
                            db.session.add(feature)
                        else:
                            logging.info(f"{row['Entry']}:{position} not in database")
            except ValueError:
                logging.info(f"Error processing {binding_site}")

        binding_site = row["DNA binding"]
        if binding_site != "":
            binding_sites = binding_site.split("DNA_BIND ")
            try:
                for site in binding_sites[1:]:
                    if ".." in site.split(";")[0]:
                        start = int(site.split(";")[0].split("..")[0])
                        end   = int(site.split(";")[0].split("..")[1])
                        positions = list(range(start, end+1))
                    else:
                        positions = [int(site.split(";")[0])]
                    #Check if residue exists
                    for position in positions:
                        residue = db.session.query(Residue).filter(Residue.uniprot==row['Entry'],Residue.position==position).first()
                        if residue:
                            feature = ResidueFeature(residue_id = residue.id, description = site, source = "UniProt DNA_BIND")
                            logging.info(f"{row['Entry']}:{position} adding feature: {site}")
                            db.session.add(feature)
                        else:
                            logging.info(f"{row['Entry']}:{position} not in database")
            except ValueError:
                logging.info(f"Error processing {binding_site}")

        mutated_site = row["Mutagenesis"]
        if mutated_site != "":
            mutated_sites = mutated_site.split("MUTAGEN ")
            try:
                for site in mutated_sites[1:]:
                    if ".." in site.split(";")[0]:
                        start = int(site.split(";")[0].split("..")[0])
                        end   = int(site.split(";")[0].split("..")[1])
                        positions = list(range(start, end+1))
                    else:
                        positions = [int(site.split(";")[0])]
                    #Check if residue exists
                    for position in positions:
                        residue = db.session.query(Residue).filter(Residue.uniprot==row['Entry'],Residue.position==position).first()
                        if residue:
                            feature = ResidueFeature(residue_id = residue.id, description = site, source = "UniProt MUTAGEN")
                            logging.info(f"{row['Entry']}:{position} adding feature: {site}")
                            db.session.add(feature)
                        else:
                            logging.info(f"{row['Entry']}:{position} not in database")
            except ValueError:
                logging.info(f"Error processing {mutated_site}")
        
    db.session.commit()

@databp.cli.command("delete_plex")
@click.option('-p', '--plexid', required=True)
@click.option('-l', '--logfile', default='discoverome_delete_plex.log')
@click.option('--dry-run/--no-dry-run', default=False)
def delete_plex(plexid, logfile, dry_run):
    """Imports discoverome structural biology information"""

    from models.chemopro import db, Plex, CompetitionRatio, IntensityReading, CompoundTreatment

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    logging.info(f"Attempting to delete plex '{plexid}'")

    plex = db.session.query(Plex).filter(Plex.id == plexid).first()
    if not plex:
        logging.info(f"Could not find plex '{plexid}' to delete - aborting")
        exit()
    
    crs = db.session.query(CompetitionRatio).filter(CompetitionRatio.plex_id == plexid).all()
    logging.info(f"Found {len(crs)} CRs to delete from '{plexid}'")

    irs = db.session.query(IntensityReading).filter(IntensityReading.plex_id == plexid).all()
    logging.info(f"Found {len(irs)} IRs to delete from '{plexid}'")

    cts = db.session.query(CompoundTreatment).filter(CompoundTreatment.plex_id == plexid).all()
    logging.info(f"Found {len(cts)} compound treatments to delete from '{plexid}'")

    logging.info(f"Deleting CRs")
    db.session.query(CompetitionRatio).filter(CompetitionRatio.plex_id == plexid).delete()
    logging.info(f"Deleting IRs")
    db.session.query(IntensityReading).filter(IntensityReading.plex_id == plexid).delete()
    logging.info(f"Deleting compound treatments")
    db.session.query(CompoundTreatment).filter(CompoundTreatment.plex_id == plexid).delete()
    logging.info(f"Deleting plex '{plexid}'")
    db.session.query(Plex).filter(Plex.id == plexid).delete()

    if not dry_run:
        logging.info("Committing changes")
        db.session.commit()
    else:
        logging.info("Dry run - no deletions made")




@databp.cli.command("external_compound_data")
@click.option('-crf', '--crfile', required=True) 
@click.option('-comf', '--comfile') 
@click.option('-l', '--logfile', default='external_compound_data.log')
def external_compound_data(crfile, comfile, logfile):
    "Import external data to populate compound treatment and competition ratio table"
    
    from models.chemopro import db, Protein, Residue, Plex, Compound, CompoundTreatment, CompetitionRatio, Experiment, CellType
    
    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )
    
    # Adding experiment ID and description to the experiment table 
    exp = 'Kuljanin_screen'
    desc = 'Reimagining high-throughput profiling of reactive cysteines for cell-based screening of large electrophile libraries' 

    # Check if kul_experiment already exists
    kul_experiment = db.session.query(Experiment).get(exp)

    if kul_experiment is None:
        # kul_experiment doesn't exist, create a new entry
        kul_experiment = Experiment(id=exp, description=desc)
        logging.info(f"Adding experiment {exp}")
        db.session.add(kul_experiment)
    else:
        logging.info(f"Experiment {exp} already exists in the database")

    # Adding Plex 
    plex_id = 'Kuljanin_plex'
    exp_id = exp

    # Check if kul_plex already exists
    kul_plex = db.session.query(Plex).get(plex_id)

    if kul_plex is None:
        # kul_plex doesn't exist, create a new entry
        kul_plex = Plex(id=plex_id, experiment_id=exp_id)
        logging.info(f"Adding plex {plex_id}")
        db.session.add(kul_plex)
    else:
        logging.info(f"Plex {plex_id} already exists in the database")

    db.session.commit()


    # Populating the compound table with Kuljaninin compounds

    if comfile:
        # open the Kuljanin compounds file as dataframe
        kul_df = pd.read_csv(comfile, sep=',')

        compounds = db.session.query(Compound).all() # Query compound table 

        for _, row in tqdm(kul_df.iterrows(), total=kul_df.shape[0]):
            kul_compounds = row['Compound']
            kul_smiles = row['SMILES']

            # Changing the prefix of the compounds from CL or AC to KUL-
            if kul_compounds.startswith("CL"):
                kul_compounds = 'KULCL' + kul_compounds[2:]
            elif kul_compounds.startswith("AC"):
                kul_compounds = 'KULAC' + kul_compounds[2:]

            # Check if Kuljanin compounds already in table, if not add it to the compounds table
            found = False
            for c in compounds:
                if c.id == kul_compounds:
                    logging.info(f"Compound {kul_compounds} already in database")
                    found = True
                    break
            
            if not found:
                # Add the Kuljanin compounds to the compound table
                add_kulCompound = Compound(id=kul_compounds, smiles=kul_smiles)
                db.session.add(add_kulCompound)
        
        db.session.commit()
      

    # open the Kuljanin competition ratio file as dataframe
    df = pd.read_csv(crfile, sep=',')

    celltypes = set()


    for _, row in  tqdm(df.iterrows(), total=df.shape[0]):


        celltypes.add(row['Cell Type'])

    for celltype in celltypes:
        logging.info(f"Searching for cell type {celltype}")
        query = db.session.query(CellType).filter(CellType.name == celltype)
        if not db.session.query(query.exists()).scalar():
            logging.info(f"Adding cell type {celltype}")
            clt = CellType(name = celltype)
            db.session.add(clt)
    db.session.commit()


    added_entries = set()


    for _, row in  tqdm(df.iterrows(), total=df.shape[0]):

        uniprot = row['Uniprot ID']
        position = row['Site Position']
        cell_type = row['Cell Type']
        compound_name = row['Compound']
        comp_ratio = row['Competition Ratio']


        if compound_name.startswith("CL"):
             compound_name = 'KULCL' + compound_name[2:] #need to remove the '-' once the scriipt is working properly
        elif compound_name.startswith("AC"):
            compound_name = 'KULAC' + compound_name[2:]

        # Check that the Cystine residue is the residue table

        Cys_position = db.session.query(Residue).filter(Residue.uniprot == uniprot, Residue.position == position).first()
        if Cys_position is  None:
            logging.info(f"Uniprot{uniprot} in position {position} not in the Residue table")


        compound = db.session.query(Compound).filter(Compound.id == compound_name).first()
        probe = db.session.query(Compound).filter(Compound.id == 'DTB-IA').first()
        celltype = db.session.query(CellType).filter(CellType.name == cell_type).first()
        plex = db.session.query(Plex).filter(Plex.id == plex_id).first()

        samplename = f'{compound_name[3:]}{cell_type}'
        entry_key = (samplename, plex_id)

        if entry_key in added_entries:
            continue
        added_entries.add(entry_key)

        kul_cpdtrt = CompoundTreatment(
            plex = plex,
            samplename = samplename,
            compound = compound,
            probe = probe,
            celltype=celltype,
            tmtchannel = 'TMT11/16',
            concentration = 25,
            concentrationunits = 'uM',
            time = 2,
            timeunits = 'h',
            probeconcentration = 500,
            probeconcentrationunits = 'uM',
            probetime = 1,
            probetimeunits = 'h',
            comments = 'Nan',
            isreference = False
        )
        try:
            
            db.session.add(kul_cpdtrt)

        except Exception as err:
            logging.error(f"Error adding compound treatment {row}\n{err}")
            db.session.rollback()
        else:
            db.session.commit()            

    try:
        with open("progress.txt", "r") as file:
            last_processed_index = int(file.readline())
    except FileNotFoundError:
        last_processed_index = 0

    # Load competition ratio table, if connection dies last_processed_index will start from the last data position.
    for index, row in tqdm(df.iterrows(), total=df.shape[0], initial=last_processed_index):
        try:
            uniprot = row['Uniprot ID']
            position = row['Site Position']
            type = row['type']
            cell_type = row['Cell Type']
            compound_name = row['Compound']
            comp_ratio = row['Competition Ratio']
            if compound_name.startswith("CL"):
                compound_name = 'KULCL' + compound_name[2:]
            elif compound_name.startswith("AC"):
                compound_name = 'KULAC' + compound_name[2:]

            residue = db.session.query(Residue).filter(Residue.position == position, Residue.uniprot == uniprot, Residue.type == type).first()
            compoundt_treatment_id = db.session.query(CompoundTreatment).filter(CompoundTreatment.compound_id == compound_name).first()
            plex = db.session.query(Plex).filter(Plex.id == 'Kuljanin_plex').first()

            kul_cpr = CompetitionRatio(
                plex=plex,
                compoundtreatment_id=compoundt_treatment_id.id,
                residue=residue,
                scan=str(uuid.uuid4()),
                cr=comp_ratio,
                control_rsd=0,
                display_flag=True,
                multimapper=False,
                group_id=str(uuid.uuid4()),
                group_cr = comp_ratio
            )

            db.session.add(kul_cpr)
            db.session.commit()

            # Update the progress
            last_processed_index = index
            with open("progress.txt", "w") as file:
                file.write(str(last_processed_index))
        except Exception as err:
            logging.error(f"Error adding compound treatment {row}\n{err}")
            db.session.rollback()


@databp.cli.command("missed_residues")
@click.option('-mr', '--mrfile', required=True) 
@click.option('-l', '--logfile', default='missed_residues.log')
def external_compound_data(mrfile, logfile):
    "Command to check the percentage of cysteine resiudes not in the database due to outdated uniprot dataset used in the Kuljanin dataset"
    
    from models.chemopro import db, Residue
    
    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    df = pd.read_csv(mrfile, sep=',')
    

    unique_combinations = df[['Uniprot ID', 'Site Position']].drop_duplicates()

    uniq_cys = len(unique_combinations)
    
    logging.info(f"Number of Unique cysteine: {uniq_cys}")

    not_in_db = 0
    in_db = 0

    for _, row in tqdm(unique_combinations.iterrows(), total=unique_combinations.shape[0]):
        uniprot = row['Uniprot ID']
        position = row['Site Position']
    
        cysteine_check = db.session.query(Residue).filter(Residue.uniprot== uniprot, Residue.position == position).first()

        if cysteine_check is None:
            not_in_db +=1
        else:
            in_db += 1

    total_unique_combinations = unique_combinations.shape[0]

    percentage_not_in_db = (not_in_db / total_unique_combinations) * 100

    logging.info(f"Number of times cysteine is not in the database: {not_in_db}")
    logging.info(f"Number of times cysteine is in the database: {in_db}")

    logging.info(f"Percentage of times cysteine is not in the database: {percentage_not_in_db}%")



@databp.cli.command("gpprotein_data_import")
@click.option('-pf', '--proteinfile', required = True)
@click.option('-og', '--organism', required = True)
@click.option('-l', '--logfile', default = 'discoverome_protein_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def gpprotein_data_import(proteinfile, organism, logfile, drop_all):
    "Import protein synonyms and the position of all residues"

    from models.globalpro import db, GpProteinSynonym, GpProtein

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    if drop_all:
        print("***Dropping and creating all while testing***")
        db.create_all(bind_key= "globalproteomics")

    df = pd.read_csv(proteinfile, sep="\t", dtype={'Entry': str, 'Gene Names (primary)': str, 'Gene Names (synonym)': str}, keep_default_na=False)
    #df['Gene Names'] = df['Gene Names'].str.split(' ') # Split the values in the synonym column using the space between each values
    #df = df.explode('Gene Names') # Creates new row for each synonym retaining the same data as the original  
    df['Function [CC]'] = df['Function [CC]'].astype(str).map(lambda x:x.lstrip("FUNCTION:")) # Remove the 'FUNCTION:' string from function column
    df.rename(columns={'Entry':'uniprot', 'Gene Names (primary)': 'symbol', 'Gene Names (synonym)': 'synonyms', 'Function [CC]' : 'description'}, inplace= True)

    for _,row in tqdm(df.iterrows(), total=df.shape[0]):

        #check if protein already exists
        protein = db.session.query(GpProtein).filter(GpProtein.uniprot==row['uniprot']).first()

        #add Protein with symbol and description
        if protein:
            protein.symbol = row['symbol']
            protein.description = row['description']
            protein.organism = organism
            db.session.commit()
        else:
            protein = GpProtein(uniprot=row['uniprot'], symbol=row['symbol'], description=row['description'], organism = organism)
            db.session.add(protein)
            db.session.commit()

        #Go through each synonym and add to synonym table
        synonyms = row['synonyms'].split(' ')
        for synonym in synonyms:
            p_synonym = GpProteinSynonym(uniprot=row['uniprot'], synonym=synonym, type='uniprot_name')
            logging.info(f"Creating new protein synonym for {row['uniprot']} - {synonym}")
            db.session.add(p_synonym)
        db.session.commit()



@databp.cli.command("gpensembl_import")
@click.option('-f', '--file', required=True) 
@click.option('-l', '--logfile', default='ensembl_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def gpensembl_import(file, logfile, drop_all):
    "Import ensembl id for uniprot ids"
    
    from models.globalpro import db, GpProtein, GpProteinSynonym
    
    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    df = pd.read_csv(file, sep='\t')
    
    for _, row in tqdm(df.iterrows(), total=df.shape[0]):
        uniprot = row['uniprot']
        ensembl = row['ensembl']

        # Skip if 'ensembl' is empty or NaN
        if pd.isna(ensembl) or ensembl.strip() == "":
            logging.warning(f"Skipping row with uniprot '{uniprot}' due to missing Ensembl ID")
            continue

        data_type = 'ensembl_ID'
        
        # Check if the ProteinSynonym row already exists
        existing_synonym = db.session.query(GpProteinSynonym).filter_by(synonym=ensembl, type=data_type).first()
        if existing_synonym is not None:
            logging.info(f"Synonym '{ensembl}' already exists for uniprot '{existing_synonym.uniprot}'")
            continue
        
        # Add the ProteinSynonym row if it does not exist
        protein = db.session.query(GpProtein).filter_by(uniprot=uniprot).first()
        if protein is None:
            logging.warning(f"No protein found with uniprot '{uniprot}'")
            continue
        
        synonym = GpProteinSynonym(uniprot=uniprot, synonym=ensembl, type=data_type)
        db.session.add(synonym)
        db.session.commit()
        
        logging.info(f"Added synonym '{ensembl}' for uniprot '{uniprot}'")



@databp.cli.command("create_tables")
def create_tables():
    "creating all tables"

    from models.globalpro import db, FoldChange, GpIntensityReading, GpCompound, GpCellType, GpExperiment, GpPlex, GpCompoundTreatment

    print("***Creating all while testing***")
    db.create_all()

 

@databp.cli.command("abpp_create_tables")
def create_tables():
    "creating all tables"

    from models.chemopro import db, Protein, Residue, CompoundTreatment, Compound, \
    CellType, Experiment, IntensityReading, Plex, TargetList, ProteinToList, \
    ProteinSynonym, Compound, CompetitionRatio, Structure, StructureChain, StructureResidue, \
    Ligand, LigandResidueDistance, PocketResidue, Pocket, ResidueList, ResidueToList, ResidueFeature, \
    Compound_cr_four, Compound_cr_fifteen

    print("***Creating all while testing***")
    db.create_all()

@databp.cli.command("gp_experiment_import")
@click.option('-d', '--datafile', required=True)
@click.option('-s', '--sampletable', required=True)
@click.option('-og', '--organism', required = True)
@click.option('-l', '--logfile', default='discoverome_data_import.log')
@click.option('--drop-all/--no-drop-all', default=False)
def gp_experiment_import(datafile, sampletable, organism, logfile, drop_all):
    """Imports GP data file and sample sheet"""

    from models.chemopro import db, GpProtein, GpCompoundTreatment, GpCompound, GpCellType, GpExperiment, GpIntensityReading, GpPlex, FoldChange

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    logging.info(f"Importing {datafile} and {sampletable}")
    
    if not exists(datafile):
        raise FileNotFoundError

    if not exists(sampletable):
        raise FileNotFoundError

    if drop_all:
        print("***Dropping and creating all tables***")
        db.metadata.drop_all(db.engine)
        db.metadata.create_all(db.engine)
   
    #Read the sample table and import each sample as a compound treatment object with associated compounds and probes
    #Now find all the compounds (and probes) and cell types referenced in the file
    sample_table = pd.read_csv(sampletable)
    compounds = set()
    celltypes = set()
    samplenames = set()
    experiment_id = set()
    plex_id = set()
    print(sample_table.head())
    for _, row, in sample_table.iterrows():
        compounds.add(row['Compound'])
        celltypes.add(row['CellType'])
        samplenames.add(row['SampleName'])
        experiment_id.add(row['Experiment'])
        plex_id.add(row['Plex'])
    print(compounds)
    logging.info(f'list of cpd: {compounds}')
    #There should only be one experiment and one plex per datafile
    if len(experiment_id) > 1 or len(plex_id) > 1:
        logging.error(f"More than one experiment or plex ID in {sampletable}")
        raise Exception(f"More than one experiment or plex ID in {sampletable}")
    else:
        experiment_id = experiment_id.pop()
        plex_id = plex_id.pop()

    #Find the experiment or create it if it is new
    #If the plex already exist then throw an error
    experiment = db.session.query(GpExperiment).filter(GpExperiment.id == experiment_id)
    if not db.session.query(experiment.exists()).scalar():
        experiment = GpExperiment(id = experiment_id, description = "Experiment description")
        db.session.add(experiment)
        db.session.commit()
    else:
        experiment = experiment.first()

    plex = db.session.query(GpPlex).filter(GpPlex.id == plex_id)
    if not db.session.query(plex.exists()).scalar():
        plex = GpPlex(id = plex_id, gpexperiment = experiment, description = "Plex description")
        db.session.add(plex)
        db.session.commit()
    else:
        logging.error(f"Found plex: {plex_id} in db already.")
        raise Exception(f"Found plex: {plex_id} in db already.")

    #Create compound objects for each one as long as it doesn't exist in the db
    for compound_id in compounds:
        logging.info(f"Searching for {compound_id}")
        query = db.session.query(GpCompound).filter(GpCompound.id == compound_id)
        if not db.session.query(query.exists()).scalar():
            logging.info(f"Adding compound {compound_id}")
            cpd = GpCompound(id = compound_id)
            db.session.add(cpd)
    db.session.commit()

    for celltype in celltypes:
        logging.info(f"Searching for cell type {celltype}")
        query = db.session.query(GpCellType).filter(GpCellType.name == celltype)
        if not db.session.query(query.exists()).scalar():
            logging.info(f"Adding cell type {celltype}")
            clt = GpCellType(name = celltype)
            db.session.add(clt)
    db.session.commit()

    #Now create the compound treatment objects (also store these in a dictionary for quick look up later)
    cpdtrts = {}
    cpdtrt_groups = {}
    refs = []
    for _, row in sample_table.iterrows():
        logging.info(f"Adding compound treatment: {row['SampleName']}")
        compound = db.session.query(GpCompound).filter(GpCompound.id == row['Compound']).first()
        celltype = db.session.query(GpCellType).filter(GpCellType.name == row['CellType']).first()
            
        # Check if concentration_two exists in the row and is not NaN
        concentration_two = str(row['Concentration_two']) if 'Concentration_two' in row and not pd.isna(row['Concentration_two']) else None

        # Check if concentration_three exists in the row and is not NaN
        concentration_three = str(row['Concentration_three']) if 'Concentration_three' in row and not pd.isna(row['Concentration_three']) else None


        
        
        cpdtrt = GpCompoundTreatment(
            gpplex = plex,
            samplename = row['SampleName'],
            gpcompound = compound,
            gpcelltype = celltype,
            tmtchannel = row['TMTChannel'],
            concentration = str(row['Concentration']),
            concentration_two = concentration_two,
            concentration_three = concentration_three,
            concentrationunits = row['ConcentrationUnits'],
            data_acq = 'DIA' if pd.isna(row['TMTChannel']) or row['TMTChannel'] == '' else 'DDA',
            time = row['Time'],
            timeunits = row['TimeUnits'],
            temperature = row['Temperature'],
            comments = row['Comments'],
            isreference = row['IsReference'] == "Y",
            organism = organism
        )
        try:
            db.session.add(cpdtrt)
            db.session.flush()
            db.session.refresh(cpdtrt)
            cpdtrts[row['SampleName']] = cpdtrt.id
            cpdtrt_groups[row['SampleName']] = \
                row['Compound']+str(row['Concentration'])+row['ConcentrationUnits']+ \
                str(row['Time'])+row['TimeUnits']+\
                row['CellType'] #Same compound treatment and probe treatment means samples share a 'group'
            if row['IsReference'] == "Y":
                refs.append(row['SampleName'])
        except Exception as err:
            logging.error(f"Error adding compound treatment {row}\n{err}")
            db.session.rollback()
        else:
            db.session.commit()

    print(refs)
    #replace the compound ids with uuids that define group membership
    # Create a set of unique values in the dictionary
    unique_cpds = set(cpdtrt_groups.values())

    # Generate a uuid for each unique value and create a mapping between old values and new uuids
    mapping = {val: str(uuid.uuid4()) for val in unique_cpds}

    # Replace each value in the original dictionary with its corresponding uuid from the mapping
    for key, val in cpdtrt_groups.items():
        cpdtrt_groups[key] = mapping[val]

    #Read the data file and import each intensity reading linking it to each compound treatment    
    data = pd.read_csv(datafile, dtype={'Position': str}, keep_default_na=False)

    #We will save the mean reference intensity for each peptide/scan from each residue in this dictionary of dictionaries
    #To use at the end to flag which peptide/scan is to be displayed
    mean_ref_intensities = {}

    for _, row in tqdm(data.iterrows(), total=data.shape[0]):

        #Errors here should back out the plex and compound treatments. Otherwise you cannot reload

        logging.info(f"Searching for protein {row['UniProt']}")

        #Find the proteins
        uniprots = row['UniProt'].split(";")

        for uniprot in uniprots:

            protein = db.session.query(GpProtein).filter(GpProtein.uniprot == uniprot, GpProtein.organism == organism).first()
            if protein == None:
                logging.warning(f"Protein not found in db: {uniprot}")
                continue

            #Add an entry to the dictionary of residues if we've not before
            protein_id = uniprot
            if not protein_id in mean_ref_intensities:
                mean_ref_intensities[protein_id] = {}

            #Now go through each sample name from the sample table and add that intensity reading
            #Find reference samples
            if len(refs) == 0:
                logging.error("No references found - check if compound treatments added?")

            # ref_values = []
            # for ref_sample in refs:
            #     ref_values.append(row[ref_sample])
    
            # ref_value_mean = stat.mean(ref_values)

            # if ref_value_mean == 0:
            #     logging.info(f"All control channels are 0 - skipping")
            #     continue

            # ref_value_rsd = (stat.stdev(ref_values) / ref_value_mean) * 100
            
            # mean_ref_intensities[protein_id][row['Scan']] = ref_value_mean

            # logging.info(f"Reference values {ref_values}")
            # logging.info(f"Adding DMSO rsd {ref_value_rsd}")
             
            replicate_group= {}

            existing_entries = set()  # Set to store existing combinations of protein_id and scan to stop duplicate entries in foldchange table

            for sample in samplenames:
                value = row[sample]
                if value == "NA" or value == '':
                    value = float('nan')
                else:
                    value = float(value)
                if value == 0:
                    value = 1  # To avoid infinity

                group_id = cpdtrt_groups[sample]
                scan = row['Scan']

                if group_id not in replicate_group:
                    replicate_group[group_id] = {}
                    replicate_group[group_id][protein_id] = {}
                    replicate_group[group_id][protein_id][scan] = [value]
                else:
                    if protein_id not in replicate_group[group_id]:
                        replicate_group[group_id][protein_id] = {}
                        replicate_group[group_id][protein_id][scan] = [value]
                    else:
                        if scan not in replicate_group[group_id][protein_id]:
                            replicate_group[group_id][protein_id][scan] = [value]
                        else:
                            replicate_group[group_id][protein_id][scan].append(value)

                intensityreading = GpIntensityReading(
                    gpplex=plex,
                    compoundtreatment_id=cpdtrts[sample],
                    protein_id=protein_id,
                    scan=scan,
                    value=value
                )

                if 'DMSO' in sample:
                    continue

                # Check if the combination of protein_id and scan already exists
                combination = (protein_id, scan, cpdtrts[sample])
                if combination in existing_entries:
                    continue  # Skip adding if the combination already exists
                else:
                    existing_entries.add(combination)  # Add the new combination to existing entries set


                if 'no_peptides' in row:
                    no_peptides=row['no_peptides']
                    no_unique_peptides=row['no_unique_peptides']

                else:
                    no_peptides = None
                    no_unique_peptides = None

                try:
                    logfc_value = float(row[f'logFC_{sample[:-2]}'])
                    p_value = float(row[f'p.mod_{sample[:-2]}'])
                except (ValueError, TypeError):
                    # Skip rows with invalid values
                    logging.warning(f"Skipping row due to invalid values: {row}")
                    continue  # Move to the next row

                foldchange = FoldChange(
                    gpplex=plex,
                    compoundtreatment_id=cpdtrts[sample],
                    protein_id=protein_id,
                    scan=scan,
                    foldchange=logfc_value,
                    p_value=p_value,
                    # control_rsd =ref_value_rsd,
                    no_peptides = no_peptides,
                    no_unique_peptides = no_unique_peptides,
                    group_id=group_id
                )
                db.session.add(intensityreading)
                db.session.add(foldchange)
            try:
                db.session.commit()
            except:
                logging.info(f"Error adding {row} rolling back")
                db.session.rollback()
                raise Exception
            
            for group_id, group_data in replicate_group.items():
                for protein_id, scans_data in group_data.items():
                    for scan, values in scans_data.items():

                        foldchange_replicates = db.session.query(FoldChange).filter(
                            FoldChange.protein_id == protein_id,
                            FoldChange.group_id == group_id,
                            FoldChange.scan == scan
                        ).all()


                        for line in foldchange_replicates:
                            line.replicate_no = len(values)
                        db.session.commit()


@databp.cli.command("gp_delete_plex")
@click.option('-p', '--plexid', required=True)
@click.option('-l', '--logfile', default='discoverome_delete_plex.log')
@click.option('--dry-run/--no-dry-run', default=False)
def delete_plex(plexid, logfile, dry_run):
    """Imports discoverome structural biology information"""

    from models.globalpro import db, GpPlex, FoldChange, GpIntensityReading, GpCompoundTreatment

    logging.basicConfig(
        filename=logfile,
        format='%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    logging.info(f"Attempting to delete plex '{plexid}'")

    plex = db.session.query(GpPlex).filter(GpPlex.id == plexid).first()
    if not plex:
        logging.info(f"Could not find plex '{plexid}' to delete - aborting")
        exit()
    
    crs = db.session.query(FoldChange).filter(FoldChange.plex_id == plexid).all()
    logging.info(f"Found {len(crs)} CRs to delete from '{plexid}'")

    irs = db.session.query(GpIntensityReading).filter(GpIntensityReading.plex_id == plexid).all()
    logging.info(f"Found {len(irs)} IRs to delete from '{plexid}'")

    cts = db.session.query(GpCompoundTreatment).filter(GpCompoundTreatment.plex_id == plexid).all()
    logging.info(f"Found {len(cts)} compound treatments to delete from '{plexid}'")

    logging.info(f"Deleting CRs")
    db.session.query(FoldChange).filter(FoldChange.plex_id == plexid).delete()
    logging.info(f"Deleting IRs")
    db.session.query(GpIntensityReading).filter(GpIntensityReading.plex_id == plexid).delete()
    logging.info(f"Deleting compound treatments")
    db.session.query(GpCompoundTreatment).filter(GpCompoundTreatment.plex_id == plexid).delete()
    logging.info(f"Deleting plex '{plexid}'")
    db.session.query(GpPlex).filter(GpPlex.id == plexid).delete()

    if not dry_run:
        logging.info("Committing changes")
        db.session.commit()
    else:
        logging.info("Dry run - no deletions made")




@databp.cli.command("gp_smiles_import")
@click.option('-sf', '--smilefile', required = True)
@click.option('-l', '--logfile', default = 'smiles_import.log')
def smiles_import(smilefile, logfile):
    "Import smiles string to compound table"

    from models.globalpro import db, GpCompoundTreatment, GpCompound

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    df_smiles = pd.read_csv(smilefile, sep=",") # Load file in as a pandas dataframe
    df_smiles= df_smiles.filter(['ID','SMILES', 'Batch ID', "Formatted ID", 'Project'], axis = 1) # Filter the table to have only 'ID' and 'Smiles' columns 

    for _,row in tqdm(df_smiles.iterrows(), total=df_smiles.shape[0]): # For the loading function 
        compound_id = row["ID"] # get the compound id's 
        smiles = row["SMILES"] # get the smiles string 
        batch_id = row["Batch ID"]
        formatted_id = row["Formatted ID"]
        project = row["Project"] if pd.notna(row["Project"]) else "library"

        compound_smiles = db.session.query(GpCompound).filter(or_(GpCompound.id == compound_id, GpCompound.id.like(f'%{batch_id}%'), GpCompound.id.like(f'%{formatted_id}%'))) # query the database to see if 'id' in table matched 'id' in dataframe
        if compound_smiles.count() == 0: # if not matches, entry into the log file
            logging.info(f"Compound{compound_id} does not exist in the compound table ")

        else: # If compound 'id' matches, add the smiles string to the database compound table
            logging.info(f"Updating compound {compound_id} with SMILES string {smiles}")
            for c in compound_smiles:
                c.smiles = smiles
                db.session.commit()

        compound_projects = (
            db.session.query(GpCompoundTreatment)
            .join(GpCompound, GpCompoundTreatment.compound_id == GpCompound.id)
            .filter(or_(GpCompoundTreatment.compound_id.like(f"%{batch_id}%"), GpCompound.id.like(f"%{formatted_id}%")))
            .all()
        )        

        if not compound_projects:
            logging.info(f"Compound {batch_id} does not exist in the CompoundTreatment table")
        else:
            logging.info(f"Updating compound {batch_id} with project id {project} for {len(compound_projects)} records")
            for c in compound_projects:
                c.project = project
            db.session.commit()







@databp.cli.command("gp_cpd_images")
@click.option('-l', '--logfile', default = 'cpd_images.log')
def smiles_import(logfile):
    "Calculate images from SMILES, save into folder and import path to compound table"
    
    from models.globalpro import db, GpCompound

    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )

    compounds = db.session.query(GpCompound).all() # Query compound table 

    for s in compounds: # loop through each compound in the table 

        if s.image is None and s.smiles is not None:
            logging.info(f"For Compound {s.id} creating new SMILES strucutre {s.smiles}")            
            mol = Chem.MolFromSmiles(s.smiles) # create a chemical strucutre from the smiles string
            img = Draw.MolToImage(mol) #draw the strucutre 
            file = f'/resources/smiles-structure/{s.id}.png'
            if os.path.exists(file):
                logging.info(f"File for {s.id} already exists")
            else:
                img.save(file)

            s.image = file
            logging.info(f"Adding image path of {s.id} to image field")
            db.session.commit()

        else:
            if s.image is not None:
                logging.info(f"Image path already exists for {s.id}")
            else:
                logging.info(f"No smiles entry for {s.id}")






@databp.cli.command("pathway_import")
@click.option('-f', '--file', required = True) 
@click.option('-db', '--dbname', type=str, required = True) # List name is required unlike description
@click.option('-og', '--organism', type=str, required = True) # List name is required unlike description
@click.option('-d', '--description', type =str )
@click.option('-l', '--logfile', default = 'pathway_list_import.log')
def pathway_import(file, dbname, organism, description, logfile): # Function takes in file, name of list and description of list
    "Import pathway list data"

    from models.globalpro import db, PathwayList, GpProtein, ProteinToPathway
    
    logging.basicConfig(
        filename = logfile,
        format = '%(asctime)s %(message)s',
        level=logging.DEBUG
    )
    

    pathway_table = pd.read_csv(file)

    for _, row in pathway_table.iterrows():
        pathway_name = row['pathways']
        protein_list = [sym.strip() for sym in row['symbol'].split(',')]  # Ensure all proteins are included and stripped of whitespace

        # Check if the pathway already exists
        existing_pathway = db.session.query(PathwayList).filter_by(name=pathway_name, database=dbname, organism=organism).first()

        if not existing_pathway:
            new_list = PathwayList(
                name=pathway_name,
                database=dbname,
                description=description,
                organism=organism
            )
            db.session.add(new_list)
            db.session.flush()  # Ensure new_list.id is available

            logging.info(f"Created new pathway - {pathway_name} ({new_list.id})")
        else:
            new_list = existing_pathway
            logging.info(f"Using existing pathway - {pathway_name} ({new_list.id}`  )")

        # Link proteins to the pathway
        for sym in protein_list:
            protein = db.session.query(GpProtein).filter(GpProtein.symbol == sym, GpProtein.organism == organism).first()

            if protein:
                # Check if the association already exists
                existing_link = db.session.query(ProteinToPathway).filter_by(protein_id=protein.uniprot, pathway_id=new_list.id).first()

                if not existing_link:
                    link = ProteinToPathway(protein_id=protein.uniprot, pathway_id=new_list.id)
                    db.session.add(link)
                    logging.info(f"Linked {sym} to pathway {pathway_name}")
                else:
                    logging.info(f"Link already exists for {sym} in {pathway_name}")
            else:
                logging.warning(f"Protein {sym} not found in database")

    db.session.commit()
    logging.info("Pathway import completed successfully")

@databp.cli.command("remove_compound_data")
@click.option( "-cpid", "--compoundid", required=True, help="Compound ID pattern")
@click.option( "-l", "--logfile", default="remove_compound.log")
@click.option( "--dry-run", is_flag=True, help="Preview deletions without committing")
def remove_compound_data(compoundid, logfile, dry_run):

    """
    Cascade remove experiment data for compound IDs matching a pattern.
    """

    from models.globalpro import ( db, Compound, CompoundTreatment, CompetitionRatio, IntensityReading, Compound_cr_fifteen, Compound_cr_four)

    logging.basicConfig(
        filename=logfile,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s"
    )

    logging.info(f"Starting deletion for pattern: {compoundid}")

    try:

        # Find matching compound IDs
        compounds = Compound.query.filter(
            Compound.id.like(f"%{compoundid}%")
        ).all()

        logging.info(f"Found {len(compounds)} compounds")

        for compound in compounds:

            logging.info(f"Processing compound: {compound.id}")

            # Find treatments
            treatments = CompoundTreatment.query.filter_by(
                compound_id=compound.id
            ).all()

            treatment_ids = [t.id for t in treatments]

            logging.info(
                f"Found {len(treatment_ids)} compound treatments"
            )

            if treatment_ids:

                # Delete Competition Ratios
                cr_deleted = CompetitionRatio.query.filter(
                    CompetitionRatio.compoundtreatment_id.in_(treatment_ids)
                ).delete(synchronize_session=False)

                logging.info(
                    f"Deleted {cr_deleted} CompetitionRatio rows"
                )

                # Delete Intensity Readings
                ir_deleted = IntensityReading.query.filter(
                    IntensityReading.compoundtreatment_id.in_(treatment_ids)
                ).delete(synchronize_session=False)

                logging.info(
                    f"Deleted {ir_deleted} IntensityReading rows"
                )

                # Delete Compound Treatments
                ct_deleted = CompoundTreatment.query.filter(
                    CompoundTreatment.id.in_(treatment_ids)
                ).delete(synchronize_session=False)

                logging.info(
                    f"Deleted {ct_deleted} CompoundTreatment rows"
                )

            

            # Delete compound itself
            db.session.delete(compound)

            logging.info(f"Deleted compound row: {compound.id}")

        # Commit or rollback
        if dry_run:

            db.session.rollback()

            logging.info(
                "DRY RUN enabled - rollback performed"
            )

            print("Dry run completed. No changes committed.")

        else:

            db.session.commit()

            logging.info("Database commit successful")

            print("Deletion completed successfully.")

    except Exception as e:

        db.session.rollback()

        logging.exception("Deletion failed")

        print(f"Error: {str(e)}")