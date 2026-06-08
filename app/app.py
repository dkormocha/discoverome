from inspect import iscoroutinefunction
from re import I
from flask import render_template, request, redirect, url_for, g, jsonify, flash, send_file, send_from_directory, session
from sqlalchemy.sql import func, or_
from sqlalchemy import desc, nullslast, not_, asc
from sqlalchemy.orm import joinedload

import numpy as np
import re
import pandas as pd
import json
import math
import plotly.graph_objects as go
import os
import requests as http_requests
import ollama

import plotly.express as px

import plotly.io as pio
from plotly.io import to_html
import plotly

import gseapy as gp
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from app import create_app
from app.commands import databp

from app.extensions import db
from app.models.models import *

from gseapy import enrichment_map
import networkx as nx
import matplotlib as mpl
from matplotlib.colors import to_rgba


app = create_app()


@app.route('/autocomplete')
def autocomplete():
    term = request.args.get('term')
    matching_proteins = db.session.query(Protein).filter(
        Protein.symbol.ilike(f'%{term}%')
    ).limit(20).all()
    matching_uniprots = db.session.query(Protein).filter(
        Protein.uniprot.ilike(f'%{term}%')
    ).limit(20).all()
    matching_synonyms = db.session.query(ProteinSynonym).filter(
        ProteinSynonym.synonym.ilike(f'%{term}%')
    ).limit(20).all()
    matching_compounds = db.session.query(Compound).filter(
        Compound.id.ilike(f'%{term}%')
    ).limit(20).all()
    matching_items = [protein.symbol for protein in matching_proteins] + [synonym.synonym for synonym in matching_synonyms] + [compound.id for compound in matching_compounds] + [ uniprot.uniprot for uniprot in matching_uniprots]

    return jsonify(matching_items)






@app.route('/', methods=['POST','GET'])
def index():
    

    if request.method == 'POST':
        symbol = request.form['sym']
        symbol = symbol.upper()
        # check if input is a protein symbol
        protein = db.session.query(Protein).filter(Protein.symbol == symbol).first()
        compound = db.session.query(Compound).filter(Compound.id == symbol).first()
        synonyms = db.session.query(ProteinSynonym).filter(ProteinSynonym.synonym == symbol).all()

        if protein is not None:
            return redirect(url_for('protein_search', symbol=protein.symbol))
        # check if input is a synonym
        if len(synonyms) > 0:
            protein = synonyms[0].protein
            return redirect(url_for('protein_search', symbol=protein.symbol))
        # check if input is a uniprot ID
        protein = db.session.query(Protein).filter(Protein.uniprot == symbol).first()
        if protein is not None:
            return redirect(url_for('protein_search', symbol=protein.symbol))
        elif compound is not None:
            return redirect(url_for('compound_search', name=compound.id))
        else:
            flash('No results found. Please check your search query and try again.', 'error')


    


    return render_template('index.html')





@app.route('/advance_search', methods=['GET'])
def advance_search():
    # --- Parse parameters ---
    params = {
        'cr_min': request.args.get('cr_min'),
        'cr_max': request.args.get('cr_max'),
        'pocket_vol_min': request.args.get('pocket_vol_min'),
        'pocket_vol_max': request.args.get('pocket_vol_max'),
        'drug_score_min': request.args.get('drug_score_min'),
        'drug_score_max': request.args.get('drug_score_max'),
        'distance_min': request.args.get('distance_min'),
        'distance_max': request.args.get('distance_max'),
        'access_min': request.args.get('access_min'),
        'access_max': request.args.get('access_max'),
    }

    # Promiscuity needs special int conversion
    promis_min = promis_max = None
    try:
        promis_min = int(float(request.args.get('promis_min')))
        promis_max = int(float(request.args.get('promis_max')))
    except (TypeError, ValueError):
        pass

    def is_set(val):
        return val is not None and val != 'None'

    def to_float(val):
        try:
            return float(val)
        except (TypeError, ValueError):
            return None

    # --- Build protein filter ---
    proteinlist = request.args.get('proteinlist')
    res_type = request.args.get('res_type')
    res_type = res_type if res_type in ('CYS', 'LYS') else 'CYS'

    protein_filter = None  # None means "no protein restriction — search all"

    if proteinlist and proteinlist.strip():
        # User typed in specific UniProt IDs
        proteins = [p.upper() for p in re.split(r',\s*|\s+', proteinlist) if p.strip()]
        if proteins:
            protein_filter = Protein.uniprot.in_(proteins)




    # Only add to base query filters if a protein filter was actually built
    base_filters = [
        Residue.type == res_type,
        CompoundTreatment.compound_id != 'DMSO',
        CompetitionRatio.control_rsd <= 30,
        CompetitionRatio.display_flag == True,
        CompetitionRatio.multimapper == False,
        CompetitionRatio.group_cr != 'NaN',
        CompetitionRatio.group_cr >= 4

    ]

    if protein_filter is not None:
        base_filters.append(protein_filter)


    # --- Base query (shared joins) ---
    base_query = (
        db.session.query(
            CompetitionRatio, Residue, Protein, CompoundTreatment,
            CellType, StructureResidue, Pocket, PocketResidue, StructureChain
        )
        .join(Residue, CompetitionRatio.residue_id == Residue.id)
        .join(Protein, Residue.uniprot == Protein.uniprot)
        .join(CompoundTreatment, CompetitionRatio.compoundtreatment_id == CompoundTreatment.id)
        .join(CellType, CompoundTreatment.celltype_id == CellType.id)
        .join(Compound_cr_four, CompoundTreatment.compound_id == Compound_cr_four.compound_id)
        .outerjoin(StructureResidue, Residue.id == StructureResidue.residue_id)
        .outerjoin(PocketResidue, StructureResidue.id == PocketResidue.structureresidue_id)
        .outerjoin(Pocket, PocketResidue.pocket_id == Pocket.id)
        .outerjoin(StructureChain, StructureResidue.chain_id == StructureChain.id)
        .outerjoin(LigandResidueDistance, StructureResidue.id == LigandResidueDistance.structureresidue_id)
        .outerjoin(Ligand, LigandResidueDistance.ligand_id == Ligand.id)
        .filter(*base_filters)          # ← use base_filters here, not the old inline list
        .with_entities(
            CompoundTreatment.compound_id, Compound_cr_four.count,
            Protein.uniprot, Protein.symbol, Residue.position,
            CellType.name, CompetitionRatio.group_cr,
            Pocket.structure_id, StructureResidue.in_disulfide,
            StructureResidue.accessibility, Pocket.pocket_volume_MC,
            Pocket.drug_score, Ligand.code,
            Ligand.structure_id.label('ligand_structure'),
            Ligand.mw, LigandResidueDistance.distance, StructureResidue.depth
        )
    )



    # --- Dynamically build optional filters ---
    filters = []

    cr_min_v   = to_float(params['cr_min'])
    cr_max_v   = to_float(params['cr_max'])
    pv_min_v   = to_float(params['pocket_vol_min'])
    pv_max_v   = to_float(params['pocket_vol_max'])
    ds_min_v   = to_float(params['drug_score_min'])
    ds_max_v   = to_float(params['drug_score_max'])
    dist_min_v = to_float(params['distance_min'])
    dist_max_v = to_float(params['distance_max'])
    acc_min_v  = to_float(params['access_min'])
    acc_max_v  = to_float(params['access_max'])

    # Competition ratio
    if cr_min_v is not None:
        filters.append(CompetitionRatio.group_cr >= cr_min_v)
    if cr_max_v is not None:
        filters.append(CompetitionRatio.group_cr <= cr_max_v)

    # Pocket volume
    if pv_min_v is not None:
        filters.append(Pocket.pocket_volume_MC >= pv_min_v)
    if pv_max_v is not None:
        filters.append(Pocket.pocket_volume_MC <= pv_max_v)

    # Drug score — NULL rows are treated as "passing" (your original logic)
    if ds_min_v is not None:
        filters.append(or_(Pocket.drug_score.is_(None), Pocket.drug_score >= ds_min_v))
    if ds_max_v is not None:
        filters.append(or_(Pocket.drug_score.is_(None), Pocket.drug_score <= ds_max_v))

    # Accessibility
    if acc_min_v is not None:
        filters.append(StructureResidue.accessibility >= acc_min_v)
    if acc_max_v is not None:
        filters.append(StructureResidue.accessibility <= acc_max_v)

    # Ligand distance
    if dist_min_v is not None:
        filters.append(LigandResidueDistance.distance >= dist_min_v)
    if dist_max_v is not None:
        filters.append(LigandResidueDistance.distance <= dist_max_v)

    # Promiscuity
    if promis_min is not None:
        filters.append(Compound_cr_four.count >= promis_min)
    if promis_max is not None:
        filters.append(Compound_cr_four.count <= promis_max)

    # Pocket confidence — apply whenever any pocket filter is active
    pocket_filters_active = any(
        v is not None for v in [pv_min_v, pv_max_v, ds_min_v, ds_max_v]
    )
    if pocket_filters_active:
        filters.append(or_(Pocket.min_confidence.is_(None), Pocket.min_confidence >= 70))

    # Ligand MW guard — apply whenever distance filter is active
    if dist_min_v is not None or dist_max_v is not None:
        filters.append(or_(Ligand.mw.is_(None), Ligand.mw >= 300))

    # Apply all collected filters in one go
    results = (
        base_query
        .filter(*filters)
        .distinct()
        .limit(20000)
        .all()
    )

    print(results)


    def _fmt(value, decimals=2):
        if value is None:
            return ''
        try:
            return round(value, decimals)
        except TypeError:
            return ''  # handles unexpected types gracefully

    rows = []
    for s in results:
        rows.append({
            'compound':         s.compound_id or '',
            'promiscuity':      s.count or '',
            'uniprot':          s.uniprot or '',
            'symbol':           s.symbol or '',
            'position':         s.position or '',
            'cell':             s.name or '',
            'cr':               _fmt(s.group_cr),
            'pocket_structure': s.structure_id or '',
            'disulfide':        bool(s.in_disulfide),
            'accessibility':    _fmt(s.accessibility) if s.accessibility is not None else 0,
            'depth':            _fmt(s.depth) if s.depth is not None else 0,
            'pocket_volume':    _fmt(s.pocket_volume_MC),
            'drug_score':       _fmt(s.drug_score),
            'ligand_structure': s.ligand_structure or '',
            'ligand_code':      s.code or '',
            'distance':         _fmt(s.distance),
            'ligand_mw':        _fmt(s.mw),
        })
    for i, row in enumerate(rows):
        if len(row) != 17:
            print(f"Row {i} has {len(row)} keys: {row}")

    print(rows)

    return render_template('advance_search.html', rows=rows, res_type=res_type)



@app.route('/images/<path:filename>')
def get_image(filename):
    project_root = os.path.dirname(app.root_path)
    return send_from_directory(project_root, filename)



@app.route('/pqr_files/<path:subpath>/<string:pocket_structure>')
def serve_pqr_file(subpath, pocket_structure):
    if len(pocket_structure) == 4:
        base_dir = '/resources/structural-biology/data/fpocket/mmCIF'
    else:
        base_dir = '/resources/structural-biology/data/fpocket/AF2'
    return send_from_directory(base_dir, subpath, as_attachment=True)


@app.route('/protein/<symbol>')
def protein_search(symbol):
    from collections import defaultdict

    celltype     = request.args.get('celltype')
    show_all     = bool(request.args.get('all_res'))
    residue_type = request.args.get('residue_type') if request.args.get('residue_type') in ('CYS', 'LYS') else 'CYS'

    protein = db.session.query(Protein).filter(Protein.symbol == symbol).first()
    if protein is None:
        return redirect(url_for('index'))

    celltypes    = db.session.query(CellType).all()
    ps           = list({s.synonym for s in db.session.query(ProteinSynonym).filter(ProteinSynonym.uniprot == protein.uniprot).all()})
    description  = protein.description
    protein_name = f'{protein.uniprot[1:3]}/{protein.uniprot}/{protein.uniprot}_pockets.pqr/{protein.uniprot}'

    residues = db.session.query(Residue).filter(
        Residue.uniprot == protein.uniprot, Residue.type == residue_type
    ).order_by(Residue.position).all()

    if not residues:
        return render_template('protein_search.html', protein=protein, description=description,
            ps=ps, data_labels={}, data={}, celltypes=celltypes, selected_residue_type=residue_type,
            selected_celltype=celltype, position_data={}, protein_structure={},
            protein_name=protein_name, btn_position={})

    residue_ids = [r.id for r in residues]

    # 1. Pocket structure list
    pocket_list = [p.Pocket.structure_id for p in (
        db.session.query(Pocket, PocketResidue, Residue, StructureResidue)
        .join(PocketResidue, Pocket.id == PocketResidue.pocket_id)
        .join(StructureResidue, PocketResidue.structureresidue_id == StructureResidue.id)
        .join(Residue, StructureResidue.residue_id == Residue.id)
        .filter(Residue.uniprot == protein.uniprot, Residue.type == residue_type)
        .distinct(Pocket.structure_id).all()
    )]

    # 2. Promiscuity maps
    cr4_map  = {r.compound_id: r.count for r in db.session.query(Compound_cr_four).all()}
    cr15_map = {r.compound_id: r.count for r in db.session.query(Compound_cr_fifteen).all()}

    # 3. All CRs for this protein in one query
    cr_q = (
        db.session.query(CompetitionRatio, CompoundTreatment, CellType, Compound, Plex)
        .filter(
            CompetitionRatio.residue_id.in_(residue_ids),
            CompetitionRatio.display_flag == True,
            CompetitionRatio.multimapper == False,
            CompetitionRatio.control_rsd <= 30,
            CompoundTreatment.id == CompetitionRatio.compoundtreatment_id,
            CompoundTreatment.isreference == False,
            CellType.id == CompoundTreatment.celltype_id,
            Compound.id == CompoundTreatment.compound_id,
            Plex.id == CompoundTreatment.plex_id,
        )
    )
    if celltype:
        cr_q = cr_q.filter(CellType.name == celltype)
    cr_rows = cr_q.all()

    cr_by_residue = defaultdict(list)
    for row in cr_rows:
        cr_by_residue[row.CompetitionRatio.residue_id].append(row)

    # 4. Best pocket per residue
    pocket_by_residue = {}
    for row in (
        db.session.query(Residue, PocketResidue, Pocket, StructureResidue, Ligand, LigandResidueDistance)
        .select_from(PocketResidue)
        .join(Pocket)
        .join(StructureResidue, PocketResidue.structureresidue_id == StructureResidue.id)
        .join(Residue, StructureResidue.residue_id == Residue.id)
        .outerjoin(LigandResidueDistance, StructureResidue.id == LigandResidueDistance.structureresidue_id)
        .outerjoin(Ligand, LigandResidueDistance.ligand_id == Ligand.id)
        .filter(Residue.id.in_(residue_ids), Residue.type == residue_type)
        .order_by(Pocket.pocket_volume_MC.desc())
        .all()
    ):
        rid = row[0].id
        if rid not in pocket_by_residue:
            pocket_by_residue[rid] = row

    # 5. Best ligand per residue (no pocket)
    ligand_by_residue = {}
    for row in (
        db.session.query(Residue, StructureResidue, Ligand, LigandResidueDistance)
        .filter(
            Residue.id.in_(residue_ids),
            Residue.type == residue_type,
            StructureResidue.residue_id == Residue.id,
            StructureResidue.id == LigandResidueDistance.structureresidue_id,
            LigandResidueDistance.ligand_id == Ligand.id,
        )
        .order_by(Ligand.mw.desc())
        .all()
    ):
        rid = row[0].id
        if rid not in ligand_by_residue:
            ligand_by_residue[rid] = row

    # 6. Residue features
    features_by_residue = defaultdict(list)
    for f in db.session.query(ResidueFeature).filter(ResidueFeature.residue_id.in_(residue_ids)).all():
        features_by_residue[f.residue_id].append(f)

    # 7. Compound cache and btn_position
    all_compound_ids = {row.CompoundTreatment.compound_id for row in cr_rows}
    compound_cache = {c.id: c for c in db.session.query(Compound).filter(
        Compound.id.in_(all_compound_ids)
    ).all()} if all_compound_ids else {}

    protein_structure = {p: (p, f'{p[1:3]}/{p}/{p}_pockets.pqr') for p in pocket_list}

    btn_position = {}
    if pocket_list:
        btn_map = defaultdict(set)
        for sr, res, sc in (
            db.session.query(StructureResidue, Residue, StructureChain)
            .join(Residue, StructureResidue.residue_id == Residue.id)
            .join(StructureChain, StructureResidue.chain_id == StructureChain.id)
            .filter(
                StructureResidue.structure_id.in_(pocket_list),
                Residue.uniprot == protein.uniprot,
                Residue.type == residue_type,
            )
            .all()
        ):
            btn_map[sr.structure_id].add((res.position, int(sr.pdb_position), sc.chain))
        btn_position = {p: sorted(btn_map[p]) for p in pocket_list}

    def parse_features(features):
        src_list, des_list = [], []
        for rf in features:
            m = re.search(r'UniProt (.+)', rf.source)
            if m:
                src_list.append(m.group(1))
            evidence = [s[11:-1] for s in re.findall(r'/evidence=".+?"', rf.description)]
            note     = [s[6:-1]  for s in re.findall(r'/note=".+?"',     rf.description)]
            li       = re.findall(r'/ligand="(.+?)";\s*/ligand_id="(.+?)";', rf.description)
            li_str   = "; ".join(f"{l[0].upper()}; {l[1]}" for l in li)
            des_list.append('; '.join(note + evidence + [li_str]))
        return (
            '~'.join(src_list) if src_list else "N/A",
            '~'.join(des_list) if des_list else "N/A",
        )

    # Main loop: zero DB queries inside
    data, data_labels, position_data = {}, {}, {}
    all_compounds = set()

    for residue in residues:
        position     = residue.position
        comp_val     = cr_by_residue.get(residue.id, [])
        pocket_row   = pocket_by_residue.get(residue.id)
        ligand_row   = ligand_by_residue.get(residue.id)
        res_features = features_by_residue.get(residue.id, [])

        res_sources, res_des = parse_features(res_features) if res_features else (None, None)

        if pocket_row:
            pkt, sres, lig, lig_dist_row = pocket_row[2], pocket_row[3], pocket_row[4], pocket_row[5]
            pocket_structure = pkt.structure_id
            p_score, d_score = pkt.pocket_score, pkt.drug_score
            pocket_vol  = round(pkt.pocket_volume_MC,   2) if pkt.pocket_volume_MC   else None
            mean_con    = round(pkt.mean_confidence,    2) if pkt.mean_confidence    else None
            medium_con  = round(pkt.median_confidence,  2) if pkt.median_confidence  else None
            min_con     = round(pkt.min_confidence,     2) if pkt.min_confidence     else None
            di_sulfide  = sres.in_disulfide
            accessibility = round(sres.accessibility,  2) if sres.accessibility     else None
            ligand_dist = round(lig_dist_row.distance,  2) if lig_dist_row           else None
            lig_struId  = lig.structure_id if lig else None
            lig_name    = lig.name         if lig else None
            lig_mw      = round(lig.mw, 2) if lig and lig.mw else None
            artefact    = lig.artefact     if lig else None
            has_ld = lig_dist_row is not None
            if has_ld and res_sources:   position_data[position] = 'Pocket Ligand and Residue Features'
            elif has_ld:                 position_data[position] = 'Pocket Ligand'
            elif res_sources:            position_data[position] = 'Pocket and Residue Features'
            else:                        position_data[position] = 'Pocket'
        else:
            pocket_structure = p_score = d_score = pocket_vol = None
            mean_con = medium_con = min_con = di_sulfide = accessibility = None
            lig_struId = ligand_dist = lig_name = lig_mw = artefact = None
            if ligand_row:
                rl = ligand_row
                lig_struId  = rl[2].structure_id
                ligand_dist = round(rl[3].distance, 2) if rl[3].distance else None
                lig_name    = rl[2].name
                lig_mw      = round(rl[2].mw, 2)      if rl[2].mw      else None
                artefact    = rl[2].artefact
                position_data[position] = 'Ligand and Residue Features' if res_sources else 'Ligand'
            elif res_sources:
                position_data[position] = 'Residue Features'

        struct_tuple = (
            pocket_structure, p_score, d_score, pocket_vol,
            mean_con, medium_con, min_con, di_sulfide, accessibility,
            lig_struId, ligand_dist, lig_name, lig_mw, artefact,
            res_sources, res_des,
        )

        for cr in comp_val:
            compound     = cr.CompoundTreatment.compound_id
            all_compounds.add(compound)
            data.setdefault(position, {})[compound] = round(cr.CompetitionRatio.group_cr, 2)
            data_labels.setdefault(position, {})[compound] = (
                protein.uniprot, cr.CellType.name,
                cr.Compound.smiles, cr.Compound.image,
                *struct_tuple,
                cr4_map.get(compound), cr15_map.get(compound),
                cr.CompetitionRatio.p_value or '',
                cr.CompetitionRatio.replicate_no or '',
            )

        if show_all:
            data.setdefault(position, {})
            data_labels.setdefault(position, {})
            for compound in (all_compounds or ['Compound']):
                if compound not in data[position]:
                    comp = compound_cache.get(compound)
                    data[position][compound] = None
                    data_labels[position][compound] = (
                        protein.uniprot, None,
                        comp.smiles if comp else None,
                        comp.image  if comp else None,
                        *struct_tuple,
                        None, None, None, None,
                    )

    if not show_all:
        for position in list(data):
            existing = next(iter(data_labels.get(position, {}).values()), None)
            for compound in all_compounds:
                if compound not in data[position]:
                    comp = compound_cache.get(compound)
                    data[position][compound] = None
                    data_labels[position][compound] = (
                        protein.uniprot, None,
                        comp.smiles if comp else None,
                        comp.image  if comp else None,
                        *(existing[4:20] if existing and len(existing) >= 20 else ()),
                    )

    return render_template('protein_search.html', protein=protein, description=description,
        ps=ps, data_labels=data_labels, data=data, celltypes=celltypes,
        selected_residue_type=residue_type, selected_celltype=celltype,
        position_data=position_data, protein_structure=protein_structure,
        protein_name=protein_name, btn_position=btn_position)



 
@app.route('/compound/<name>')
def compound_search(name):

    data = []
    max_cr_dict = {}
    compound = db.session.query(Compound).filter(Compound.id == name).first()
    celltypes = db.session.query(CellType).all()

    celltype = request.args.get('celltype')


    selected_residue_type = request.args.get('residue_type')  # Convert to uppercase for consistency
    
    # Check if the user selected 'CYS' or 'LYS', otherwise set the default to 'CYS'
    if selected_residue_type in ['CYS', 'LYS']:
        residue_type = selected_residue_type
    else:
        residue_type = 'CYS'  # Default to 'CYS' if the selection is not valid


        


    if compound is None:
        return redirect(url_for('index'))
    
    else:       
        comp_image = compound.image
        comp_smiles = compound.smiles
        
        if celltype:
            comp_proteins = db.session.query(CompoundTreatment, CompetitionRatio, Residue, Protein, CellType).filter(
                CompoundTreatment.compound_id == name
            ).filter(
                CellType.id == CompoundTreatment.celltype_id, CellType.name == celltype
            ).filter(
                CompetitionRatio.compoundtreatment_id == CompoundTreatment.id, 
                CompetitionRatio.residue_id == Residue.id, CompoundTreatment.isreference == False
            ).filter(
                Residue.type == residue_type,
                CompetitionRatio.control_rsd <= 30,
                CompetitionRatio.display_flag == True,
                CompetitionRatio.multimapper == False
            ).filter(
                Protein.uniprot == Residue.uniprot
            ).order_by(func.random())
        else:
            comp_proteins = db.session.query(CompoundTreatment, CompetitionRatio, Residue, Protein, CellType).filter(
                CompoundTreatment.compound_id == name
            ).filter( CellType.id == CompoundTreatment.celltype_id,
                CompetitionRatio.compoundtreatment_id == CompoundTreatment.id, 
                CompetitionRatio.residue_id == Residue.id, CompoundTreatment.isreference == False
            ).filter(
                Residue.type == residue_type,
                CompetitionRatio.control_rsd <= 30,
                CompetitionRatio.display_flag == True,
                CompetitionRatio.multimapper == False


            ).filter(
                Protein.uniprot == Residue.uniprot
            ).order_by(func.random())


        for comp in comp_proteins: 
            protein_symbol = comp.Protein.symbol
            if math.isnan(comp.CompetitionRatio.group_cr):
                data.append({"x": protein_symbol, "y": 0, "position": comp.Residue.position })
            else:
                if protein_symbol in max_cr_dict and comp.CompetitionRatio.group_cr <= max_cr_dict[protein_symbol]:
                    continue
                max_cr_dict[protein_symbol] = comp.CompetitionRatio.group_cr
                data.append({"x": protein_symbol, "y": round(min(comp.CompetitionRatio.group_cr, 20), 2), "position": comp.Residue.position})


        cr_filtered = comp_proteins.filter(CompetitionRatio.group_cr > 4).order_by(CompetitionRatio.group_cr.desc()).limit(1000)




        protein_tbl = {
            'Protein': [],
            'Symbol': [],
            'Position': [],
            'Cell Type': [],
            'Concentration': [],
            'Concentration Unit': [],
            'Time': [],
            'Time Unit': [],
            'CR': []
        }

        for p in cr_filtered:
            protein_tbl['Protein'].append(p.Protein.uniprot)
            protein_tbl['Symbol'].append(p.Protein.symbol)
            protein_tbl['Position'].append(p.Residue.position)
            protein_tbl['Cell Type'].append(p.CellType.name)
            protein_tbl['Concentration'].append(p.CompoundTreatment.concentration)
            protein_tbl['Concentration Unit'].append(p.CompoundTreatment.concentrationunits)
            protein_tbl['Time'].append(p.CompoundTreatment.time)
            protein_tbl['Time Unit'].append(p.CompoundTreatment.timeunits)
            protein_tbl['CR'].append(round(p.CompetitionRatio.group_cr,2))


        protein_tbl['Concentration'] = [f"{c} {cu}" for c, cu in zip(protein_tbl['Concentration'], protein_tbl['Concentration Unit'])]
        protein_tbl['Time'] = [f"{t} {tu}" for t, tu in zip(protein_tbl['Time'], protein_tbl['Time Unit'])]



        df = pd.DataFrame(protein_tbl)

        # Drop duplicate rows based on all columns
        df_no_duplicates = df.drop_duplicates()

        # Update the original dictionary with the cleaned data
        protein_tbl = {
            'Protein': df_no_duplicates['Protein'].tolist(),
            'Symbol': df_no_duplicates['Symbol'].tolist(),
            'Position': df_no_duplicates['Position'].tolist(),
            'Cell Type': df_no_duplicates['Cell Type'].tolist(),
            'Concentration': df_no_duplicates['Concentration'].tolist(),
            'Time': df_no_duplicates['Time'].tolist(),
            'CR': df_no_duplicates['CR'].tolist(),
        }



    return render_template('compound_search.html', compound= compound, comp_image=comp_image,comp_smiles=comp_smiles, 
    data = data, celltypes=celltypes, selected_celltype=celltype, protein_tbl=protein_tbl, selected_residue_type= residue_type)



@app.route('/all_compounds')
def all_compounds():
    # Get all compounds in one query
    all_compounds = db.session.query(Compound.id, Compound.image, Compound.smiles) \
        .all()

    data = {
        'Compound': [],
        'Image': [],
        'CYS hits': [],
        'LYS hits': [],
        'Smiles': []
    }

    # Prepare a list of compound IDs for later use
    compound_ids = [compound.id for compound in all_compounds]

    # Create a single query to fetch the hit count for each compound
    hit_query = db.session.query(
            Compound.id, Residue.type, db.func.sum(db.case(
                    (CompetitionRatio.group_cr >= 4, 1),
                    else_=0
                )
            ).label('hit_count')
        ) \
        .select_from(CompetitionRatio) \
        .join(Residue, CompetitionRatio.residue_id == Residue.id) \
        .join(CompoundTreatment, CompetitionRatio.compoundtreatment_id == CompoundTreatment.id) \
        .join(Compound, CompoundTreatment.compound_id == Compound.id) \
        .filter(
            Compound.id.in_(compound_ids),
            CompetitionRatio.display_flag == True,
            CompetitionRatio.multimapper == False,
            CompetitionRatio.control_rsd <= 30
        ) \
        .group_by(Compound.id, Residue.type) \
        .distinct(Compound.id, Residue.type).subquery()

    # Fetch the results from the subquery
    hit_results = db.session.query(
        hit_query.c.id, hit_query.c.type, hit_query.c.hit_count
    ).all()

    # Populate the data dictionary with compound information
    for compound in all_compounds:
        compound_id = compound.id
        cys_hits = next((row.hit_count for row in hit_results if row.id == compound_id and row.type == 'CYS'), '-')
        lys_hits = next((row.hit_count for row in hit_results if row.id == compound_id and row.type == 'LYS'), '-')

        data['Compound'].append(compound_id)
        data['Image'].append(compound.image)
        data['CYS hits'].append(cys_hits)
        data['LYS hits'].append(lys_hits)
        data['Smiles'].append(compound.smiles)

    return render_template('all_compounds.html', data=data)





@app.route('/global_proteomics')
def global_proteomics():



    pathway_data = {}

    pathway_query = db.session.query(PathwayList) \
        .join(ProteinToPathway, PathwayList.id == ProteinToPathway.pathway_id) \
        .join(GpProtein, ProteinToPathway.protein_id == GpProtein.uniprot) \
        .with_entities(PathwayList.name, GpProtein.uniprot, GpProtein.symbol, GpProtein.organism) \
        .order_by(PathwayList.name)

    for path in pathway_query:
        pathway_name, uniprot, symbol, organism = path

        if pathway_name not in pathway_data:
            pathway_data[pathway_name] = []

        pathway_data[pathway_name].append({'protein_id': uniprot, 'symbol': symbol, 'organism': organism })



    # pathway_list = [path[1] for path in pathway_list]

    # print(pathway_list)


    cpd_table = {
        'Compound': [],
        'concentration': [],
        'concentration 2': [],
        'concentration Units': [],
        'time' : [],
        'organism': [],
        'cell_type': [],
        'project': [],
        'temperature': []
        
    }

    compound_list = db.session.query(GpCompoundTreatment).filter(GpCompoundTreatment.compound_id.notin_(['DMSO', 'DMSO-2h', 'DMSO-24h'])).\
        join(GpCellType, GpCompoundTreatment.celltype_id == GpCellType.id).\
        distinct(GpCompoundTreatment.compound_id)
    
    compound_list = sorted(compound_list, key=lambda x: x.compound_id, reverse=True)


    # compound_table = pd.DataFrame(compound_list)


    for cpd in compound_list:
        cpd_table['Compound'].append(cpd.compound_id)
        cpd_table['concentration'].append(f'{cpd.concentration}{cpd.concentrationunits}')
        cpd_table['concentration 2'].append(f'{cpd.concentration_two}')
        cpd_table['time'].append(f'{cpd.time}{cpd.timeunits}')
        cpd_table['organism'].append(cpd.organism)
        cpd_table['cell_type'].append(cpd.gpcelltype.name)
        cpd_table['project'].append(cpd.project)        
        cpd_table['temperature'].append(cpd.temperature)


    selected_cpd = request.args.get('cpd')
    data = None
    compound = None
    if selected_cpd:
        # Query database to get updated data for the selected compound
        vol_plot = db.session.query(FoldChange).\
                    join(GpCompoundTreatment, FoldChange.compoundtreatment_id == GpCompoundTreatment.id).\
                    join(GpProtein, FoldChange.protein_id == GpProtein.uniprot).\
                    filter(GpCompoundTreatment.compound_id == selected_cpd).\
                    distinct(FoldChange.protein_id, FoldChange.scan, FoldChange.foldchange, FoldChange.p_value).\
                    all()

        # Process the data
        updated_data = []
        for f in vol_plot:
            p_value = f.p_value
            if p_value != 0:
                log_p_value = -math.log10(p_value)
            else:
                log_p_value = float('inf')
            updated_data.append({
                'foldchange': f.foldchange,
                'p_value': log_p_value,
                'protein_id': f.protein_id,
                'symbol': f.gpprotein.symbol
            })



        comp_treat = db.session.query(GpCompoundTreatment).\
                join(GpCellType, GpCompoundTreatment.celltype_id == GpCellType.id).\
                filter(GpCompoundTreatment.compound_id == selected_cpd).\
                first()

        compound_treatment = {
            'Compound': '',
            'Concentration': '',
            'Concentration 2': '',
            'Concentration Units': '',
            'Plex': '',
            'Time': '',
            'Temperature': '',
            'Project': '',
            'Organism': '',
            'Celltype': '',
            'Data_acq': '',
            'Comment': ''


        }

        compound_treatment['Compound'] = comp_treat.compound_id
        compound_treatment['Concentration']= comp_treat.concentration
        compound_treatment['Concentration 2'] = comp_treat.concentration_two
        compound_treatment['Concentration Units']=comp_treat.concentrationunits
        compound_treatment['Plex']=comp_treat.plex_id
        compound_treatment['Time']=str(comp_treat.time) + comp_treat.timeunits
        compound_treatment['Temperature']=comp_treat.temperature
        compound_treatment['Organism']=comp_treat.organism
        compound_treatment['Celltype']=comp_treat.gpcelltype.name
        compound_treatment['Project']=comp_treat.project
        compound_treatment['Data_acq']=comp_treat.data_acq
        compound_treatment['Comment']=comp_treat.comments




        # Return the updated data as JSON
        return json.dumps({'updated_data': updated_data, 'compound_treatment': compound_treatment}) 
    else:
        compound = db.session.query(GpCompound).filter(GpCompound.id.notin_(['DMSO', 'DMSO-2h', 'DMSO-24h'])).order_by(asc(GpCompound.id)).first()
        compound = compound.id
        vol_plot = db.session.query(FoldChange).\
                join(GpCompoundTreatment, FoldChange.compoundtreatment_id == GpCompoundTreatment.id).\
                join(GpProtein, FoldChange.protein_id == GpProtein.uniprot).\
                join(GpCellType, GpCompoundTreatment.celltype_id == GpCellType.id).\
                filter(GpCompoundTreatment.compound_id == compound).\
                distinct(FoldChange.protein_id, FoldChange.scan, FoldChange.foldchange, FoldChange.p_value).\
                all()


        comp_treat = db.session.query(GpCompoundTreatment).\
                join(GpCellType, GpCompoundTreatment.celltype_id == GpCellType.id).\
                filter(GpCompoundTreatment.compound_id == compound).\
                first()

        compound_treatment = {
            'Compound': '',
            'Concentration': '',
            'Concentration 2': '',
            'Concentration Units': '',
            'Plex': '',
            'Time': '',
            'Temperature': '',
            'Organism': '',
            'Project':'',
            'Celltype': '',
            'Data_acq': '',
            'Comment': ''


        }

        compound_treatment['Compound'] = comp_treat.compound_id
        compound_treatment['Concentration']= comp_treat.concentration
        compound_treatment['Concentration 2'] = comp_treat.concentration_two
        compound_treatment['Concentration Units']=comp_treat.concentrationunits
        compound_treatment['Plex']=comp_treat.plex_id
        compound_treatment['Time']=str(comp_treat.time) + comp_treat.timeunits
        compound_treatment['Temperature']=comp_treat.temperature
        compound_treatment['Organism']=comp_treat.organism
        compound_treatment['Celltype']=comp_treat.gpcelltype.name
        compound_treatment['Project']=comp_treat.project
        compound_treatment['Data_acq']=comp_treat.data_acq
        compound_treatment['Comment']=comp_treat.comments



        data = []
        
        for f in vol_plot:
            p_value = f.p_value
            if p_value != 0:  # Avoid taking log(0)
                log_p_value = -math.log10(p_value)
            else:
                log_p_value = float('inf')  # Handling log(0) which goes to infinity
            data.append({
                'foldchange': f.foldchange,
                'p_value': log_p_value,
                'protein_id': f.protein_id,
                'symbol': f.gpprotein.symbol
            })



        
    return render_template('global_proteomics.html', cpd_table=cpd_table, compound=compound, data=data, 
    compound_treatment=compound_treatment, pathway_data=pathway_data) 

def nudge(pos, x_shift, y_shift):
    return {n:(x + x_shift, y + y_shift) for n,(x,y) in pos.items()}

def network_analysis(df_significant, compound):
    # Ensure folder exists
    plot_dir = os.path.join('static', 'images', 'network')
    os.makedirs(plot_dir, exist_ok=True)

    plot_filename = f'{compound}_network.png'
    plot_path = os.path.join(plot_dir, plot_filename)

    nodes, edges = enrichment_map(df_significant, column ='NOM p-val', cutoff= 0.05, top_term = 25)

    G = nx.from_pandas_edgelist(edges,
                                source='src_idx',
                                target='targ_idx',
                                edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

    fig, ax = plt.subplots(figsize=(30, 20))

    # Init node coordinates
    pos = nx.spring_layout(G, weight='jaccard_coef', seed=6)


    # Prepare node sizes and colors
    node_sizes = [nodes.loc[n, 'Hits_ratio'] * 1000 for n in G.nodes]
    node_colors = [nodes.loc[n, 'NES'] for n in G.nodes]

    # Draw nodes with colormap
    cmap = plt.cm.RdYlBu_r
    nodes_drawn = nx.draw_networkx_nodes(
        G,
        pos=pos,
        node_color=node_colors,
        node_size=node_sizes,
        cmap=cmap,
        ax=ax
    )

    pos_nodes = nudge(pos, 0, 0.04) 

    # Draw node labels
    labels = {n: nodes.loc[n, 'Term'] for n in G.nodes}
    nx.draw_networkx_labels(G, pos=pos_nodes, labels=labels, ax=ax, font_size=10, 
        horizontalalignment= 'center', font_weight = 'heavy')

    # Draw edges
    edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
    nx.draw_networkx_edges(
        G,
        pos=pos,
        width=[w * 10 for w in edge_weight],
        edge_color='#CDDBD4',
        ax=ax
    )

    # Add colorbar as legend for NES values
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min(node_colors), vmax=max(node_colors)))
    sm.set_array([])  
    cbar = plt.colorbar(sm, ax=ax, shrink=0.8, pad=0.005)  
    cbar.set_label('Normalised Enrichment Score (NSE)', fontsize=14)

    plt.axis('off')  # Hide axes    

    plot_filename = f'{compound}_network.png'

    plot_path = os.path.join('static', 'images/network/', plot_filename)

    fig.savefig(plot_path, bbox_inches='tight', pad_inches=0.1)
    plt.close(fig)
    fig.tight_layout()

    network_url = url_for('static', filename=f'images/network/{plot_filename}')
    return network_url

@app.route('/global_proteomics/image', methods=['POST'])
def kegg_analysis():

    # Get JSON data from the request
    path_data = request.json.get("path_analysis")
    compound = request.json.get("compound")

    organism = request.json.get("org_in")

    BASE_DIR = os.path.dirname(os.path.abspath(__file__))  

    human_gene_set = os.path.join(BASE_DIR, "human_kegg_pathway.gmt")
    mouse_gene_set = os.path.join(BASE_DIR, "mouse_kegg_pathway.gmt")

    path_data = pd.DataFrame(path_data)

        
    path_data = path_data.rename(columns = {'symbol': 'Gene','foldchange':'log2FoldChange', 'p_value':'padj'})
    path_data['padj'] = 10 ** -path_data['padj']
    
    path_data['Rank'] = -np.log10(path_data.padj)*path_data.log2FoldChange

    path_data = path_data.sort_values('Rank', ascending = False).reset_index(drop = True)

    path_data.dropna()

    ranking = path_data[['Gene', 'Rank']]
    

    if organism == 'human':
        pre_res = gp.prerank(rnk=ranking, gene_sets=human_gene_set, seed=6, permutation_num=1000)
        
        df_data = pd.DataFrame(pre_res.res2d)

        # Filter significant pathways
        df_significant = df_data[(df_data['NOM p-val'] <= 0.05) & (df_data['FDR q-val'] <= 0.25)]

        if df_significant.empty:
            print("No significant pathways found.")
        else:
            # Sort pathways by NES
            df_plot = df_significant.reindex(df_significant['NES'].abs().sort_values(ascending=False).index)
            df_plot = df_plot.head(50)
            df_plot = df_plot.sort_values(by='NES', ascending=False)

            # Normalize p-values for coloring
            cmap = plt.cm.RdYlBu_r  # Reverse color map for better contrast
            norm = mcolors.Normalize(vmin=df_plot['NOM p-val'].min(), vmax=df_plot['NOM p-val'].max())

            # Map colors efficiently
            df_plot["color"] = df_plot["NOM p-val"].map(lambda p: cmap(norm(p)))

            # Create figure and axis
            fig, ax = plt.subplots(figsize=(12, 10))

            net_url = network_analysis(df_significant, compound)
            
    else:

        pre_res = gp.prerank(rnk = ranking, gene_sets = mouse_gene_set, seed = 6, permutation_num = 1000)
    
        df_data= pd.DataFrame(pre_res.res2d)

        # Filter significant pathways
        df_significant = df_data[(df_data['NOM p-val'] <= 0.05) & (df_data['FDR q-val'] <= 0.25)]
        if df_significant.empty:
            print("No significant pathways found.")
            
        else:
            # Sort pathways by NES
            df_plot = df_significant.reindex(df_significant['NES'].abs().sort_values(ascending=False).index)
            df_plot = df_plot.head(50)
            df_plot = df_plot.sort_values(by='NES', ascending=False)
            
            # Normalize p-values for coloring
            cmap = plt.cm.RdYlBu_r  # Reverse color map for better contrast
            norm = mcolors.Normalize(vmin=df_plot['NOM p-val'].min(), vmax=df_plot['NOM p-val'].max())

            # Map colors efficiently
            df_plot["color"] = df_plot["NOM p-val"].map(lambda p: cmap(norm(p)))

            # Create figure and axis
            fig, ax = plt.subplots(figsize=(12, 10))

            net_url = network_analysis(df_significant, compound)


    # Plot barplot
    bars = ax.barh(df_plot["Term"], df_plot["NES"], color=df_plot["color"])

    # Labels and styling
    ax.axvline(0, color='black', linewidth=1)
    ax.set_xlabel("Normalised Enrichment Score (NES)")
    ax.set_ylabel("Pathway")
    ax.set_title("Two-Sided Barplot of NES (Significant Pathways Only)")

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Normalised p-value")

    plot_dir = os.path.join('static', 'images', 'gsea')
    os.makedirs(plot_dir, exist_ok=True)

    plot_filename = f'{compound}.png'
    plot_path = os.path.join(plot_dir, plot_filename)

    fig.savefig(plot_path, bbox_inches='tight')
    plt.close(fig)

    plot_url = url_for('static', filename=f'images/gsea/{plot_filename}')

    df_significant = df_significant

    df_significant = df_significant[['Term', 'ES', 'NES', 'NOM p-val', 
    'FDR q-val', 'Tag %', 'Gene %', 'Lead_genes']]

    df_significant['Lead_genes']=  df_significant['Lead_genes'].str.replace(';', ', ')

    df_dict = df_significant.to_dict(orient='records')

    

    return jsonify({'plot_url': plot_url, 'gsea_data': df_dict, 'net_url': net_url})
 


@app.route('/gp_compound_images/<filename>')
def gp_compound_images(filename):
    image_dir = '/resources/smiles-structure/'
    return send_from_directory(image_dir, filename)




@app.route('/gp_autocomplete')
def gp_autocomplete():
    gp_protein = request.args.get('gp_protein')  # changed from 'target_id' to 'term'
    matching_proteins = db.session.query(GpProtein).filter(
        GpProtein.symbol.ilike(f'%{gp_protein}%')
    ).limit(20).all()
    matching_uniprots = db.session.query(GpProtein).filter(
        GpProtein.uniprot.ilike(f'%{gp_protein}%')
    ).limit(20).all()
    matching_synonyms = db.session.query(GpProteinSynonym).filter(
        GpProteinSynonym.synonym.ilike(f'%{gp_protein}%')
    ).limit(20).all()

    gp_matching_items = [protein.symbol for protein in matching_proteins] \
        + [synonym.synonym for synonym in matching_synonyms] \
        + [uniprot.uniprot for uniprot in matching_uniprots]
    gp_matching_items = list(set(item.upper() for item in gp_matching_items))
    print(gp_matching_items)
    return jsonify(gp_matching_items)



@app.route('/protein_compound_gp')
def protein_compound_gp():
    target_protein = request.args.get('gp_auto').upper() if request.args.get('gp_auto') is not None else None
    data = []



    protein = db.session.query(GpProtein).filter(func.upper(GpProtein.symbol) == target_protein).first()

    uniprot = db.session.query(GpProtein).filter(func.upper(GpProtein.uniprot) == target_protein).first()

    synonyms = (
        db.session.query(GpProteinSynonym.synonym, func.upper(GpProtein.symbol))
        .join(GpProtein, GpProteinSynonym.uniprot == GpProtein.uniprot)
        .filter(func.upper(GpProteinSynonym.synonym) == target_protein)
        .first()
    )


    if protein is not None:
        target_protein = protein.symbol
    
    elif uniprot is not None:
        target_protein = uniprot.symbol

    elif synonyms is not None:
        target_protein = synonyms.synonym
    plot_url = ''
    if target_protein:

        vol_plot = db.session.query(
            GpProtein.uniprot,
            GpProtein.symbol,
            GpProtein.organism,
            GpCompoundTreatment.compound_id,
            GpCompoundTreatment.concentration,
            GpCompoundTreatment.concentrationunits,
            GpCompoundTreatment.time,
            GpCompoundTreatment.timeunits,
            GpCompoundTreatment.temperature, 
            GpCompoundTreatment.plex_id,
            FoldChange.foldchange,
            FoldChange.p_value,
            GpCellType.name,
            GpCompound.smiles,
            GpCompoundTreatment.project
        ).join(FoldChange, GpProtein.uniprot == FoldChange.protein_id
        ).join(GpCompoundTreatment, FoldChange.compoundtreatment_id == GpCompoundTreatment.id
        ).join(GpCellType, GpCompoundTreatment.celltype_id == GpCellType.id
        ).join(GpCompound, GpCompoundTreatment.compound_id == GpCompound.id
        ).filter(GpProtein.symbol == target_protein
        ).distinct(FoldChange.protein_id, FoldChange.foldchange, FoldChange.p_value
        ).all()

        for f in vol_plot:
            p_value = f.p_value
            log_p_value = -math.log10(p_value) if p_value > 0 else float('inf')
            data.append({
                'foldchange': round(f.foldchange, 2),
                'p_value': round(log_p_value, 2),
                'uniprot': f.uniprot,
                'symbol': f.symbol,
                'compound': f.compound_id,
                'cpd_image': url_for('gp_compound_images', filename=f'{f.compound_id}.png'),
                'smiles': f.smiles,
                'concentration': f'{f.concentration} {f.concentrationunits}',
                'time': f'{f.time} {f.timeunits}',
                'temperature': f.temperature,
                'celltype': f.name,
                'project': f.project,
                'organism': f.organism,
                'experiment': f.plex_id
            })




        return json.dumps({'updated_data': data})

    return render_template('protein_compound_gp.html', data=data, protein_id=target_protein)











#DiscoveromeGPT

@app.route('/discoverome_gpt')
def discoverome_gpt():
    session.setdefault('chat_history', [])
    return render_template('discoverome_gpt.html')


def _find_proteins_in_message(message):
    """Match tokens in message against protein symbols, UniProt IDs, and synonyms."""
    tokens = list({t.upper() for t in re.findall(r'[A-Za-z][A-Za-z0-9\-]{1,}', message)})
    if not tokens:
        return []
    proteins = (
        db.session.query(Protein)
        .filter(or_(Protein.symbol.in_(tokens), Protein.uniprot.in_(tokens)))
        .limit(5)
        .all()
    )
    if not proteins:
        syns = (
            db.session.query(ProteinSynonym)
            .filter(ProteinSynonym.synonym.in_(tokens))
            .limit(5)
            .all()
        )
        proteins = [s.protein for s in syns if s.protein]
    return list({p.uniprot: p for p in proteins}.values())


def _get_compound_hits(uniprot, limit=15):
    """Best compound hits for a protein with group_cr >= 4, one row per compound."""
    rows = (
        db.session.query(
            CompoundTreatment.compound_id,
            Residue.position,
            Residue.type,
            CellType.name.label('celltype'),
            func.max(CompetitionRatio.group_cr).label('max_cr'),
        )
        .join(Residue, CompetitionRatio.residue_id == Residue.id)
        .join(Protein, Residue.uniprot == Protein.uniprot)
        .join(CompoundTreatment, CompetitionRatio.compoundtreatment_id == CompoundTreatment.id)
        .join(CellType, CompoundTreatment.celltype_id == CellType.id)
        .filter(
            Protein.uniprot == uniprot,
            CompetitionRatio.group_cr >= 4,
            CompetitionRatio.display_flag == True,
            CompetitionRatio.multimapper == False,
            CompetitionRatio.control_rsd <= 30,
            CompoundTreatment.compound_id != 'DMSO',
        )
        .group_by(CompoundTreatment.compound_id, Residue.position, Residue.type, CellType.name)
        .order_by(func.max(CompetitionRatio.group_cr).desc())
        .limit(limit)
        .all()
    )
    return [
        {
            'compound': r.compound_id,
            'position': r.position,
            'residue_type': r.type,
            'celltype': r.celltype,
            'max_cr': round(r.max_cr, 2),
        }
        for r in rows
    ]


def _fetch_uniprot(uniprot_id):
    """Fetch protein function, disease associations, and subcellular location from UniProt."""
    try:
        r = http_requests.get(
            f'https://rest.uniprot.org/uniprotkb/{uniprot_id}',
            params={'format': 'json'},
            timeout=6,
        )
        if r.status_code != 200:
            return {}
        d = r.json()

        name = (
            d.get('proteinDescription', {})
            .get('recommendedName', {})
            .get('fullName', {})
            .get('value', '')
        )
        functions = [
            c['texts'][0]['value']
            for c in d.get('comments', [])
            if c.get('commentType') == 'FUNCTION' and c.get('texts')
        ]
        diseases = [
            c.get('disease', {}).get('diseaseName', {}).get('value', '')
            for c in d.get('comments', [])
            if c.get('commentType') == 'DISEASE'
        ]
        subcell = list({
            loc.get('location', {}).get('value', '')
            for c in d.get('comments', [])
            if c.get('commentType') == 'SUBCELLULAR LOCATION'
            for loc in c.get('subcellularLocations', [])
        })
        # Ensembl gene ID for Open Targets
        ensembl_gene = next(
            (
                prop['value']
                for ref in d.get('uniProtKBCrossReferences', [])
                if ref.get('database') == 'Ensembl'
                for prop in ref.get('properties', [])
                if prop.get('key') == 'GeneId'
            ),
            None,
        )
        return {
            'name': name,
            'function': ' '.join(functions)[:700],
            'diseases': [dd for dd in diseases if dd][:6],
            'subcellular_location': [s for s in subcell if s][:4],
            'ensembl_gene': ensembl_gene,
        }
    except Exception:
        return {}


def _fetch_opentargets(ensembl_gene_id):
    """Fetch disease associations and tractability from Open Targets using Ensembl gene ID."""
    if not ensembl_gene_id:
        return {}
    query = """
    query targetInfo($id: String!) {
      target(ensemblId: $id) {
        id
        approvedName
        approvedSymbol
        functionDescriptions
        associatedDiseases(page: {index: 0, size: 5}) {
          rows {
            disease { id name }
            score
          }
        }
        tractability {
          value
          modality
          label
        }
      }
    }
    """
    try:
        r = http_requests.post(
            'https://api.platform.opentargets.org/api/v4/graphql',
            json={'query': query, 'variables': {'id': ensembl_gene_id}},
            timeout=8,
        )
        if r.status_code != 200:
            return {}
        t = r.json().get('data', {}).get('target') or {}
        if not t:
            return {}
        diseases = [
            f"{row['disease']['name']} (score: {round(row['score'], 2)})"
            for row in t.get('associatedDiseases', {}).get('rows', [])
        ]
        tractability = [
            f"{tr['modality']}: {tr['label']}"
            for tr in (t.get('tractability') or [])
            if tr.get('value')
        ]
        return {
            'ot_name': t.get('approvedName', ''),
            'functions': (t.get('functionDescriptions') or [])[:2],
            'diseases': diseases[:5],
            'tractability': tractability[:5],
        }
    except Exception:
        return {}


@app.route('/api/chat', methods=['POST'])
def api_chat():
    data = request.get_json(silent=True) or {}
    user_msg = (data.get('message') or '').strip()
    if not user_msg:
        return jsonify({'error': 'empty message'}), 400

    history = session.get('chat_history', [])

    proteins = _find_proteins_in_message(user_msg)

    context_parts = []
    sql_shown = None

    for protein in proteins[:3]:
        hits = _get_compound_hits(protein.uniprot)
        uniprot_data = _fetch_uniprot(protein.uniprot)
        ot_data = _fetch_opentargets(uniprot_data.get('ensembl_gene'))

        block = [f"### Protein: {protein.symbol} ({protein.uniprot})"]
        if protein.description:
            block.append(f"Description: {protein.description}")
        if uniprot_data.get('name'):
            block.append(f"Full name: {uniprot_data['name']}")
        if uniprot_data.get('function'):
            block.append(f"Function (UniProt): {uniprot_data['function']}")
        if uniprot_data.get('subcellular_location'):
            block.append(f"Subcellular location: {', '.join(uniprot_data['subcellular_location'])}")
        if uniprot_data.get('diseases'):
            block.append(f"UniProt disease links: {', '.join(uniprot_data['diseases'])}")
        if ot_data.get('diseases'):
            block.append(f"Open Targets disease associations: {'; '.join(ot_data['diseases'])}")
        if ot_data.get('tractability'):
            block.append(f"Open Targets tractability: {'; '.join(ot_data['tractability'])}")

        if hits:
            if sql_shown is None:
                sql_shown = (
                    f"SELECT ct.compound_id, r.position, r.type, c.name AS celltype,\n"
                    f"       MAX(cr.group_cr) AS max_cr\n"
                    f"FROM chemoproteomics.competitionratios cr\n"
                    f"JOIN chemoproteomics.residues r ON cr.residue_id = r.id\n"
                    f"JOIN core.proteins p ON r.uniprot = p.uniprot\n"
                    f"JOIN core.compoundtreatments ct ON cr.compoundtreatment_id = ct.id\n"
                    f"JOIN core.celltypes c ON ct.celltype_id = c.id\n"
                    f"WHERE p.uniprot = '{protein.uniprot}'\n"
                    f"  AND cr.group_cr >= 4\n"
                    f"  AND cr.display_flag = TRUE\n"
                    f"  AND cr.multimapper = FALSE\n"
                    f"  AND cr.control_rsd <= 30\n"
                    f"  AND ct.compound_id != 'DMSO'\n"
                    f"GROUP BY ct.compound_id, r.position, r.type, c.name\n"
                    f"ORDER BY max_cr DESC\n"
                    f"LIMIT 15;"
                )
            block.append(f"\nTop compound hits from Discoverome (group_cr >= 4, {len(hits)} results):")
            for h in hits[:12]:
                block.append(
                    f"  • {h['compound']}  CR={h['max_cr']}  "
                    f"residue={h['residue_type']}{h['position']}  cell={h['celltype']}"
                )
        else:
            block.append("No compound hits with group_cr >= 4 found in Discoverome database.")

        context_parts.append('\n'.join(block))

    system_prompt = (
        "You are DiscoveromeGPT, an expert AI assistant embedded in Discoverome — a chemoproteomics and "
        "global proteomics platform for drug discovery. "
        "Answer questions about proteins, compounds, binding data, disease associations, and druggability. "
        "When database context is provided, summarise the key findings clearly: highlight the most potent "
        "compounds (highest CR), relevant disease context, and Open Targets tractability insights. "
        "CR (Competition Ratio) >= 4 means significant covalent compound binding at a residue. "
        "Use concise markdown. If no database context is present, answer from general biochemistry knowledge."
    )

    llm_messages = [{'role': 'system', 'content': system_prompt}]
    llm_messages += history[-10:]

    user_content = user_msg
    if context_parts:
        user_content = f"{user_msg}\n\n---\n**Database & API context:**\n\n{''.join(context_parts)}"

    llm_messages.append({'role': 'user', 'content': user_content})

    ollama_host = os.getenv('OLLAMA_HOST', 'http://localhost:11434')
    try:
        response = ollama.Client(host=ollama_host).chat(
            model='llama3',
            messages=llm_messages,
            options={'temperature': 0.3},
        )
        if hasattr(response, 'message'):
            answer = response.message.content
        else:
            answer = response['message']['content']
    except Exception as exc:
        answer = f"LLM unavailable: {exc}. Check that Ollama is running with `ollama serve`."

    history.append({'role': 'user', 'content': user_msg})
    history.append({'role': 'assistant', 'content': answer})
    session['chat_history'] = history[-30:]

    return jsonify({'response': answer, 'sql': sql_shown})


@app.route('/api/chat/reset', methods=['POST'])
def api_chat_reset():
    session['chat_history'] = []
    return jsonify({'ok': True})




if (__name__ == '__main__'):
    app.run()




