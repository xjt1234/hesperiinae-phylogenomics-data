#!/usr/bin/env python3
"""Rebuild an auditable three-state table under the historical 102 tip labels.

This is a provenance reconstruction, not a new specimen identification or screen.
Only exact joins and the frozen, explicit historical source chains below are used.
Old imputed zeros become unknown. Input files are never modified. Existing output
files are accepted only when byte-identical; otherwise this script refuses them.
Python 3 standard library only.
"""
import argparse
import collections
import csv
import hashlib
import io
import json
from pathlib import Path
import re

DEFAULT_BASE = Path('/home/data/t200301/xjt/Wolchbia')
UNKNOWN = {
    'Potanthus_ganda', 'Polytremis_nascens', 'Thymelicus_sylvestris',
    'Oarisma_garita', 'Polites_origenes_origenes', 'Wallengrenia_otho',
}
# Exact source-record IDs reviewed against infection_status.tsv and the original
# summary filenames. These are historical provenance links, NOT a claim that
# the dated-tree specimen and screened specimen are the same individual.
SOURCE_OVERRIDES = {
    'Aeromachus_piceus': 'Pedesta_xiaoqingae',
    'Ampittia_trimacula': 'Ampittia_trimacula_SRR28513893',
    'Pedesta_masuriensis': 'Pedesta_masuriensis_raw',
    'Astictopterus_jama': 'Astictopterus-jama_jama_SRR7174401',
    'Cephrenes_acalle_oceanica': 'Cephrenes_acalle_oceanica_raw',
    'Erionota_thrax': 'Erionota_thrax_SRR7174519',
    'Udaspes_foluss': 'Udaspes_folus',
    'Iton_watsonii': 'Iton_watsoni',
    'Gegenes_nostrodamus': 'Gegenes_nostrodamus_ERR13800478',
    'Pelopidas_agna': 'Pelopidas_agna_SRR28430792',
    'Agathymus_mariae_mariae': 'Agathymus_mariae_mariae_SRR7174357',
    'Megathymus_yuccae_yuccae': 'Megathymus_yuccae_yuccae_SRR25297859',
    'Megathymus_violae': 'Megathymus_violae_SRR25297858',
    'Megathymus_ursus_ursus': 'Megathymus_ursus_deserti_SRR25297841',
    'Carystus_phorcus': 'Carystus_phorcus_SRR7174425',
    'Molo_mango': 'Molo_mango_SRR7174583',
    'Calpodes_ethlius': 'Calpodes_ethlius_SRR7174423',
    'Thymelicus_lineola': 'Thymelicus_lineola_ERR10851543',
    'Copaeodes_aurantiaca': 'Copaeodes_aurantiaca_SRR7174406',
    'Oarisma_powesheik': 'Oarisma_powesheik_SRR19395065',
    'Lerema_liris': 'Lerema_liris_SRR23000799',
    'Asbolis_capucinus': 'Asbolis_capucinus_SRR7174410',
    'Ochlodes_sylvanus': 'Ochlodes_sylvanus_ERR6054647',
    'Ochlodes_thibetana': 'Ochlodes_subhyalina',
    'Hesperia_balcones': 'Hesperia_balcones_SRR23000750',
    'Hesperia_meskei_straton': 'Hesperia_meskei_straton_SRR13833996',
    'Hesperia_nevada_nevada': 'Hesperia_nevada_nevada_SRR13834000',
    'Hesperia_comma': 'Hesperia_comma_ERR6054630',
    'Hesperia_colorado_sublima': 'Hesperia_colorado_sublima_SRR13833921',
    'Hylephila_phyleus': 'Hylephila_phyleus_SRR15257228',
    'Polites_vibex_praeceps': 'Polites_vibex_praeceps_SRR15257222',
    'Polites_baracoa_baracoa': 'Polites_baracoa_baracoa_SRR15257214',
    'Polites_carus': 'Polites_carus_SRR15257215',
    'Polites_themistocles_turneri': 'Polites_themistocles_turneri_SRR15257208',
    'Polites_peckius_peckius': 'Polites_peckius_peckius_SRR15257209',
    'Polites_draco': 'Polites_draco_SRR15257210',
    'Polites_sabuleti_sabuleti': 'Polites_sabuleti_sabuleti_SRR15257211',
}
COMPOSITE_ROWS = {
    'Aeromachus_piceus': 'Aeromachus_piceus.Pedesta_xiaoqingae',
    'Astictopterus_jama': 'Astictopterus_jama.Astictopterus-jama_jama_SRR7174401',
    'Iton_watsonii': 'Iton_watsonii.Iton_watsoni',
    'Megathymus_ursus_ursus': 'Megathymus_ursus_ursus.Megathymus_ursus_deserti_SRR25297841',
    'Ochlodes_thibetana': 'Ochlodes_thibetana.Ochlodes_subhyalina',
    'Udaspes_foluss': 'Udaspes_foluss.Udaspes_folus',
}
OUTSIDE = {'Anthoptus_insignis_SRR7174402', 'Parnara_batta_SRR28430794'}
NAME_WARNINGS = {
    'Aeromachus_piceus': 'Historical cross-genus chain to Pedesta_xiaoqingae; taxon and specimen identity unresolved.',
    'Astictopterus_jama': 'Historical hyphen/subspecies spelling chain; specimen identity unresolved.',
    'Iton_watsonii': 'Historical watsonii/watsoni spelling chain; specimen identity unresolved.',
    'Megathymus_ursus_ursus': 'Historical ursus/deserti subspecies chain; not confirmed same taxon or specimen. Composite raw record has NA PASS count.',
    'Ochlodes_thibetana': 'Historical thibetana/subhyalina species chain; taxon and specimen identity unresolved.',
    'Udaspes_foluss': 'Historical foluss/folus spelling chain; specimen identity unresolved.',
}


def read_tsv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    for line, row in enumerate(rows, 2):
        row['_line'] = str(line)
    return rows


def index_unique(rows, key):
    result = {}
    for row in rows:
        if row[key] in result:
            raise ValueError('Duplicate exact key: ' + row[key])
        result[row[key]] = row
    return result


def sha256(path):
    if not path.is_file():
        return 'NA'
    return hashlib.sha256(path.read_bytes()).hexdigest()


def integer(value):
    return None if value in (None, '', 'NA') else int(value)


def summary_evidence(path):
    if not path.is_file():
        return {'sample': 'NA', 'pass_contigs': None, 'putative_contigs': None}
    text = path.read_text(encoding='utf-8')
    def one(pattern):
        values = re.findall(pattern, text, flags=re.MULTILINE)
        if len(values) != 1:
            raise ValueError('Expected exactly one summary field: ' + str(path) + ' / ' + pattern)
        return values[0]
    return {
        'sample': one(r'^Sample:\s*(.+)$').strip(),
        'pass_contigs': int(one(r'^\s*PASS contigs[^:]*:\s*(\d+)\s*$')),
        'putative_contigs': int(one(r'^\s*Putative contigs[^:]*:\s*(\d+)\s*$')),
    }


def encoded(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'))


def render_tsv(rows):
    stream = io.StringIO(newline='')
    fields = list(dict.fromkeys(key for row in rows for key in row))
    writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows({k: 'NA' if v is None else v for k, v in row.items()} for row in rows)
    return stream.getvalue().encode('utf-8')


def build(base, script_path):
    mapping_path = base / '26.3.24final/Hesperiinae_dated_pruned_102_mapping.tsv'
    old_path = base / 'host_infection_analysis_20260222/DATA/tip_infection_table.tsv'
    raw_path = base / 'infection_status.tsv'
    mapping = read_tsv(mapping_path)
    index_unique(mapping, 'final_taxon')
    old = index_unique(read_tsv(old_path), 'tree_tip')
    raw = index_unique(read_tsv(raw_path), 'host_tip')
    if len(mapping) != 102 or len(old) != 107 or len(raw) != 104:
        raise ValueError('Historical dataset sizes changed; re-audit before rebuilding.')
    manifest = []
    for role, path in [('exact_102_mapping', mapping_path), ('historical_tip_states', old_path),
                       ('historical_raw_states', raw_path), ('reconstruction_script', script_path)]:
        manifest.append(dict(record_type='input_file', input_role=role, file_path=str(path),
                             file_sha256=sha256(path)))
    states, assignments, evidence = [], {}, {}
    for mapping_row in mapping:
        taxon = mapping_row['final_taxon']
        historical = old[taxon]  # Exact label join only, including the old imputation flag.
        row = dict(final_taxon=taxon, dated_tree_tip_original=mapping_row['dated_tree_tip_original'],
                   mapping_file=str(mapping_path), mapping_line=mapping_row['_line'],
                   historical_tip_file=str(old_path), historical_tip_line=historical['_line'],
                   historical_infected=historical['infected'], historical_source=historical['source'],
                   historical_match_key=historical['match_key'],
                   screening_state='unknown', infected_tri='NA', threshold_pass_contigs=10,
                   raw_record_ids='[]', raw_record_lines='[]', raw_pass_counts='{}',
                   provenance_join='exact_final_taxon_to_historical_tree_tip',
                   source_file='NA', source_sha256='NA', summary_sample='NA',
                   summary_pass_contigs='NA', pass_ids_file='NA', pass_ids_sha256='NA',
                   pass_ids_count='NA', evidence_qc='no_screening_record_in_audited_sources',
                   conflicts='[]', identity_status='unverified_specimen_correspondence',
                   identity_note='Dated-tip naming and screening labels do not establish the same individual specimen.')
        if historical['source'] == 'imputed_missing_as_0':
            if taxon not in UNKNOWN or historical['infected'] != '0':
                raise ValueError('Unexpected historical imputation: ' + taxon)
            row['identity_note'] = 'No original screening record located in audited sources; never treat this old imputed zero as below threshold.'
            states.append(row)
            continue
        if historical['source'] != 'matched' or historical['infected'] not in ('0', '1'):
            raise ValueError('Unexpected historical state: ' + taxon)
        primary = SOURCE_OVERRIDES.get(taxon, taxon)
        ids = [primary] + ([COMPOSITE_ROWS[taxon]] if taxon in COMPOSITE_ROWS else [])
        source_record = raw[primary]  # Exact record ID, no normalization or token fallback.
        for record_id in ids:
            if record_id in assignments:
                raise ValueError('Raw record assigned to more than one tip: ' + record_id)
            assignments[record_id] = taxon
        path = Path(source_record['source_file'])
        item = summary_evidence(path)
        ids_path = Path(str(path).removesuffix('.wolbachia_summary.txt') + '.wol.pass.ids')
        ids_count = sum(bool(line.strip()) for line in ids_path.read_text().splitlines()) if ids_path.is_file() else None
        conflicts = []
        if not path.is_file():
            conflicts.append('missing_primary_summary')
        if item['sample'] != primary:
            conflicts.append('summary_sample_differs_from_exact_primary_raw_record')
        if item['pass_contigs'] is not None and int(item['pass_contigs'] >= 10) != int(historical['infected']):
            conflicts.append('summary_threshold_disagrees_with_historical_state')
        for record_id in ids:
            record = raw[record_id]
            if record['infected'] != historical['infected']:
                conflicts.append('raw_state_disagreement:' + record_id)
            n = integer(record['pass_contigs'])
            if n is not None and n != item['pass_contigs']:
                conflicts.append('raw_summary_PASS_disagreement:' + record_id)
        if ids_count is not None and ids_count != item['pass_contigs']:
            conflicts.append('summary_PASS_ids_disagreement')
        row.update(screening_state='threshold_positive' if historical['infected'] == '1' else 'below_threshold',
                   infected_tri=historical['infected'], raw_record_ids=encoded(ids),
                   raw_record_lines=encoded([raw[i]['_line'] for i in ids]),
                   raw_pass_counts=encoded({i: raw[i]['pass_contigs'] for i in ids}),
                   provenance_join='explicit_frozen_historical_source_chain' if taxon in SOURCE_OVERRIDES else 'exact_raw_host_tip',
                   source_file=str(path), source_sha256=sha256(path), summary_sample=item['sample'],
                   summary_pass_contigs=item['pass_contigs'], pass_ids_file=str(ids_path),
                   pass_ids_sha256=sha256(ids_path), pass_ids_count=ids_count,
                   evidence_qc='CONFLICT' if conflicts else ('summary_and_PASS_ids_agree' if ids_count is not None else 'summary_numeric_PASS_verified_ids_file_absent'),
                   conflicts=encoded(conflicts), identity_note=NAME_WARNINGS.get(taxon, row['identity_note']))
        evidence[taxon] = row
        states.append(row)
    if set(raw) - set(assignments) != OUTSIDE:
        raise ValueError('Unexpected raw records outside the frozen source chains: ' + repr(set(raw) - set(assignments)))
    for record_id, raw_row in raw.items():
        taxon = assignments.get(record_id)
        support = evidence.get(taxon, {})
        declared_path = Path(raw_row['source_file'])
        effective_path = Path(support.get('source_file', raw_row['source_file']))
        item = summary_evidence(effective_path)
        manifest.append(dict(record_type='historical_raw_record', input_role='infection_status_row',
                             file_path=str(raw_path), file_sha256=sha256(raw_path),
                             raw_record_line=raw_row['_line'], raw_host_tip=record_id,
                             raw_infected=raw_row['infected'], raw_pass_contigs=raw_row['pass_contigs'],
                             final_taxon=taxon or 'NA',
                             scope='outside_historical_102' if taxon is None else ('composite_duplicate_provenance_row' if record_id == COMPOSITE_ROWS.get(taxon) else 'selected_primary_source_row'),
                             declared_source_file=raw_row['source_file'], declared_source_sha256=sha256(declared_path),
                             effective_evidence_file=str(effective_path), effective_evidence_sha256=sha256(effective_path),
                             source_relation='explicit_compound_filename_chain_not_specimen_confirmation' if taxon and record_id == COMPOSITE_ROWS.get(taxon) else 'exact_source_file_field',
                             summary_sample=item['sample'], summary_pass_contigs=item['pass_contigs'],
                             pass_ids_file=support.get('pass_ids_file', 'NA'),
                             pass_ids_sha256=support.get('pass_ids_sha256', 'NA'),
                             pass_ids_count=support.get('pass_ids_count', 'NA'),
                             evidence_qc=support.get('evidence_qc', 'outside_target_set'),
                             conflicts=support.get('conflicts', '[]'),
                             identity_status='unverified_specimen_correspondence'))
    counts = collections.Counter(row['screening_state'] for row in states)
    if counts != {'threshold_positive': 43, 'below_threshold': 53, 'unknown': 6}:
        raise ValueError('Reconstructed historical counts changed: ' + repr(counts))
    if {r['final_taxon'] for r in states if r['infected_tri'] == 'NA'} != UNKNOWN:
        raise ValueError('Unknown states were lost or imputed.')
    flow = []
    def add(stage, category, count, note):
        flow.append(dict(stage=stage, category=category, n=count, interpretation=note))
    add('historical_raw_table', 'all_rows_not_unique_specimens', 104, '47 positive and 57 zero rows; contains duplicate provenance labels.')
    add('historical_raw_table', 'positive_rows', 47, 'Two duplicate compound labels and two outside-target positive labels explain 47 to 43.')
    add('historical_raw_table', 'zero_rows', 57, 'Four duplicate compound zero labels explain 57 to 53; specimen identities remain unresolved.')
    add('historical_raw_table', 'positive_compound_duplicate_rows', 2, 'Astictopterus compound label; Ochlodes compound label.')
    add('historical_raw_table', 'zero_compound_duplicate_rows', 4, 'Aeromachus/Pedesta; Iton; Megathymus; Udaspes compound labels.')
    add('historical_raw_table', 'outside_target_positive_rows', 2, 'Anthoptus_insignis_SRR7174402; Parnara_batta_SRR28430794.')
    add('historical_host_tree', 'all_tips', 107, '43 positive; 53 matched zero; 11 missing imputed zero by the old script.')
    add('historical_host_tree', 'old_imputed_zero_tips', 11, 'Five outgroups and six ingroup tips; none are confirmed below-threshold by imputation.')
    add('exact_dated_ingroup_mapping', 'all_tips', 102, 'Exact final_taxon join; 102 is target tree size, not evidence of 102 screened individuals.')
    for category in ('threshold_positive', 'below_threshold', 'unknown'):
        add('reconstructed_historical_labels', category, counts[category], 'Historical-label provenance only; specimen correspondence not confirmed.')
    add('reconstructed_historical_labels', 'numeric_evidence_available', 96, '43 plus 53 under historical source chains; conditional on specimen identity audit.')
    add('reconstructed_historical_labels', 'summary_PASS_ids_agree', sum(r['evidence_qc'] == 'summary_and_PASS_ids_agree' for r in states), 'Independent line-count comparison with original PASS ids when present.')
    add('reconstructed_historical_labels', 'evidence_conflicts', sum(r['conflicts'] != '[]' for r in states), 'No conflict is suppressed; any conflict causes nonzero exit after writing auditable results.')
    return {'screening_status_102.tsv': states, 'screening_source_manifest.tsv': manifest, 'screening_flow.tsv': flow}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', type=Path, default=DEFAULT_BASE)
    parser.add_argument('--outdir', type=Path, default=Path(__file__).resolve().parent.parent / '02_QC')
    parser.add_argument('--check-existing', action='store_true', help='Read-only deterministic reproduction check; write nothing.')
    args = parser.parse_args()
    tables = build(args.base.resolve(), Path(__file__).resolve())
    rendered = {name: render_tsv(rows) for name, rows in tables.items()}
    # Check every destination before creating any output; never overwrite a different file.
    for name, content in rendered.items():
        path = args.outdir / name
        if path.exists() and path.read_bytes() != content:
            raise SystemExit('Refusing to overwrite different existing output: ' + str(path))
        if args.check_existing and not path.is_file():
            raise SystemExit('Missing output for reproduction check: ' + str(path))
    if not args.check_existing:
        args.outdir.mkdir(parents=True, exist_ok=True)
        for name, content in rendered.items():
            path = args.outdir / name
            if not path.exists():
                with path.open('xb') as handle:
                    handle.write(content)
    counts = collections.Counter(r['screening_state'] for r in tables['screening_status_102.tsv'])
    conflicts = [r['final_taxon'] for r in tables['screening_status_102.tsv'] if r['conflicts'] != '[]']
    print(json.dumps({'mode': 'verified_existing_byte_identical' if args.check_existing else 'reconstructed',
                      'counts': counts, 'raw_manifest_records': 104, 'conflicting_tips': conflicts,
                      'output_files': [str(args.outdir / name) for name in rendered]}, ensure_ascii=False))
    if conflicts:
        raise SystemExit(2)


if __name__ == '__main__':
    main()
