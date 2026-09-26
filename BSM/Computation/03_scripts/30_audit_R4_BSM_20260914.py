#!/usr/bin/env python3
"""Read-only gates and independent flat-record audit for the new R4 BSM only."""
from __future__ import annotations
import csv
import hashlib
import itertools
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

JOB = Path(__file__).resolve().parents[1]
RUN_REL = '05_bsm/M1_R4_BSM_BRIDGE_20260914'
CONTROL_REL = '08_qa/conditional_bridge_R4_20260914'
CONTRACT_REL = '02_config/conditional_bridge_bsm_R4_20260914.json'
AREAS = ['AF','AUS','CAM','ENA','EPA','IND','MDG','ORI','SAM','WNA','WPA']
BOUNDS = [0.,3.3,13.9,23.03,33.9,45.]
EPOCHS = ['0-3.3','3.3-13.9','13.9-23.03','23.03-33.9','33.9-45']
FIT_SHA = 'bb855fe8bfc2919d5eed94224b177463d517dc78fbe958cc34faf44fcca1bff4'
POSTFIT_SHA = '566b86c007d123be0fccf1064abf076ad26e6c082ceac8133cbed447ef52e08d'
GEO_SHA = '3aa227e4600629eb77178ab68a30bb7cf408189bbcb9f8989ae1f184ce834ca4'
MANIFEST_SHA = 'f7e121a3ce4aff304767eef6df39eba6f9ef473242aad48411a3b2c204d4a306'
FORCED = ['Error_in_stochastic_simulation; no success on branch', 'manually devising a history to force-fit', 'program is in the manually-sorted events section', 'Attempted manual history sampled disallowed state']

class NotReady(RuntimeError):
    """Expected when a fresh finite BSM is not at its terminal checkpoint."""

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), f'Expected regular non-symlink file: {path}')
    before = path.stat()
    with path.open('rb') as handle:
        h = hashlib.sha256()
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    after = path.stat()
    require((before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns), f'File changed during audit: {path}')
    return h.hexdigest()

def rows(path):
    with Path(path).open(encoding='utf-8',newline='') as handle:
        yield from csv.DictReader(handle,delimiter='\t')

def one(path):
    values=list(rows(path)); require(len(values)==1,f'Expected single row: {path}')
    return values[0]

def number(row,key):
    value=float(row[key]); require(math.isfinite(value),f'Nonfinite {key}')
    return value

def integer(row,key):
    value=number(row,key); require(value==int(value),f'Noninteger {key}')
    return int(value)

def boolean(value):
    require(str(value).upper() in ('TRUE','FALSE'),f'Invalid boolean: {value}')
    return str(value).upper()=='TRUE'

def epoch(age):
    require(math.isfinite(age) and 0<=age<45,'Event age outside declared strata')
    return next(EPOCHS[i] for i in range(5) if BOUNDS[i]<=age<BOUNDS[i+1])

def context():
    control_path=JOB/CONTROL_REL/'STATUS.json'
    if not control_path.is_file():
        raise NotReady('BSM controller terminal status does not exist; no final analysis is permitted.')
    control=json.loads(control_path.read_text())
    if control.get('status')!='COMPUTE_COMPLETE_CONDITIONAL_REVIEW_REQUIRED':
        raise NotReady('BSM is not compute-complete: '+str(control.get('status')))
    require(control.get('scientific_acceptance')=='NONE' and control.get('original_KKT1') is False,'Conditional status changed')
    launch=json.loads((JOB/CONTROL_REL/'launch.json').read_text())
    contract_path=JOB/CONTRACT_REL; contract=json.loads(contract_path.read_text())
    require(sha(contract_path)==launch['contract_sha256'],'Contract changed after launch')
    require(contract['model']=='M1' and contract['j']==0 and contract['w']==1 and contract['KKT1'] is False and contract['scientific_acceptance']=='NONE','Model or diagnostic scope changed')
    require(contract['n_tips']==417 and contract['area_order']==AREAS and contract['seed_base']==202609140,'Input dimensions/seed scope changed')
    require(contract['fit_sha256']==FIT_SHA and contract['postfit_sha256']==POSTFIT_SHA and contract['input_manifest_sha256']==MANIFEST_SHA,'Unexpected fit/postfit/input identity')
    require(not contract.get('bridge_native_cache'),'Historical preparation cache is forbidden')
    pins={**contract['frozen_inputs'],**contract['frozen_helpers'],contract['fit_rds']:FIT_SHA,contract['postfit_rds']:POSTFIT_SHA,**launch['source_pins']}
    pins['01_inputs/frozen/geog_scenarioA417_analysis_order.LagrangePHYLIP']=GEO_SHA
    pins['01_inputs/frozen/input_manifest.json']=MANIFEST_SHA
    for rel,expected in pins.items():
        path=JOB/rel
        require(path.resolve().is_relative_to(JOB) and path.resolve()==path,'Source escapes current R4 root')
        require(sha(path)==expected,'Frozen file hash mismatch: '+rel)
    stages=control['completed_stages']; require(len(stages)>=2,'Pilot and production stage evidence missing')
    for stage in stages:
        require(stage['exit_code']==0,'Nonzero completed stage')
        status_path=Path(stage['status_path'])
        require(status_path.resolve().is_relative_to(JOB/RUN_REL/'invocations'),'Invocation path escapes new run')
        require(sha(status_path)==stage['status_sha256'],'Invocation status changed')
        require(sha(JOB/CONTROL_REL/f"stage_{stage['stage_index']:02d}.log")==stage['log_sha256'],'Stage log changed')
        require(json.loads(status_path.read_text())==stage['result'],'Embedded and on-disk result differ')
    final=stages[-1]['result']; n=int(final['successful'])
    require(final['status']=='MONTE_CARLO_STABLE_CONDITIONAL_ONLY' and final['ledger_pass'] is True,'Final outcome not stable with valid ledger')
    require(n in (200,300,400,500),'Invalid final checkpoint')
    require(final['failed']==0 and final['attempted']==n,'Failed/rejected attempts require manual review before automatic figure finalization; preserved BSM must not be rerun')
    require(final['conditional_contract_sha256']==launch['contract_sha256'] and final['fit_sha256']==FIT_SHA and final['postfit_sha256']==POSTFIT_SHA,'Final result source binding mismatch')
    run=JOB/RUN_REL; summary=run/'summaries'/f'n_{n:04d}'
    status=one(summary/'summary_status.tsv')
    require(integer(status,'n_maps')==n and integer(status,'successful')==n and integer(status,'failed')==0,'Summary population mismatch')
    require(status['status']=='STABLE' and boolean(status['overall_stable']),'Final summary not stable')
    for field in ('conservation_pass','native_totals_crosscheck_pass','seed_schedule_pass','attempt_ledger_pass'):
        require(boolean(status[field]),'Final summary gate failed: '+field)
    return dict(job=JOB,run=run,summary=summary,n=n,contract=contract,contract_sha=launch['contract_sha256'],control=control,source_pins=pins)

def audit(ctx):
    run,summary,n,c=ctx['run'],ctx['summary'],ctx['n'],ctx['contract']
    inventory={}; checks=[]; map_audits=[]
    def hashed(path):
        digest=sha(path); inventory[str(path.relative_to(JOB))]=digest; return digest
    def passed(name): checks.append(dict(check=name,status='PASS'))
    manifest=json.loads((run/'bsm_run_manifest.json').read_text())
    technical=json.loads((run/'technical_validation.json').read_text())
    identity=json.loads((run/'stochastic_mapping_inputs_identity.json').read_text())
    require(manifest['conditional_contract_sha256']==ctx['contract_sha'] and manifest['scientific_acceptance']=='NONE','Manifest contract/status mismatch')
    require(manifest['source_optimizer_fit_sha256']==FIT_SHA and manifest['postfit_ancestral_states_sha256']==POSTFIT_SHA and manifest['input_manifest_sha256']==MANIFEST_SHA,'Manifest identity mismatch')
    require(technical['immutable_source_KKT1'] is False and technical['states']==562 and technical['time_strata']==5 and technical['tips']==417 and technical['internal_nodes']==416,'Technical dimensions/diagnostic mismatch')
    require(abs(float(technical['lnL_delta']))<=1e-8 and float(technical['node_probability_row_error'])<=1e-8,'Marginal/likelihood technical gate failed')
    require(identity['cache_origin']=='FRESH_NATIVE_PREPARATION_FROM_NEW_R4_POSTFIT_ONLY','Cache not newly prepared')
    require(identity['fit_sha256']==FIT_SHA and identity['postfit_sha256']==POSTFIT_SHA and identity['input_manifest_sha256']==MANIFEST_SHA,'Cache scientific ownership mismatch')
    require(identity['technical_validation_sha256']==manifest['technical_validation_sha256'] and hashed(run/'stochastic_mapping_inputs.rds')==identity['sha256'],'Cache digest/technical ownership mismatch')
    for name in ('bsm_run_manifest.json','technical_validation.json','stochastic_mapping_inputs_identity.json'):
        hashed(run/name)
    state_path=JOB/'04_runs/M1_R4_postfit_v1/state_dictionary.tsv'
    states={integer(r,'state_index'):r['internal_state'].replace('+','') for r in rows(state_path)}
    expected=['NULL']+[''.join(x) for k in range(1,5) for x in itertools.combinations('ABCDEFGHIJK',k)]
    require([states[i] for i in range(1,563)]==expected,'State dictionary changed')
    fit=one(JOB/'04_runs/M1_R4_final_v1/fit_summary.tsv')
    require(integer(fit,'max_range')==4 and integer(fit,'j')==0 and fit['kkt1']=='FALSE','Fit scope mismatch')
    require(abs(number(fit,'d')-c['d'])<1e-14 and abs(number(fit,'e')-c['e'])<1e-14,'Fit parameter mismatch')
    passed('new_frozen_fit_geography_marginals_and_cache')
    commits=sorted((run/'commits').iterdir()); maps=sorted((run/'maps').iterdir()); ledgers=sorted((run/'attempts').iterdir())
    require([x.name for x in commits]==[f'attempt_{i:06d}' for i in range(1,n+1)],'Commit sequence mismatch')
    require([x.name for x in maps]==[f'map_{i:04d}.rds' for i in range(1,n+1)],'Map sequence mismatch')
    require([x.name for x in ledgers]==[f'attempt_{i:06d}_SUCCESS.tsv' for i in range(1,n+1)],'Ledger sequence mismatch')
    require(not any((run/'work').iterdir()),'Uncommitted attempt work remains')
    jump_counts={}; max_tail=0.; max_union=0.; qhash=set()
    for i,commit in enumerate(commits,1):
        rec=one(commit/'record.tsv')
        require(integer(rec,'attempt_id')==integer(rec,'map_id')==i and rec['status']=='SUCCESS','Map/attempt identity mismatch')
        require(integer(rec,'warnings')==0 and rec['error']=='','Accepted warning/error')
        require(integer(rec,'map_seed')==c['seed_base']+i and integer(rec,'source_seed')==c['seed_base']+1000000+i,'Seed schedule mismatch')
        require(hashed(commit/'record.tsv')==hashed(ledgers[i-1]),'Ledger projection mismatch')
        map_sha=hashed(commit/'map.rds'); require(map_sha==hashed(maps[i-1]),'Map projection mismatch')
        for f in ('commit_seal.rds','raw_native_map.rds','timing_and_warnings.json','stdout.log'):
            hashed(commit/f)
        timing=json.loads((commit/'timing_and_warnings.json').read_text())
        require(timing['attempt_id']==i and timing['forced_history_detected'] is False and timing['warnings']==[],'Accepted warning/force-fit evidence')
        require(not any(marker in (commit/'stdout.log').read_text() for marker in FORCED),'Forced-history marker in accepted log')
        calls=0; real_total=0; tail_sum=0.
        for calls,row in enumerate(rows(commit/'bridge_branch_audit.tsv'),1):
            require(integer(row,'call_id')==calls and integer(row,'n_states')==562,'Bridge call/state mismatch')
            tail=number(row,'relative_tail_bound'); require(0<=tail<=c['bridge_relative_tail_tolerance'],'Bridge tail tolerance exceeded')
            a,b=integer(row,'a'),integer(row,'b'); require(1<=a<=562 and 1<=b<=562,'Bridge endpoint outside dictionary')
            require(row['start_range'].replace('_','NULL')==states[a] and row['end_range'].replace('_','NULL')==states[b],'Bridge state semantics mismatch')
            t,old,young=(number(row,k) for k in ('t','older_endpoint_age','younger_endpoint_age'))
            require(t>0 and young>=-1e-8 and old<=45+1e-8 and abs(old-young-t)<1e-8,'Bridge time inconsistency')
            require(1<=integer(row,'master_node')<=833 and integer(row,'master_node')!=418,'Invalid bridge master node')
            real,virtual,loops=(integer(row,k) for k in ('real_jumps','virtual_events','self_loops'))
            distance=len((set() if a==1 else set(states[a])).symmetric_difference(set() if b==1 else set(states[b])))
            require(min(real,virtual,loops)>=0 and virtual==real+loops and real>=distance and (real-distance)%2==0,'Impossible bridge counts/endpoints')
            require(integer(row,'last_series_term')<=10000 and number(row,'lambda')<=1000,'Bridge finite bound exceeded')
            require(row['Q_sha256']==row['Qeffective_sha256'],'Q effective identity mismatch')
            qhash.add(row['Qeffective_sha256']); max_tail=max(max_tail,tail); tail_sum+=tail; real_total+=real
        require(calls==c['bridge_expected_branch_calls'],'Bridge branch coverage mismatch')
        hashed(commit/'bridge_branch_audit.tsv'); jump_counts[i]=real_total; max_union=max(max_union,tail_sum)
        map_audits.append(dict(map_id=i,map_sha256=map_sha,branch_calls=calls,real_jumps=real_total,tail_sum=tail_sum,elapsed_seconds=timing['elapsed_seconds']))
    passed('all_commits_maps_ledger_seeds_warnings_and_branch_records')
    counts=Counter(); periods=Counter(); routes=Counter(); ext=Counter(); event_ids=defaultdict(list); code=dict(zip('ABCDEFGHIJK',AREAS))
    for row in rows(summary/'anagenetic_events_long.tsv'):
        i=integer(row,'map_id'); require(1<=i<=n,'Event map id outside population'); event_ids[i].append(integer(row,'event_id'))
        typ=row['event_type']; ep=epoch(number(row,'age_ma')); require(ep==row['epoch'] and typ in ('d','e'),'Event type/epoch mismatch')
        old,new=set(row['current_range']),set(row['new_range'])
        require(old<=set(code) and new<=set(code) and 1<=len(old)<=4 and 1<=len(new)<=4,'Event range outside R4')
        gain,lost=new-old,old-new
        if typ=='d':
            require(len(gain)==1 and not lost and row['to_code'] in gain and row['from_code'] in old,'Expansion/source semantics invalid')
            require(row['from_area']==code[row['from_code']] and row['to_area']==code[row['to_code']],'Direction label mismatch')
            routes[(i,row['from_code'],row['to_code'])]+=1
        else:
            require(len(lost)==1 and not gain and row['affected_code'] in lost,'Contraction semantics invalid'); ext[(i,row['affected_code'])]+=1
        require(row['affected_area']==code[row['affected_code']],'Affected-area mismatch'); counts[(i,typ)]+=1; periods[(i,ep,typ)]+=1
    for i in range(1,n+1):
        require(event_ids[i]==list(range(1,jump_counts[i]+1)),'Long-event IDs/bridge counts differ')
    clado=defaultdict(list)
    for row in rows(summary/'cladogenetic_events_long.tsv'):
        i=integer(row,'map_id'); require(1<=i<=n and epoch(number(row,'age_ma'))==row['epoch'] and '(j)' not in row['event_type'],'Cladogenetic identity/epoch/j mismatch')
        clado[i].append(integer(row,'node'))
    for i in range(1,n+1): require(sorted(clado[i])==list(range(418,834)),'Expected 416 true unique cladogenetic nodes')
    for name,keys,expected,count_expected in [('dispersal_routes_by_map.tsv',('map_id','from_code','to_code'),routes,n*110),('extinction_by_area_by_map.tsv',('map_id','area_code'),ext,n*11),('per_map_period_counts.tsv',('map_id','epoch','event'),periods,n*10)]:
        seen=set()
        for row in rows(summary/name):
            key=(integer(row,'map_id'),)+tuple(row[k] for k in keys[1:]); require(key not in seen,'Duplicate summary grid row'); seen.add(key)
            require(integer(row,'count')==expected[key],'Independent event recount differs: '+name)
        require(len(seen)==count_expected,'Incomplete summary grid: '+name)
    pm=list(rows(summary/'per_map_counts.tsv')); require([integer(r,'map_id') for r in pm]==list(range(1,n+1)),'Per-map sequence mismatch')
    for row in pm:
        i=integer(row,'map_id'); require(integer(row,'d')==counts[(i,'d')] and integer(row,'e')==counts[(i,'e')] and integer(row,'j')==0 and integer(row,'cladogenetic_events')==416,'Per-map count invariant failed')
    passed('independent_event_semantics_counts_grids_and_cladogenesis')
    conv=list(rows(summary/'bsm_convergence.tsv')); require([integer(r,'n_maps') for r in conv]==list(range(100,n+1,100)),'Checkpoint sequence mismatch')
    require([boolean(r['overall_stable']) for r in conv]==[False]*(len(conv)-1)+[True],'Not the first passing checkpoint')
    require(all(boolean(v) for k,v in conv[-1].items() if k.endswith('_pass') or k.endswith('_stable')) and number(conv[-1],'metric_threshold')==.05,'Prespecified stopping gates failed')
    require(boolean(one(summary/'attempt_counts_audit.tsv')['all_pass']),'Producer attempt ledger QA failed')
    for path in sorted(summary.iterdir()): require(path.is_file(),'Unexpected summary directory'); hashed(path)
    for rel in ctx['source_pins']: hashed(JOB/rel)
    hashed(JOB/CONTRACT_REL); hashed(JOB/CONTROL_REL/'STATUS.json'); hashed(JOB/CONTROL_REL/'launch.json')
    passed('first_passing_prespecified_checkpoint_and_input_hash_inventory')
    report=dict(status='TECHNICAL_PASS_CONDITIONAL_ONLY',scientific_acceptance='NONE',original_KKT1=False,n_maps=n,failed=0,forced_histories=0,warnings=0,source_fit_sha256=FIT_SHA,postfit_sha256=POSTFIT_SHA,contract_sha256=ctx['contract_sha'],total_d=sum(v for (i,t),v in counts.items() if t=='d'),total_e=sum(v for (i,t),v in counts.items() if t=='e'),bridge_rows_checked=n*c['bridge_expected_branch_calls'],maximum_relative_tail=max_tail,maximum_per_map_tail_sum=max_union,effective_Q_hashes=sorted(qhash),all_checks_pass=True,checked_utc=datetime.now(timezone.utc).isoformat(),new_histories_generated=0,new_model_fits=0,limitations=['Conditional on fixed tree, geography, M1, R4, fitted parameters and model-weighted source allocation.','KKT1 FALSE is retained; Monte Carlo stability is not global optimization or model robustness.','Python hashes opaque RDS seals and raw maps; the producer decodes and validates seals and map payloads on every rebuild. This independent audit does not reconstruct every raw branch event chain.','Rejected attempts are preserved, and any nonzero failure count requires manual review before this automatic finalizer.','R6 is a separate joint geography/range scenario, not a pure range-cap comparison; no R6 data are included.'])
    return dict(report=report,checks=checks,map_audits=map_audits,inventory=inventory)

def write_audit(result,out):
    out.mkdir(parents=True,exist_ok=False)
    for name,key in [('checks.tsv','checks'),('map_audit.tsv','map_audits')]:
        values=result[key]
        with (out/name).open('x',newline='') as handle:
            writer=csv.DictWriter(handle,fieldnames=list(values[0]),delimiter='\t'); writer.writeheader(); writer.writerows(values)
    (out/'audit_report.json').write_text(json.dumps(result['report'],indent=2)+'\n')
    (out/'artifact_sha256.json').write_text(json.dumps(result['inventory'],indent=2)+'\n')

if __name__=='__main__':
    print(json.dumps(audit(context())['report'],indent=2))
