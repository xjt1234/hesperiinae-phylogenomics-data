#!/usr/bin/env python3
"""BSM-only bounded supervisor. Default is read-only preflight; --run is explicit."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone

JOB = Path(__file__).resolve().parents[1]
CONTRACT = JOB / '02_config/conditional_completion_contract_20260908.json'
CONTRACT_SHA = 'f3336caf6d5d398062c128d832f9a520a937cd7cf1995746dac9b7a19a9ff714'
USER_DECISION = JOB / '09_SUMMARIES/USER_DECISION_20260908_NO_R5_CONTINUE_BSM.md'
USER_DECISION_SHA = '851d4fd92e13569633f75fe8b4347ca1762e30c18cd31cb2069e6bec8750d889'
RUNNER = JOB / '03_scripts/69_run_conditional_native_bsm_20260908.R'
WRAPPER = JOB / '03_scripts/r44_env.sh'
OUT = JOB / '08_qa/conditional_completion_20260908'
BSM = JOB / '05_bsm/M1_BSM_CONDITIONAL_20260908'
UNIT = 'bgb-r24-completion-20260908.service'
THREADS = ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS',
           'VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS')
STAGES = (('pilot',1,10,False),('run',100,200,True),('run',500,1000,True))
WALLTIME = 172800

def now():
    return datetime.now(timezone.utc).isoformat()

def log_size(path):
    try:return path.stat().st_size
    except FileNotFoundError:return 0  # A valid atomic work-to-commit move can race a snapshot.

def sha(path):
    if path.is_symlink() or not path.is_file():
        raise ValueError(f'Required regular, non-symlink file: {path}')
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()

def exact_path(path):
    if path.resolve(strict=False) != path:
        raise ValueError(f'Unexpected symlink/path redirection: {path}')
    return path

def write_json(path, value, replace=False):
    temp = path.with_name(path.name + '.tmp')
    with temp.open('x') as handle:
        json.dump(value,handle,indent=2,sort_keys=True);handle.write('\n')
        handle.flush();os.fsync(handle.fileno())
    if path.exists() and not replace:
        raise ValueError(f'Refusing overwrite: {path}')
    os.replace(temp,path)

def checked_inputs(runner_sha, fresh=True):
    if not re.fullmatch('[a-f0-9]{64}',runner_sha):
        raise ValueError('Expected caller-pinned runner SHA256')
    for path in (CONTRACT,RUNNER,WRAPPER,OUT,BSM):exact_path(path)
    if sha(CONTRACT)!=CONTRACT_SHA or sha(RUNNER)!=runner_sha:
        raise ValueError('Frozen contract/runner SHA256 mismatch')
    c=json.loads(CONTRACT.read_text())
    expected={'status':'APPROVED_CONDITIONAL_CONTINUATION','scientific_acceptance':'NONE',
              'backend':'native','model':'M1','KKT1':False,'pilot_successful_maps':1,
              'minimum_successful_maps':100,'first_stability_checkpoint':200,
              'max_successful_maps':500,'max_attempts':1000,'seed_base':202609080,
              'no_success_failure_stop':10,'maxtries':40000,'map_batch_size':100}
    if any(c.get(k)!=v for k,v in expected.items()):raise ValueError('Contract field mismatch')
    if c['resource_limits']['overall_walltime_seconds']!=WALLTIME or c['resource_limits']['bsm_bgb_cores']!=1:
        raise ValueError('BSM resource/time contract changed')
    pins={**c['frozen_inputs'],**c['frozen_helpers'],c['fit_rds']:c['fit_sha256'],c['postfit_rds']:c['postfit_sha256'],
          c['numerical_audit']['path']:c['numerical_audit']['sha256']}
    pins[str(USER_DECISION.relative_to(JOB))]=USER_DECISION_SHA
    for rel,digest in pins.items():
        path=exact_path(JOB/rel)
        if not path.is_relative_to(JOB) or sha(path)!=digest:raise ValueError(f'Source pin mismatch: {rel}')
    if fresh and (OUT.exists() or BSM.exists()):raise ValueError('New supervisor and BSM namespaces must not exist')
    return c,pins

def command(stage):
    action,target,attempts,resume=stage
    return ['bash',str(WRAPPER),'--file='+str(RUNNER),'--args','--action',action,
            '--contract-json',str(CONTRACT),'--contract-sha',CONTRACT_SHA,'--outdir',str(BSM),
            '--target-maps',str(target),'--max-attempts',str(attempts),'--resume',str(resume).lower()]

def resource_guard():
    rel=[x.split(':',2)[2] for x in Path('/proc/self/cgroup').read_text().splitlines() if x.startswith('0::')]
    if len(rel)!=1 or not rel[0].endswith('/'+UNIT):raise ValueError('Must run in the dedicated completion systemd unit')
    base=Path('/sys/fs/cgroup')/rel[0].lstrip('/')
    limits={k:(base/k).read_text().strip() for k in ('memory.max','memory.swap.max','pids.max','cpu.max')}
    quota,period=limits['cpu.max'].split()
    if not (0<int(limits['memory.max'])<=180000000000 and limits['memory.swap.max']=='0'
            and 0<int(limits['pids.max'])<=48 and quota!='max' and 0<int(quota)/int(period)<=40):
        raise ValueError(f'Unsafe resource limits: {limits}')
    return base,limits

def transactions(final=False):
    rows=[];commits=BSM/'commits'
    if commits.exists():
        for path in sorted(commits.iterdir()):
            exact_path(path)
            if not path.is_dir() or not re.fullmatch('attempt_[0-9]{6}',path.name):raise ValueError('Unknown commit entry')
            seal=path/'commit_seal.rds';sha(seal)
            records=[p for p in (path/'record.tsv',path/'recovery_record.tsv') if p.exists()]
            if len(records)!=1:raise ValueError('Ambiguous committed attempt record')
            with records[0].open() as handle:record=list(csv.DictReader(handle,delimiter='\t'))
            if len(record)!=1:raise ValueError('Attempt record must contain one row')
            r=record[0];aid=int(r['attempt_id'])
            if path.name!=f'attempt_{aid:06d}' or r['status'] not in ('SUCCESS','FAILED'):raise ValueError('Commit ID/status mismatch')
            if r['status']=='SUCCESS':
                mid=int(r['map_id']);mapsha=sha(path/'map.rds')
                if final and sha(BSM/'maps'/f'map_{mid:04d}.rds')!=mapsha:raise ValueError('Map projection mismatch')
            if final:
                view=BSM/'attempts'/f'attempt_{aid:06d}_{r["status"]}.tsv'
                if sha(view)!=sha(records[0]):raise ValueError('Attempt projection mismatch')
            rows.append(r)
    if [int(r['attempt_id']) for r in rows]!=list(range(1,len(rows)+1)):raise ValueError('Noncontiguous committed attempts')
    success=[int(r['map_id']) for r in rows if r['status']=='SUCCESS']
    if success!=list(range(1,len(success)+1)):raise ValueError('Noncontiguous successful map IDs')
    return {'committed_transactions':len(rows),'successful_transactions':len(success),'failed_transactions':len(rows)-len(success)}

def validate_result(result,stage,pid,code,c,runner_sha,counts):
    action,target,attempts,resume=stage
    expected={'scientific_acceptance':'NONE','invocation_pid':pid,'action':action,'resume':resume,
              'requested_target':target,'invocation_attempt_ceiling':attempts,'maximum_successful_maps':500,
              'global_maximum_attempts':1000,'ledger_pass':True,'conditional_contract_sha256':CONTRACT_SHA,
              'runner_sha256':runner_sha,'fit_sha256':c['fit_sha256'],'postfit_sha256':c['postfit_sha256'],
              'outdir':str(BSM),'no_success_failure_stop':10}
    if any(result.get(k)!=v for k,v in expected.items()):raise ValueError('Current invocation identity/ledger mismatch')
    if result.get('downstream_actions_launched') not in ([],{}):raise ValueError('Unexpected downstream action')
    for k,ck in [('successful','successful_transactions'),('attempted','committed_transactions'),('failed','failed_transactions')]:
        if type(result.get(k)) is not int or result[k]!=counts[ck]:raise ValueError('Result/committed transaction count mismatch')
    if not (0<=result['successful']<=500 and result['attempted']<=attempts):raise ValueError('Attempt/map cap exceeded')
    if not result.get('invocation_started') or not result.get('finished'):raise ValueError('Missing invocation timestamps')
    status=result.get('status');success=result['successful']
    if code==0 and target==1 and success==1 and status=='PILOT_COMPLETE_CONDITIONAL_ONLY':return 'CONTINUE'
    if code==0 and target==100 and success==100 and status=='TARGET_REACHED_NEEDS_EXTENSION':return 'CONTINUE'
    if code==0 and target==500 and success in (200,300,400,500) and status=='MONTE_CARLO_STABLE_CONDITIONAL_ONLY':return 'COMPUTE_COMPLETE_CONDITIONAL_REVIEW_REQUIRED'
    if code!=0 and target==500 and success==500 and status=='INCOMPLETE_AT_500':return 'INCOMPLETE_CONDITIONAL_REVIEW_REQUIRED'
    raise ValueError(f'Stage did not complete as authorized: exit={code}, status={status}, successful={success}')

def selftest(c,runner_sha):
    # Pure in-memory fake status records; no subprocesses, simulations or files.
    def fixture(stage,success,status):
        action,target,attempts,resume=stage
        return dict(scientific_acceptance='NONE',invocation_pid=123,action=action,resume=resume,requested_target=target,
          invocation_attempt_ceiling=attempts,maximum_successful_maps=500,global_maximum_attempts=1000,ledger_pass=True,
          conditional_contract_sha256=CONTRACT_SHA,runner_sha256=runner_sha,fit_sha256=c['fit_sha256'],postfit_sha256=c['postfit_sha256'],
          outdir=str(BSM),no_success_failure_stop=10,downstream_actions_launched=[],successful=success,attempted=success,
          failed=0,invocation_started='2026-09-08T00:00:00Z',finished='2026-09-08T01:00:00Z',status=status)
    def check(r,stage,code):
        return validate_result(r,stage,123,code,c,runner_sha,dict(committed_transactions=r['attempted'],successful_transactions=r['successful'],failed_transactions=r['failed']))
    cases=[(STAGES[0],1,'PILOT_COMPLETE_CONDITIONAL_ONLY',0,'CONTINUE'),
           (STAGES[1],100,'TARGET_REACHED_NEEDS_EXTENSION',0,'CONTINUE'),
           (STAGES[2],200,'MONTE_CARLO_STABLE_CONDITIONAL_ONLY',0,'COMPUTE_COMPLETE_CONDITIONAL_REVIEW_REQUIRED'),
           (STAGES[2],500,'INCOMPLETE_AT_500',1,'INCOMPLETE_CONDITIONAL_REVIEW_REQUIRED')]
    for stage,n,status,code,wanted in cases:assert check(fixture(stage,n,status),stage,code)==wanted
    bad={'invocation_pid':124,'ledger_pass':False,'conditional_contract_sha256':'0'*64,'resume':True,
         'outdir':str(JOB),'runner_sha256':'0'*64,'fit_sha256':'0'*64,'scientific_acceptance':'ACCEPTED'}
    for key,value in bad.items():
        r=fixture(STAGES[0],1,'PILOT_COMPLETE_CONDITIONAL_ONLY');r[key]=value
        try:check(r,STAGES[0],0)
        except ValueError:pass
        else:raise AssertionError('Must reject '+key)
    try:check(fixture(STAGES[0],1,'PILOT_COMPLETE_CONDITIONAL_ONLY'),STAGES[0],1)
    except ValueError:pass
    else:raise AssertionError('Nonzero exit accepted')
    return {'status':'PASS','positive_tests':4,'negative_tests':9,'files_created':0,'subprocesses_launched':0}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--runner-sha',required=True)
    mode=parser.add_mutually_exclusive_group();mode.add_argument('--run',action='store_true');mode.add_argument('--selftest',action='store_true')
    args=parser.parse_args();c,pins=checked_inputs(args.runner_sha)
    if args.selftest:print(json.dumps(selftest(c,args.runner_sha),indent=2));return 0
    launch_cmd=['systemd-run','--user','--unit='+UNIT,'--property=CPUQuota=4000%','--property=TasksMax=48',
      '--property=MemoryMax=180000000000','--property=MemorySwapMax=0','--property=RuntimeMaxSec=48h',
      '--property=KillMode=control-group','--property=Restart=no','--property=RemainAfterExit=yes',
      *['--setenv='+k+'=1' for k in THREADS],sys.executable,str(Path(__file__).resolve()),'--runner-sha',args.runner_sha,'--run']
    if not args.run:
        print(json.dumps({'status':'READ_ONLY_PREFLIGHT_PASS','scientific_acceptance':'NONE','files_created':0,
          'computations_launched':0,'scope':'BSM only; no R5/R4/optimization/diagnostics','stages':[command(s) for s in STAGES],
          'suggested_systemd_command_not_executed':shlex.join(launch_cmd),'contract_sha256':CONTRACT_SHA,
          'latest_user_scope_decision_sha256':USER_DECISION_SHA},indent=2));return 0
    cg,limits=resource_guard();env=os.environ.copy();env.update({k:'1' for k in THREADS})
    OUT.mkdir();start=time.monotonic();deadline=start+WALLTIME;completed=[];child=None
    write_json(OUT/'launch.json',{'started_utc':now(),'supervisor_pid':os.getpid(),'supervisor_sha256':sha(Path(__file__)),
      'runner_sha256':args.runner_sha,'contract_sha256':CONTRACT_SHA,'source_pins':pins,'resource_limits':limits,
      'scientific_acceptance':'NONE','scope':'BSM_ONLY','retry_policy':'NONE','deadline_seconds':WALLTIME})
    def interrupted(signum,frame):raise RuntimeError(f'Supervisor interrupted by signal {signum}')
    signal.signal(signal.SIGTERM,interrupted);signal.signal(signal.SIGINT,interrupted)
    last_progress=now();previous=None
    try:
        for index,stage in enumerate(STAGES,1):
            checked_inputs(args.runner_sha,fresh=False)
            if time.monotonic()>=deadline:raise TimeoutError('48-hour total deadline reached before next stage')
            before={p.name for p in (BSM/'invocations').glob('status_*.json')}
            log=OUT/f'stage_{index:02d}.log';stage_start=time.monotonic();started=now()
            with log.open('x') as handle:
                child=subprocess.Popen(command(stage),cwd=JOB,stdin=subprocess.DEVNULL,stdout=handle,stderr=subprocess.STDOUT,env=env,start_new_session=True)
                write_json(OUT/f'stage_{index:02d}.launch.json',{'pid':child.pid,'started_utc':started,'argv':command(stage),'runner_sha256':args.runner_sha})
                while True:
                    code=child.poll();counts=transactions()
                    native_logs=list((BSM/'work').glob('attempt_*/stdout.log'))+list((BSM/'commits').glob('attempt_*/stdout.log'))
                    log_bytes=log.stat().st_size;native_bytes=sum(log_size(p) for p in native_logs)
                    observation=(*counts.values(),log_bytes,native_bytes)
                    progress=previous is not None and any(a>b for a,b in zip(observation,previous))
                    if progress:last_progress=now()
                    previous=observation
                    state={'status':'RUNNING_CONDITIONAL_BSM','timestamp_utc':now(),'stage_index':index,'child_pid':child.pid,
                      'child_exit_code':code,'elapsed_seconds':time.monotonic()-start,'stage_elapsed_seconds':time.monotonic()-stage_start,
                      **counts,'child_log_bytes':log_bytes,'native_attempt_log_bytes':native_bytes,'progress_since_previous':progress,
                      'last_observed_progress_utc':last_progress,'progress_definition':'Committed attempts/maps or log byte growth, not process existence',
                      'scientific_acceptance':'NONE','original_KKT1':False,'completed_stages':completed}
                    state['resources']={k:(cg/k).read_text().strip() for k in ('memory.current','memory.peak','memory.events','pids.current','pids.events')}
                    write_json(OUT/'STATUS.json',state,replace=True)
                    with (OUT/'heartbeat.jsonl').open('a') as hb:hb.write(json.dumps(state,sort_keys=True)+'\n')
                    if code is not None:break
                    if time.monotonic()>=deadline:raise TimeoutError('48-hour total deadline reached')
                    time.sleep(min(30,max(0.1,deadline-time.monotonic())))
            counts=transactions(final=True)
            statuses=[p for p in (BSM/'invocations').glob(f'status_*_{child.pid}.json') if p.name not in before]
            if len(statuses)!=1:raise ValueError('Expected exactly one new current-PID invocation status')
            result=json.loads(statuses[0].read_text());outcome=validate_result(result,stage,child.pid,code,c,args.runner_sha,counts)
            record={'stage_index':index,'pid':child.pid,'exit_code':code,'finished_utc':now(),'status_path':str(statuses[0]),
                    'status_sha256':sha(statuses[0]),'log_sha256':sha(log),'result':result,'validated_outcome':outcome}
            write_json(OUT/f'stage_{index:02d}.exit.json',record);completed.append(record);child=None
            if outcome!='CONTINUE':
                write_json(OUT/'STATUS.json',{'status':outcome,'timestamp_utc':now(),'scientific_acceptance':'NONE',
                  'original_KKT1':False,'completed_stages':completed,'downstream_actions_launched':[]},replace=True)
                return 0 if outcome=='COMPUTE_COMPLETE_CONDITIONAL_REVIEW_REQUIRED' else 1
        raise ValueError('Stage schedule ended without terminal result')
    except BaseException as exc:
        if child is not None and child.poll() is None:
            try:os.killpg(child.pid,signal.SIGTERM)
            except ProcessLookupError:pass
            try:child.wait(timeout=10)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait(timeout=10)
        failure={'status':'INCOMPLETE_TIMEOUT_REVIEW_REQUIRED' if isinstance(exc,TimeoutError) else 'FAILED_CONDITIONAL_REVIEW_REQUIRED',
          'timestamp_utc':now(),'error':repr(exc),'scientific_acceptance':'NONE','original_KKT1':False,
          'completed_stages':completed,'last_child_pid':None if child is None else child.pid,
          'last_child_exit_code':None if child is None else child.poll(),'retry_policy':'NONE','downstream_actions_launched':[]}
        write_json(OUT/'controller_error.json',failure);write_json(OUT/'STATUS.json',failure,replace=True)
        raise

if __name__=='__main__':
    try:sys.exit(main())
    except Exception as exc:print(f'ERROR: {exc}',file=sys.stderr);sys.exit(1)
