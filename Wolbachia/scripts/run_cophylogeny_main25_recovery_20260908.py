#!/usr/bin/env python3
"""Recovery with the audited existing user R library; original failed attempt preserved."""
import argparse, concurrent.futures, datetime, fcntl, hashlib, json, os, subprocess, threading, time
from pathlib import Path
PLAN_SHA = '96fc5349d49f7cb73d758a6c094d8c815f34181098e40ac60a0e369cc5ca6094'
R_SHA = 'df207f4e1cf052f2a2e86490f272aa73689f2b68eb92e80f103b5c7578a7c982'
R_EXE = '/home/data/t200301/miniconda3/envs/woltree/bin/Rscript'
USER_R_LIBRARY = '/home/data/t200301/R/x86_64-pc-linux-gnu-library/4.4'
SITE_R_LIBRARY = '/home/data/t200301/miniconda3/envs/woltree/lib/R/library'

def r_environment():
    env=os.environ.copy();env.pop('R_HOME',None)
    env.update({k:'1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS']})
    env.update({'R_LIBS_USER':USER_R_LIBRARY,'R_LIBS_SITE':SITE_R_LIBRARY})
    return env

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(1048576), b''): h.update(b)
    return h.hexdigest()

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def write_new(path, value):
    with path.open('x') as f: json.dump(value, f, ensure_ascii=False, indent=2)
def assert_outputs(cp):
    assert cp['status']=='complete'
    assert cp['outputs'], 'Completed record has no output hashes'
    for rec in cp['outputs']:
        assert Path(rec['path']).is_absolute(), rec['path']
        assert sha(rec['path'])==rec['sha256'], rec['path']

def assert_r_completion(path, sid):
    inner=json.loads(path.read_text());assert_outputs(inner)
    assert inner['data_kind']=='real' and inner['scenario_count']==1
    assert inner['scenarios'][0]['scenario_id']==sid and inner['scenarios'][0]['nperm']==9999
    assert inner['source_hashes'], 'R completion has no frozen source hashes'
    for rec in inner['source_hashes']:
        assert Path(rec['path']).is_absolute(), rec['path']
        assert sha(rec['path'])==rec['sha256'], 'Frozen source changed: '+rec['path']
    return inner

class Driver:
    def __init__(self,args):
        self.run=args.run_dir.resolve();self.threads=args.threads
        self.plan_path=self.run/'01_PROVENANCE/cophylogeny_main25_16scenarios_20260908.json'
        self.r=self.run/'scripts/frozen/cophylogeny_exact_inputs_R26_df207f4e1cf052f2.R'
        assert sha(self.plan_path)==PLAN_SHA and sha(self.r)==R_SHA
        self.plan=json.loads(self.plan_path.read_text());assert len(self.plan['scenarios'])==16
        self.root=self.run/'04_COPHYLOGENY/main25_R26_recovery_20260908'
        self.root.mkdir(exist_ok=True)
        for name in ['inputs','analyses','logs']: (self.root/name).mkdir(exist_ok=True)
        self.config={'plan_sha256':PLAN_SHA,'R_sha256':R_SHA,'wrapper_sha256':sha(__file__),'max_concurrent_R_processes':self.threads,'R_compute_threads_each':1,'R_LIBS_USER':USER_R_LIBRARY,'R_LIBS_SITE':SITE_R_LIBRARY}
        configfile=self.root/'execution_config.json'
        if configfile.exists(): assert json.loads(configfile.read_text())==self.config
        else: write_new(configfile,self.config)
        self.guard=threading.Lock();self.active={};self.complete=[];self.failure=None;self.failed_cases={};self.finished=False
        self.env=r_environment()

    def progress(self):
        with self.guard:
            data={'updated_utc':now(),'driver_pid':os.getpid(),'status':'failed' if self.failure else ('complete' if self.finished else ('running' if self.active else 'between_cases')),'active':dict(self.active),'completed_scenarios':list(self.complete),'planned_scenarios':16,'failure':self.failure,'failed_cases':dict(self.failed_cases),'config':self.config}
            dest=self.run/'checkpoints/cophylogeny_main25_recovery_progress.json'
            temp=dest.with_name(dest.name+'.tmp.'+str(os.getpid()))
            temp.write_text(json.dumps(data,indent=2)+'\n');os.replace(temp,dest)

    def prepare_case_manifests(self):
        for s in self.plan['scenarios']:
            manifest=dict(self.plan,scenarios=[s]);path=self.root/'inputs'/(s['id']+'.json')
            if path.exists():assert json.loads(path.read_text())==manifest
            else:write_new(path,manifest)

    def case(self,s):
        sid=s['id']
        try:
            return self._case(s)
        except BaseException as error:
            with self.guard:
                previous=self.active.pop(sid,None)
                self.failure=str(error)
                self.failed_cases[sid]={'status':'failed','failed_utc':now(),'reason':str(error),'exception_type':type(error).__name__,'log':str(self.root/'logs'/(sid+'.log')),'previous_process':previous}
            self.progress()
            raise

    def _case(self,s):
        sid=s['id'];manifest=self.root/'inputs'/(sid+'.json');out=self.root/'analyses'/sid
        cpfile=self.run/'checkpoints'/('cophylogeny_main25_recovery_case_'+str(s['scenario_index'])+'.complete.json')
        cmd=[R_EXE,'--vanilla',str(self.r),'--inputs-json',str(manifest),'--out-dir',str(out)]
        if cpfile.exists():
            cp=json.loads(cpfile.read_text());assert cp['config']==self.config and cp['command']==cmd
            assert_outputs(cp)
            innerfile=out/'completed.json'
            assert any(Path(rec['path'])==innerfile for rec in cp['outputs']), 'Case completion is not bound to R completion'
            inner=assert_r_completion(innerfile,sid)
            assert cp['summary']==inner['scenarios'][0], 'Case summary differs from R completion'
            with self.guard:self.complete.append(sid)
            self.progress();return cp
        assert not out.exists(), 'Incomplete prior case exists; preserve and review before separate recovery: '+str(out)
        assert sha(self.plan_path)==PLAN_SHA and sha(self.r)==R_SHA
        logpath=self.root/'logs'/(sid+'.log')
        with logpath.open('x') as log:
            p=subprocess.Popen(cmd,env=self.env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
            with self.guard:self.active[sid]={'R_pid':p.pid,'started_utc':now(),'log':str(logpath),'command':cmd}
            try:
                while True:
                    self.progress()
                    try:rc=p.wait(timeout=15);break
                    except subprocess.TimeoutExpired:pass
            except BaseException:
                p.terminate()
                try:p.wait(timeout=10)
                except subprocess.TimeoutExpired:p.kill();p.wait()
                raise
        assert rc==0, 'R case failed: '+sid+'; see '+str(logpath)
        innerfile=out/'completed.json';inner=assert_r_completion(innerfile,sid)
        assert sha(self.plan_path)==PLAN_SHA and sha(self.r)==R_SHA
        records=inner['outputs']+[{'path':str(innerfile),'sha256':sha(innerfile)},{'path':str(logpath),'sha256':sha(logpath)},{'path':str(manifest),'sha256':sha(manifest)}]
        cp={'status':'complete','scenario_id':sid,'finished_utc':now(),'config':self.config,'command':cmd,'R_exit_code':rc,'summary':inner['scenarios'][0],'outputs':records}
        write_new(cpfile,cp)
        with self.guard:self.active.pop(sid,None);self.complete.append(sid)
        self.progress();return cp

    def execute(self):
        existing=self.run/'checkpoints/cophylogeny_main25_recovery.complete.json'
        if existing.exists():
            done=json.loads(existing.read_text());assert done['config']==self.config;assert_outputs(done)
            assert done['scenario_count']==16 and len(done['scenarios'])==16
            assert {x['scenario_id'] for x in done['scenarios']}=={s['id'] for s in self.plan['scenarios']}
            checked=[]
            for s in self.plan['scenarios']:
                casecp=self.run/'checkpoints'/('cophylogeny_main25_recovery_case_'+str(s['scenario_index'])+'.complete.json')
                assert casecp.exists(), 'Stage complete but case checkpoint missing: '+str(casecp)
                checked.append(self.case(s))
            checked.sort(key=lambda x:x['summary']['scenario_index'])
            assert done['scenarios']==[x['summary'] for x in checked], 'Stage summaries differ from validated cases'
            self.finished=True;self.progress();return
        self.prepare_case_manifests();self.progress()
        # The designated primary case runs first. Independent sensitivities follow.
        scenarios=self.plan['scenarios'];assert scenarios[0]['id']==self.plan['primary_scenario']
        results=[self.case(scenarios[0])]
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as pool:
            futures={pool.submit(self.case,s):s for s in scenarios[1:]}
            for future in concurrent.futures.as_completed(futures):
                try:results.append(future.result())
                except BaseException as error:
                    with self.guard:self.failure=str(error)
                    self.progress()
                    for pending in futures:pending.cancel()
                    raise
        assert len(results)==16 and len(set(self.complete))==16
        results.sort(key=lambda x:x['summary']['scenario_index'])
        records={}
        for cp in results:
            assert_outputs(cp)
            for rec in cp['outputs']:records[rec['path']]=rec['sha256']
        summary={'status':'complete','stage':'cophylogeny_main25_recovery','finished_utc':now(),'config':self.config,'scenario_count':16,'scenarios':[x['summary'] for x in results],'outputs':[{'path':p,'sha256':h} for p,h in sorted(records.items())],'visual_review':'pending','bootstrap_topology_sensitivity':'not_run_in_this_stage','reconciliation_and_timewindows':'not_run_in_this_stage'}
        write_new(self.run/'checkpoints/cophylogeny_main25_recovery.complete.json',summary)
        self.finished=True;self.progress();print(json.dumps({'status':'complete','scenario_count':16,'primary_scenario':self.plan['primary_scenario']}),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-dir',type=Path,required=True);parser.add_argument('--threads',type=int,choices=range(1,5),default=4);parser.add_argument('--stage',choices=['all'],default='all');args=parser.parse_args()
    with (args.run_dir/'checkpoints/cophylogeny_main25_recovery_driver.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);driver=Driver(args)
        try:driver.execute()
        except BaseException as error:
            driver.failure=str(error);driver.progress();raise

if __name__=='__main__':main()
