#!/usr/bin/env python3
"""Run only the 100 frozen AHE/main25 UFBoot sensitivities with at most three R workers."""
import argparse, concurrent.futures, datetime, fcntl, hashlib, json, os, subprocess, threading, time
from pathlib import Path
PLAN_SHA = '9b492d9e87f75c74af94f38c0f61f19ce3ca10f7f63799f37370b7f033f1a03f'
DERIVED_FROM_WRAPPER_SHA = '5eafd0fd5714a61a744a836138d19a708b9829b2e55cc7889fd3c8fa795a91c7'
SELECTION_SHA = 'd04ab1a2fb60f6f34950de672f67a15ba304bef2370d4dea68b6eead0fabb943'
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
        self.run=args.run_dir.resolve();self.requested_threads=args.threads;self.threads=min(3,args.threads)
        self.plan_path=self.run/'01_PROVENANCE/cophylogeny_bootstrap100_100scenarios_20260908.json'
        self.r=self.run/'scripts/frozen/cophylogeny_exact_inputs_R26_df207f4e1cf052f2.R'
        assert sha(self.plan_path)==PLAN_SHA and sha(self.r)==R_SHA
        self.plan=json.loads(self.plan_path.read_text());assert len(self.plan['scenarios'])==100
        self.assert_plan_inputs()
        self.root=self.run/'04_COPHYLOGENY/bootstrap100_AHE_R26_20260908'
        self.root.mkdir(exist_ok=True)
        for name in ['inputs','analyses','logs']: (self.root/name).mkdir(exist_ok=True)
        self.config={'plan_sha256':PLAN_SHA,'R_sha256':R_SHA,'wrapper_sha256':sha(__file__),'requested_supervisor_threads':self.requested_threads,'max_concurrent_R_processes':self.threads,'reserved_CPU_for_other_stage':1,'derived_from_wrapper_sha256':DERIVED_FROM_WRAPPER_SHA,'selection_sha256':SELECTION_SHA,'R_compute_threads_each':1,'R_LIBS_USER':USER_R_LIBRARY,'R_LIBS_SITE':SITE_R_LIBRARY}
        configfile=self.root/'execution_config.json'
        if configfile.exists(): assert json.loads(configfile.read_text())==self.config
        else: write_new(configfile,self.config)
        self.guard=threading.Lock();self.active={};self.complete=[];self.failure=None;self.failed_cases={};self.finished=False
        self.env=r_environment()

    def assert_plan_inputs(self):
        p=self.plan
        assert p['schema']==1 and p['status']=='complete' and p['data_kind']=='real'
        assert p['nperm']==9999 and p['seed']==20260908 and p['correction']=='cailliez' and p['symmetric'] is False
        assert p['all_cases_are_bootstrap_sensitivities'] is True and 'primary_scenario' not in p
        assert p['source_selection']['sha256']==SELECTION_SHA and sha(p['source_selection']['path'])==SELECTION_SHA
        selection=json.loads(Path(p['source_selection']['path']).read_text())
        assert [s['bootstrap_source_tree_1based_index'] for s in p['scenarios']]==selection['one_based_tree_indices']
        assert len({s['id'] for s in p['scenarios']})==100
        assert len({s['scenario_index'] for s in p['scenarios']})==100
        coverage={}
        for rec in p['analysis_checkpoints']:
            assert sha(rec['path'])==rec['sha256']
            cp=json.loads(Path(rec['path']).read_text());assert cp['status']=='complete'
            records=list(cp.get('outputs',[]))+[{'path':x,'sha256':h} for x,h in cp.get('output_sha256',{}).items()]
            for r in records:
                assert Path(r['path']).is_absolute() and sha(r['path'])==r['sha256']
                coverage[r['path']]=r['sha256']
        evidence=p['association_evidence'];assert coverage[evidence['path']]==evidence['sha256']
        for s in p['scenarios']:
            assert s['role']=='bootstrap_sensitivity' and s['n_links']==25
            assert s['scenario_index']==1000+s['bootstrap_source_tree_1based_index']
            assert s['expected_R_seed']==p['seed']+s['scenario_index']
            assert s['id'].startswith('AHE_A495_dated__main_mfp_ufboot_')
            for r in s['inputs'].values():assert coverage[r['path']]==r['sha256']
        self.source_records=[{'path':str(self.plan_path),'sha256':PLAN_SHA},{'path':str(self.r),'sha256':R_SHA},p['source_selection']]+list(p['analysis_checkpoints'])
        self.source_records += [{'path':x,'sha256':h} for x,h in sorted(coverage.items())]

    def assert_case_summary(self,inner,scenario):
        summary=inner['scenarios'][0]
        assert summary['scenario_index']==scenario['scenario_index']
        assert summary['seed']==scenario['expected_R_seed']
        assert summary['n_links']==25 and summary['correction']=='cailliez' and summary['symmetric'] is False
        assert inner['script_sha256']==R_SHA
        for key,rec in scenario['inputs'].items():assert summary['input_hashes'][key+'_sha256']==rec['sha256']

    def progress(self):
        with self.guard:
            data={'updated_utc':now(),'driver_pid':os.getpid(),'status':'failed' if self.failure else ('complete' if self.finished else ('running' if self.active else 'between_cases')),'active':dict(self.active),'completed_scenarios':list(self.complete),'planned_scenarios':100,'failure':self.failure,'failed_cases':dict(self.failed_cases),'config':self.config}
            dest=self.run/'checkpoints/cophylogeny_bootstrap100_progress.json'
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
        sid=s['id']
        with self.guard:
            assert self.failure is None, 'No new cases after failure: '+str(self.failure)
        manifest=self.root/'inputs'/(sid+'.json');out=self.root/'analyses'/sid
        cpfile=self.run/'checkpoints'/('cophylogeny_bootstrap100_case_'+str(s['scenario_index'])+'.complete.json')
        cmd=[R_EXE,'--vanilla',str(self.r),'--inputs-json',str(manifest),'--out-dir',str(out)]
        if cpfile.exists():
            cp=json.loads(cpfile.read_text());assert cp['config']==self.config and cp['command']==cmd
            assert_outputs(cp)
            innerfile=out/'completed.json'
            assert any(Path(rec['path'])==innerfile for rec in cp['outputs']), 'Case completion is not bound to R completion'
            inner=assert_r_completion(innerfile,sid);self.assert_case_summary(inner,s)
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
        innerfile=out/'completed.json';inner=assert_r_completion(innerfile,sid);self.assert_case_summary(inner,s)
        assert sha(self.plan_path)==PLAN_SHA and sha(self.r)==R_SHA
        records=inner['outputs']+[{'path':str(innerfile),'sha256':sha(innerfile)},{'path':str(logpath),'sha256':sha(logpath)},{'path':str(manifest),'sha256':sha(manifest)}]
        cp={'status':'complete','scenario_id':sid,'finished_utc':now(),'config':self.config,'command':cmd,'R_exit_code':rc,'summary':inner['scenarios'][0],'outputs':records}
        write_new(cpfile,cp)
        with self.guard:self.active.pop(sid,None);self.complete.append(sid)
        self.progress();return cp

    def execute(self):
        existing=self.run/'checkpoints/cophylogeny_bootstrap100.complete.json'
        if existing.exists():
            done=json.loads(existing.read_text());assert done['config']==self.config;assert_outputs(done)
            assert done['scenario_count']==100 and len(done['scenarios'])==100
            assert {x['scenario_id'] for x in done['scenarios']}=={s['id'] for s in self.plan['scenarios']}
            checked=[]
            for s in self.plan['scenarios']:
                casecp=self.run/'checkpoints'/('cophylogeny_bootstrap100_case_'+str(s['scenario_index'])+'.complete.json')
                assert casecp.exists(), 'Stage complete but case checkpoint missing: '+str(casecp)
                checked.append(self.case(s))
            checked.sort(key=lambda x:x['summary']['scenario_index'])
            assert done['scenarios']==[x['summary'] for x in checked], 'Stage summaries differ from validated cases'
            self.finished=True;self.progress();return
        self.prepare_case_manifests();self.progress()
        # The first preselected replicate checks execution before bounded parallel dispatch.
        # Its order has no primary biological or statistical interpretation.
        scenarios=self.plan['scenarios'];assert scenarios[0]['id']==self.plan['first_execution_scenario']
        results=[self.case(scenarios[0])]
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as pool:
            remaining=iter(scenarios[1:]);pending={}
            for _ in range(self.threads):
                s=next(remaining,None)
                if s is not None:pending[pool.submit(self.case,s)]=s
            while pending:
                done,_=concurrent.futures.wait(pending,return_when=concurrent.futures.FIRST_COMPLETED)
                try:
                    for future in done:
                        pending.pop(future);results.append(future.result())
                except BaseException as error:
                    with self.guard:self.failure=str(error)
                    self.progress()
                    for future in pending:future.cancel()
                    raise
                # Submit only after all currently completed futures passed validation.
                for _ in done:
                    if self.failure:raise RuntimeError(self.failure)
                    s=next(remaining,None)
                    if s is not None:pending[pool.submit(self.case,s)]=s
        assert len(results)==100 and len(set(self.complete))==100
        results.sort(key=lambda x:x['summary']['scenario_index'])
        records={}
        for cp in results:
            assert_outputs(cp)
            for rec in cp['outputs']:records[rec['path']]=rec['sha256']
        self.assert_plan_inputs()
        summary={'status':'complete','stage':'cophylogeny_bootstrap100','finished_utc':now(),'config':self.config,'scenario_count':100,'scenarios':[x['summary'] for x in results],'source_hashes':self.source_records,'outputs':[{'path':p,'sha256':h} for p,h in sorted(records.items())],'visual_review':'pending','bootstrap_topology_sensitivity':'100 preselected main25 UFBoot trees; all cases are sensitivities','source_selection_sha256':SELECTION_SHA,'reference_primary_scenario':self.plan['reference_primary_scenario'],'first_execution_scenario':self.plan['first_execution_scenario'],'reconciliation_and_timewindows':'not_run_in_this_stage'}
        write_new(self.run/'checkpoints/cophylogeny_bootstrap100.complete.json',summary)
        self.finished=True;self.progress();print(json.dumps({'status':'complete','scenario_count':100,'all_cases_are_bootstrap_sensitivities':True}),flush=True)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run-dir',type=Path,required=True);parser.add_argument('--threads',type=int,choices=range(1,5),default=4);parser.add_argument('--stage',choices=['all'],default='all');args=parser.parse_args()
    with (args.run_dir/'checkpoints/cophylogeny_bootstrap100_driver.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB);driver=Driver(args)
        try:driver.execute()
        except BaseException as error:
            driver.failure=str(error);driver.progress();raise

if __name__=='__main__':main()
