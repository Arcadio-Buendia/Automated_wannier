import copy,importlib.util,json,os,subprocess,sys,tempfile,threading,unittest
from http.server import BaseHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import wannier_workflow as wf
import wannier_analysis as wa
from fake_program import xml

class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(prefix='Wannier workflow contract ')
        self.addCleanup(self.temp.cleanup);self.root=Path(self.temp.name)
        save=self.root/'original scratch'/'X.save';save.mkdir(parents=True)
        (save/'data-file-schema.xml').write_text(xml(6,[[0,0,0],[.5,0,0]]))
        (save/'charge-density.dat').write_bytes(b'original density')
        self.save=save;self.original_xml=(save/'data-file-schema.xml').read_bytes()
        (self.root/'scf.out').write_text('convergence has been achieved\nJOB DONE.\n')
        (self.root/'scf.in').write_text("&CONTROL\n prefix='X', outdir='original scratch', pseudo_dir='pseudo', calculation='scf'\n/\n&SYSTEM\n ibrav=1, celldm(1)=10, nat=1, ntyp=1, ecutwfc=30\n/\n&ELECTRONS\n/\nATOMIC_SPECIES\nX 1 X.upf\nATOMIC_POSITIONS crystal\nX .2 .3 .4\nK_POINTS automatic\n2 1 1 0 0 0\n")
        fake=Path(__file__).with_name('fake_program.py').resolve()
        self.config={'qe_input':'scf.in','dft_output':'scf.out','save_dir':'original scratch/X.save','mesh':[2,1,1],
                     'work_dir':'result','initial_nbnd':6,'target_relative_ev':[-1.,1.], 'max_attempts':12,
                     'commands':{k:[sys.executable,str(fake),v] for k,v in [('pw','pw'),('projwfc','projwfc'),('pw2wannier','pw2wannier'),('wannier','wannier')]}}
        self.config_path=self.root/'wannier_config.json'
        self.rules=json.loads((wf.HOME/'wannier_rules.json').read_text())

    def server(self,selector):
        requests=[]
        class Handler(BaseHTTPRequestHandler):
            def do_POST(handler):
                data=json.loads(handler.rfile.read(int(handler.headers['Content-Length'])))
                requests.append(data);choice=selector(data,len(requests))
                body=json.dumps({'answers':{'decision':{'choice':choice}}}).encode()
                handler.send_response(200);handler.send_header('Content-Type','application/json');handler.end_headers();handler.wfile.write(body)
            def log_message(self,*a):pass
        server=ThreadingHTTPServer(('127.0.0.1',0),Handler)
        thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
        self.addCleanup(server.server_close);self.addCleanup(server.shutdown)
        self.config['typesafe_url']=f'http://127.0.0.1:{server.server_port}/'
        return requests

    def run_pipeline(self, selector, flag=None, mode=None):
        requests=self.server(selector);self.config_path.write_text(json.dumps(self.config))
        if flag:
            (self.root/'result').mkdir();(self.root/'result'/flag).touch()
        env=dict(os.environ,TYPESAFE_API_KEY='contract-test-key',MPLCONFIGDIR=str(self.root/'mpl'))
        command=[sys.executable,str(wf.HOME/'wannier_workflow.py'),'--config',str(self.config_path)]
        if mode:command.append(mode)
        result=subprocess.run(command,env=env,capture_output=True,text=True,timeout=90)
        return result,requests

    def test_full_pipeline_local_repair_and_matrix_reuse(self):
        def choose(data,index):
            if index==1:return next(c['id'] for c in data['state']['candidates'] if c['id'].startswith('start_scdm_erfc'))
            ids=[c['id'] for c in data['state']['candidates']]
            self.assertFalse(any(i.startswith('switch_') for i in ids))
            return 'widen_outer' if 'widen_outer' in ids else 'reduce_trial_step'
        result,requests=self.run_pipeline(choose,'fail_first')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        history=json.loads((self.root/'result/history.json').read_text())
        self.assertEqual(len(history),2);self.assertEqual(history[0]['method'],history[1]['method'])
        events=(self.root/'result/events.txt').read_text()
        self.assertEqual(events.count('pw2wannier '),1,events)
        self.assertTrue((self.root/'result/bands_0006/fatbands.png').is_file())
        self.assertEqual((self.save/'data-file-schema.xml').read_bytes(),self.original_xml)
        self.assertEqual((self.save/'charge-density.dat').read_bytes(),b'original density')
        nscf=(self.root/'result/bands_0006/nscf.in').read_text()
        self.assertIn('ibrav = 0',nscf);self.assertNotIn('celldm',nscf);self.assertIn('X 0.00000000000000',nscf)
        self.assertEqual(requests[0]['model'],'jev-latest')

    def test_unsupported_atomic_switches_to_manual(self):
        def choose(data,index):
            prefix='start_atomic' if index==1 else 'switch_manual'
            return next(c['id'] for c in data['state']['candidates'] if c['id'].startswith(prefix))
        result,_=self.run_pipeline(choose,'unsupported_atomic')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        history=json.loads((self.root/'result/history.json').read_text())
        self.assertEqual([h['method'] for h in history],['atomic','manual'])
        self.assertIn('begin projections',(self.root/'result/attempt_002/wannier.win').read_text())

    def test_untrusted_api_choice_is_not_executed(self):
        result,_=self.run_pipeline(lambda *_:'run_arbitrary_shell_command')
        self.assertNotEqual(result.returncode,0)
        self.assertFalse((self.root/'result/attempt_001').exists())
        self.assertFalse((self.root/'result/.runner.lock').exists())

    def test_attempt_limit_is_failure_not_localized(self):
        self.config['max_attempts']=2
        def choose(data,index):
            return data['state']['candidates'][0]['id']
        result,_=self.run_pipeline(choose,'never_converge')
        self.assertEqual(result.returncode,2,result.stdout+result.stderr)
        self.assertFalse(json.loads((self.root/'result/status.json').read_text())['localized'])

    def test_analysis_only_makes_no_api_requests(self):
        result,requests=self.run_pipeline(lambda *_:'stop',mode='--analyze-only')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr);self.assertEqual(requests,[])

    def test_jev_can_request_initial_band_expansion(self):
        self.config['initial_nbnd']=3
        self.config['target_relative_ev']=[-1.,3.]
        def choose(data,index):
            candidates=data['state']['candidates']
            if index==1:return next(c['id'] for c in candidates if 'increase_bands' in c)
            return next(c['id'] for c in candidates if c['id'].startswith('start_scdm_erfc'))
        result,requests=self.run_pipeline(choose)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(len(requests),2)
        history=json.loads((self.root/'result/history.json').read_text())
        self.assertEqual(history[0]['action'],'increase_bands')
        requested=history[0]['requested_nbnd']
        self.assertGreaterEqual(requested,11)
        self.assertTrue((self.root/f'result/bands_{requested:04d}/fatbands.png').exists())

    def test_band_expansion_after_failure_rebuilds_matrices(self):
        self.config['target_relative_ev']=[-1.,10.]
        def choose(data,index):
            candidates=data['state']['candidates']
            if index==1:return next(c['id'] for c in candidates if c['id'].startswith('start_scdm_erfc'))
            if index==2:return next(c['id'] for c in candidates if 'increase_bands' in c)
            return candidates[0]['id']
        result,_=self.run_pipeline(choose,'fail_first')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        events=(self.root/'result/events.txt').read_text()
        self.assertEqual(events.count('pw '),2,events)
        self.assertEqual(events.count('pw2wannier '),2,events)

    def test_bloch_and_gaussian_end_to_end(self):
        for method in ('bloch_phases','scdm_gaussian'):
            with self.subTest(method=method):
                self.config['work_dir']='result_'+method
                def choose(data,index):
                    return next(c['id'] for c in data['state']['candidates'] if c['id'].startswith('start_'+method))
                result,_=self.run_pipeline(choose)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                status=json.loads((self.root/self.config['work_dir']/'status.json').read_text())
                self.assertEqual(status['plan']['method'],method)
                self.assertTrue(status['localized'])

    def test_lock_does_not_modify_an_active_run(self):
        self.config_path.write_text(json.dumps(self.config))
        work=self.root/'result';work.mkdir();lock=work/'.runner.lock';lock.mkdir();(lock/'pid').write_text('other-process')
        (work/'status.json').write_text('{"status":"active"}')
        with self.assertRaises(ValueError):wf.Workflow(self.config_path).run()
        self.assertEqual(json.loads((work/'status.json').read_text())['status'],'active')
        self.assertTrue(lock.exists())

    def test_neighbor_shell_repair_keeps_initialization(self):
        def choose(data,index):
            if index==1:return next(c['id'] for c in data['state']['candidates'] if c['id'].startswith('start_manual'))
            self.assertEqual([c['id'] for c in data['state']['candidates']],['expand_neighbor_search'])
            return 'expand_neighbor_search'
        result,_=self.run_pipeline(choose,'fail_neighbor')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        plan=json.loads((self.root/'result/attempt_002/plan.json').read_text())
        self.assertEqual(plan['method'],'manual');self.assertEqual(plan['search_shells'],200)

    def enable_existing_nscf(self):
        (self.root/'nscf.in').write_text((self.root/'scf.in').read_text().replace("calculation='scf'","calculation='nscf'"))
        (self.root/'nscf.out').write_text('JOB DONE.\n')
        (self.save/'wfc1.dat').write_bytes(b'saved NSCF wavefunctions')
        self.config.pop('initial_nbnd',None)
        self.config.pop('mesh',None)
        self.config['existing_nscf']={'folder':str(self.root)}

    def test_existing_nscf_skips_pw_and_preserves_source(self):
        self.enable_existing_nscf()
        self.config['commands']['pw']=['this_pw_must_never_be_called']
        def choose(data,index):
            return next(c['id'] for c in data['state']['candidates'] if c['id'].startswith('start_scdm_erfc'))
        result,_=self.run_pipeline(choose)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        events=(self.root/'result/events.txt').read_text()
        self.assertNotIn('pw ',events)
        self.assertIn('projwfc ',events)
        self.assertEqual((self.save/'wfc1.dat').read_bytes(),b'saved NSCF wavefunctions')
        self.assertEqual((self.save/'data-file-schema.xml').read_bytes(),self.original_xml)
        provenance=json.loads((self.root/'result/provenance.json').read_text())
        self.assertTrue(provenance['reused_nscf']);self.assertEqual(provenance['config']['mesh'],[2,1,1])

    def test_existing_nscf_band_expansion_runs_pw_only_when_needed(self):
        self.enable_existing_nscf();self.config['target_relative_ev']=[-1,10]
        def choose(data,index):
            candidates=data['state']['candidates']
            if index==1:return next(c['id'] for c in candidates if c['id'].startswith('start_scdm_erfc'))
            if index==2:return next(c['id'] for c in candidates if 'increase_bands' in c)
            return candidates[0]['id']
        result,_=self.run_pipeline(choose,'fail_first')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        events=(self.root/'result/events.txt').read_text()
        self.assertEqual(events.count('pw '),1,events)
        self.assertEqual(events.count('projwfc '),2,events)
        self.assertEqual((self.save/'data-file-schema.xml').read_bytes(),self.original_xml)

    def test_existing_nscf_missing_wavefunctions_is_rejected(self):
        self.enable_existing_nscf();(self.save/'wfc1.dat').unlink()
        self.config_path.write_text(json.dumps(self.config))
        with self.assertRaisesRegex(ValueError,'no saved wavefunctions'):
            wf.Workflow(self.config_path).run(analyze_only=True)
        self.assertFalse((self.root/'result/scratch').exists())

    def test_saved_mesh_accepts_reordering_and_shift(self):
        points=[[.75,0,0],[.25,0,0]]
        path=self.root/'shifted.xml';path.write_text(xml(6,points))
        data=wa.read_qe_xml(path)
        mesh,actual=wf.saved_mesh_points(path,data,[2,1,1])
        self.assertEqual(mesh,[2,1,1]);np.testing.assert_allclose(actual,points)
        with self.assertRaisesRegex(ValueError,'differs from requested'):
            wf.saved_mesh_points(path,data,[4,1,1])
        path.write_text(xml(6,[[0,0,0],[.3,0,0]]));data=wa.read_qe_xml(path)
        with self.assertRaisesRegex(ValueError,'not uniformly spaced'):
            wf.saved_mesh_points(path,data)

    def test_scdm_requires_explicit_gamma(self):
        data=wa.read_qe_xml(self.save/'data-file-schema.xml')
        data['kcart']=data['kcart']+.25
        # This check precedes orbital selection: a shifted mesh without Gamma is ineligible.
        self.assertIsNone(wf.prepare_plan('scdm_erfc','unused',{},data,self.rules['defaults']))
        self.assertIsNone(wf.prepare_plan('scdm_gaussian','unused',{},data,self.rules['defaults']))

    def test_shipped_reuse_example_selects_existing_nscf(self):
        commands=self.config['commands']
        self.enable_existing_nscf()
        self.config=json.loads((wf.HOME/'wannier_config.reuse.example.json').read_text())
        self.config['existing_nscf']['folder']=str(self.root)
        self.config['work_dir']='result'
        self.config['commands']=commands
        result,requests=self.run_pipeline(lambda *_:'stop',mode='--analyze-only')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(requests,[])
        self.assertNotIn('pw ',(self.root/'result/events.txt').read_text())
        self.assertTrue(json.loads((self.root/'result/provenance.json').read_text())['reused_nscf'])

    def test_configuration_missing_inputs_is_actionable(self):
        self.config.pop('qe_input');self.config.pop('dft_output')
        self.config_path.write_text(json.dumps(self.config))
        with self.assertRaises(ValueError) as caught:
            wf.Workflow(self.config_path)
        message=str(caught.exception)
        for hint in ('qe_input','dft_output','existing_nscf','--config',str(self.config_path)):
            self.assertIn(hint,message)
        self.assertFalse((self.root/'result').exists())

    def test_configuration_root_must_be_object(self):
        self.config_path.write_text('[]')
        with self.assertRaisesRegex(ValueError,'must contain a JSON object'):
            wf.Workflow(self.config_path)
        self.assertFalse((self.root/'result').exists())

    def test_existing_nscf_needs_no_fresh_calculation_fields(self):
        self.enable_existing_nscf()
        self.config.pop('qe_input');self.config.pop('dft_output');self.config.pop('save_dir')
        result,requests=self.run_pipeline(lambda *_:'stop',mode='--analyze-only')
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertEqual(requests,[])
        self.assertNotIn('pw ',(self.root/'result/events.txt').read_text())

    def test_window_capacity_checks_each_k(self):
        eig=np.array([[-2.,-.2,.2,4.],[-1.,.1,.3,5.]])
        plan={'num_wann':1,'num_bands':4,'exclude_bands':[],'outer':[-3,6],'frozen':[-.5,.5]}
        wf.repair_windows(plan,eig,0.)
        self.assertLessEqual(wa.window_counts(eig,plan)['maximum_frozen_states'],1)
        self.assertGreaterEqual(wa.window_counts(eig,plan)['minimum_outer_states'],1)

    def test_final_state_alone_is_not_convergence(self):
        folder=self.root/'attempt';folder.mkdir()
        (folder/'x.wout').write_text('Final State\nWF centre and spread 1 (0,0,0) 1.0\n')
        for suffix in ('.chk','_hr.dat'):(folder/('x'+suffix)).write_text('nonempty')
        r=wa.localization_result(folder,'x',{'num_wann':1,'num_bands':1})
        self.assertFalse(r['success'])

    def test_spinor_and_lsda_xml(self):
        path=self.root/'spinor.xml';path.write_text(xml(6,[[0,0,0]],True))
        d=wa.read_qe_xml(path);self.assertTrue(d['spinors']);self.assertAlmostEqual(d['eigenvalues'][0,0],-.8)
        text=xml(6,[[0,0,0]]).replace('<lsda>false</lsda>','<lsda>true</lsda><nbnd_up>3</nbnd_up><nbnd_dw>3</nbnd_dw>')
        path.write_text(text)
        with self.assertRaises(ValueError):wa.read_qe_xml(path)
        d=wa.read_qe_xml(path,'down');self.assertEqual(d['nbnd'],3);self.assertGreater(d['eigenvalues'][0,0],4)

    def test_bloch_phases_not_offered_for_entanglement(self):
        data=wa.read_qe_xml(self.save/'data-file-schema.xml')
        summary={'target_absolute_ev':[-1,1],'groups':[{'id':'s','manual_count':1,'atomic_count':1}],
                 'orbital_candidates':[{'id':'a','groups':['s']}], 'isolated_candidate':{'supported':False}}
        self.assertIsNone(wf.prepare_plan('bloch_phases','a',summary,data,self.rules['defaults']))
        with self.assertRaises(ValueError):
            wf.validate_plan({'method':'bloch_phases','num_wann':1,'num_bands':6,'exclude_bands':[],'outer':[-100,100],'frozen':None},data)

    def test_switch_not_offered_with_local_repairs(self):
        data=wa.read_qe_xml(self.save/'data-file-schema.xml')
        summary={'target_absolute_ev':[-1,1],'groups':[{'id':'s','manual_count':1,'atomic_count':1}],
                 'orbital_candidates':[{'id':'a','groups':['s']}], 'isolated_candidate':{'supported':False},
                 'top_band_below_target_top_at_some_k':False}
        plan=wf.prepare_plan('scdm_erfc','a',summary,data,self.rules['defaults'])
        candidates=wf.recovery_candidates(plan,{'failure_kind':'localization'},summary,data,self.rules,{wf.fingerprint(plan)},[])
        self.assertTrue(candidates);self.assertFalse(any(c['id'].startswith('switch') for c in candidates))
        records=[{'method':'scdm_erfc','action':c.get('family',c['id'])} for c in candidates]
        tried={wf.fingerprint(plan)}|{wf.fingerprint(c['plan']) for c in candidates if 'plan' in c}
        candidates=wf.recovery_candidates(plan,{'failure_kind':'localization'},summary,data,self.rules,tried,records)
        self.assertTrue(any(c['id'].startswith('switch') for c in candidates))

    def test_generated_modes_and_spinor_count(self):
        data=wa.read_qe_xml(self.save/'data-file-schema.xml');data['spinors']=True
        group={'id':'s','manual_count':2,'atomic_count':2,'indices':[1,2],'atom':1,'l':0}
        summary={'target_absolute_ev':[-1,1],'groups':[group], 'orbital_candidates':[{'id':'a','groups':['s']}],
                 'isolated_candidate':{'supported':False}}
        for method in ('manual','atomic','scdm_erfc','scdm_gaussian'):
            folder=self.root/method;folder.mkdir();plan=wf.prepare_plan(method,'a',summary,data,self.rules['defaults'])
            self.assertEqual(plan['num_wann'],2)
            wf.generate_inputs(folder,'x',plan,data,summary,{'mesh':[2,1,1]},self.root,'X',wf.mesh_points([2,1,1]))
            win=(folder/'x.win').read_text();pp=(folder/'x.pw2wan').read_text()
            if method=='manual':self.assertIn('(u,d)',win);self.assertNotIn('auto_projections',win)
            else:self.assertIn('auto_projections',win);self.assertNotIn('begin projections',win)
            if method=='atomic':self.assertIn('atom_proj = .true.',pp)
            if method.startswith('scdm'):self.assertIn('scdm_proj = .true.',pp)

if __name__=='__main__':unittest.main()
