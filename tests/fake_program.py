"""Contract-test stand-ins, never used by the production runner."""
import json,re,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"workflow"))
from wannier_workflow import get_qe, clean_input
from wannier_analysis import HA_EV

def xml(nb, points, spinors=False):
    spectra=[]
    for i,k in enumerate(points):
        energies=[-.8+.1*i,.3+.1*i,2+.1*i]+[5+j+.1*i for j in range(max(0,nb-3))]
        spectra.append('<ks_energies><k_point>'+ ' '.join(map(str,k))+'</k_point><eigenvalues>'+ ' '.join(str(e/HA_EV) for e in energies[:nb])+'</eigenvalues><occupations>'+ ' '.join(['1']+['0']*(nb-1))+'</occupations></ks_energies>')
    return '<espresso><output><atomic_structure alat="10"><cell><a1>10 0 0</a1><a2>0 10 0</a2><a3>0 0 10</a3></cell><atomic_positions><atom name="X">0 0 0</atom></atomic_positions></atomic_structure><band_structure><noncolin>'+str(spinors).lower()+'</noncolin><lsda>false</lsda><fermi_energy>0</fermi_energy><nbnd>'+str(nb)+'</nbnd>'+''.join(spectra)+'</band_structure></output></espresso>'

if __name__=='__main__':
    tool=sys.argv[1]; cwd=Path.cwd()
    with (cwd.parent/'events.txt').open('a') as log: log.write(tool+' '+cwd.name+'\n')
    if tool=='pw':
        inp=clean_input(sys.stdin.read()); nb=int(get_qe(inp,'SYSTEM','nbnd'))
        tail=inp.split('K_POINTS crystal')[-1].strip().splitlines(); count=int(tail[0]);points=[list(map(float,l.split()[:3])) for l in tail[1:count+1]]
        save=Path(get_qe(inp,'CONTROL','outdir'))/(get_qe(inp,'CONTROL','prefix')+'.save')
        (save/'data-file-schema.xml').write_text(xml(nb,points))
        print('JOB DONE.')
    elif tool=='projwfc':
        inp=clean_input(sys.stdin.read());save=Path(get_qe(inp,'PROJWFC','outdir'))/(get_qe(inp,'PROJWFC','prefix')+'.save')
        import xml.etree.ElementTree as ET
        tree=ET.parse(save/'data-file-schema.xml')
        print('state # 1: atom 1 (X ), wfc 1 (l=0 m=1)')
        for idx in range(2,5): print(f'state # {idx}: atom 1 (X ), wfc 2 (l=1 m={idx-1})')
        for ks in tree.findall('.//ks_energies'):
            print('k = '+ks.findtext('k_point'))
            for band,e in enumerate(ks.findtext('eigenvalues').split(),1):
                print(f'==== e( {band}) = {float(e)*HA_EV:.5f} eV ====')
                if band==1: print('psi = 0.850*[# 1]+0.050*[# 2]+0.050*[# 3]+0.050*[# 4]')
                else: print('psi = 0.100*[# 1]+0.300*[# 2]+0.300*[# 3]+0.300*[# 4]')
                print('|psi|^2 = 1.000')
        print('JOB DONE.')
    elif tool=='pw2wannier':
        inp=clean_input(sys.stdin.read()); seed=get_qe(inp,'INPUTPP','seedname'); plan=json.loads((cwd/'plan.json').read_text())
        if (cwd.parent/'unsupported_atomic').exists() and plan['method']=='atomic':
            print('Error reading INPUTPP: unknown variable atom_proj'); sys.exit(1)
        for ext in ('.mmn','.amn'):
            size=1 if ext=='.mmn' else plan['num_wann']
            extra=f' {plan["scdm_mu"]:.6f} {plan["scdm_sigma"]:.6f}' if ext=='.amn' and plan['method'].startswith('scdm') else ''
            (cwd/(seed+ext)).write_text(f'contract fixture\n{plan["num_bands"]} 2 {size}{extra}\n')
        (cwd/(seed+'.eig')).write_text('1 1 0.0\n')
        print('JOB DONE.')
    elif tool=='wannier':
        seed=sys.argv[-1]
        if '-pp' in sys.argv:
            if (cwd.parent/'fail_neighbor').exists() and cwd.name=='attempt_001':
                (cwd/(seed+'.werr')).write_text('Unable to find b-vectors: not enough shells')
                sys.exit(1)
            (cwd/(seed+'.nnkp')).write_text('contract fixture\n');sys.exit(0)
        plan=json.loads((cwd/'plan.json').read_text()); nw=plan['num_wann']
        if (cwd.parent/'fail_first').exists() and cwd.name=='attempt_001':
            (cwd/(seed+'.wout')).write_text('Disentanglement convergence criteria satisfied\nFinal State\n')
            sys.exit(0)
        if (cwd.parent/'never_converge').exists():
            (cwd/(seed+'.wout')).write_text('Final State\n');sys.exit(0)
        content='Disentanglement convergence criteria satisfied\nWannierisation convergence criteria satisfied\nFinal State\n'
        content+=''.join(f'WF centre and spread {i} ( 0.0, 0.0, 0.0 ) 1.234\n' for i in range(1,nw+1))
        (cwd/(seed+'.wout')).write_text(content);(cwd/(seed+'.chk')).write_text('checkpoint')
        hr=f'contract fixture\n{nw}\n1\n1\n'
        for i in range(1,nw+1):
            for j in range(1,nw+1): hr+=f'0 0 0 {i} {j} {float(i==j)*.3} 0.0\n'
        (cwd/(seed+'_hr.dat')).write_text(hr)
