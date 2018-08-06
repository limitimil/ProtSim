import core
import TemplateBasedDock as TBD

import os,re,json
import tempfile
import pandas as pd
from shutil import copy
from multiprocessing import Process
class Stage2(core.Stage, TBD.TemplateBasedDock):
    outdirs = [
        'alignments'
    ]
    memberfiles = [
        'models.csv'
    ]
    memberrepos = [
        'models'
    ]
    def __init__(self, base = 'wildtype', template_base= 'pdbs'):
        assert(template_base)
        core.Stage.__init__(self, base = base)
        TBD.TemplateBasedDock.__init__(self)

        self.outdir = 'alignments'
        self.models = 'models'

        self.PDBrepo = os.path.join(self.base, self.models)
        self.TEMPLATErepo = os.path.abspath(template_base)
        
        self.INTERPRED_TMP = tempfile.mkdtemp()
        self.INTERPRED_OUT = tempfile.mkdtemp()
    def _env_set(self):
        TBD.TemplateBasedDock._env_set(self)
        os.environ['TMP'] = self.INTERPRED_TMP
        os.environ['OUT'] = self.INTERPRED_OUT
    def grabPDBpath(self, target):
        files = os.listdir(self.PDBrepo)
        files = filter(
            lambda f: re.match('{}.pdb'.format(target),f),
            files
        )
        assert(len(files) == 1)
        return map(lambda f: os.path.join(self.PDBrepo, f),files)
    def run(self):
        def func2(row):
            result = row[2]
            pn = os.path.join(self.base, 'alignments')
            if len(result['models_files']) != 1:
                self.errlog('zero or multiple models generated for\n {}'.format(
                    row.to_string()))
            model = os.path.join(pn, os.path.basename(result['models_files'][0]))
            copy(result['models_files'][0], model)
            predict = os.path.join(pn, os.path.basename(result['predictions']))
            copy(result['predictions'], predict)
            od = os.path.join(self.base, self.outdir)
            return pd.Series({'complex':os.path.basename(result['models_files'][0]) ,
            'IPscore': os.path.basename(result['predictions'])})
        self.note = pd.read_csv(os.path.join(self.base,'models.csv'))
        res = self.note.apply(self.func, 1)
        #save alignment result into json file
        gn = os.path.join(self.base,'alignments.json')
        json.dump(res.values.tolist(), open(gn,'w'),indent=4)
        #arrange the repository of script/pdb file/prediction score
        self.complexes = res.apply(func2, 1)
        #save complexes configuration
        cn = os.path.join(self.base, 'complexes.csv')
        ret = pd.concat([self.note, self.complexes], 1)
        ret.to_csv(cn, index=False)
        return (gn,cn)
    def func(self, row):
        def preprocess(s):
            return os.path.splitext(os.path.basename(s))[0]
        ret = []
        pdb1 = self.grabPDBpath(preprocess(row['modelA']))[0] 
        pdb1 = os.path.abspath(pdb1)
        pdb2 = self.grabPDBpath(preprocess(row['modelB']))[0] 
        pdb2 = os.path.abspath(pdb2)
        template_path = os.path.join(self.TEMPLATErepo, row['path'])
        pdbt =[
            template_path,row['codeA'],row['codeB']] #
        try:
            print '\n[TBD] working on {} + {} = {}'.format(
                os.path.basename(pdb1), os.path.basename(pdb2),
                os.path.basename(pdbt[0]))
            ret = self.make_model(pdb1, pdb2, pdbt)
        except Exception as e:
            self._env_restore()
            self.errlog(str(e), stdout=True)
            ret = []
        return ret

class Stage2_parallel(object):
    template_base = None 
    def __init__(self, joblist=[]):
        if not os.path.isdir(self.template_base):
            raise Exception('{} is not a valid path for PDB repository'.format(
                self.template_base))
        self.jobs = map(lambda k: Stage2(base = k,template_base = self.template_base), joblist)
    def progress(self, p):
        print 'Progress: %s/%s' %(p, len(self.jobs))
    def run(self, cpu=32):
        for i in xrange(0, len(self.jobs), cpu):
            self.progress(i)
            plist = []
            for s in self.jobs[i:i+cpu]:
                print s.base
                proc = Process(target = s.run)
                proc.start()
                plist.append(proc)
            for p in plist:
                p.join()
        print 'done'
