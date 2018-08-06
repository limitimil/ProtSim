import os
import pandas as pd
import myModeller as mm
import datetime
from multiprocessing import Process
import core


class Stage1(core.Stage):
    outdirs = [
        'models'
    ]
    pdbs = 'pdbs'
    def __init__(self, base = 'wildtype', template_base = 'pdbs'):
        assert(os.path.isdir(template_base))
        super(Stage1, self).__init__(base= base)

        self.outdir = 'models'
        self.pdbs = os.path.abspath(template_base)
        self.fastas = 'fastas'
    def run(self):
        self.df =\
        pd.read_csv(os.path.join(self.base,'template_pair_note.csv.seg'))
        context = os.getcwd()
        os.chdir(os.path.join(self.base,self.outdir))
        result = self.df.apply(self.func, 1)
        os.chdir(context)
        self.result = pd.concat([self.df, result], 1)
        self.result.to_csv(os.path.join(self.base,'models.csv'), index=False)
    def func(self, row):
        def modelit(row):
            prepare = mm.prepare_modeling(Fasta = row['fasta'], PDB =
            row['pdb'], chain=row['chain'])
            mhandler = mm.modeller(pdbdir=self.pdbs)
            mhandler.pushPir(prepare['alnfile'])
            mhandler.pushTemplateFile(row['pdb'])
            if not mhandler.validate():
                raise Exception('parameter invalid for:'+str(row))
            ret = mhandler.run(prepare['target'])
            return ret
        def abspathp(path):
            return\
                os.path.join(self.pdbs, path)
        def abspathf(path):
            return\
                os.path.join(self.base, self.fastas,path)
        ###start here ###
        try:
            rowA = pd.Series({
                'fasta':abspathf(row['fastaA']),
                'pdb': abspathp(row['path']),
                'chain': row['codeA']
                })
            modelA = modelit(rowA)
        except Exception as e:
            self.errlog(str(e), stdout=True)
            modelA = None
        try:
            rowB = pd.Series({
                'fasta':abspathf(row['fastaB']),
                'pdb': abspathp(row['path']),
                'chain': row['codeB']
                })
            modelB = modelit(rowB)
        except Exception as e:
            self.errlog(str(e), stdout=True)
            modelB = None
        return pd.Series({
            'modelA': os.path.basename(modelA),
            'modelB': os.path.basename(modelB),
        })
class Stage1_parallel(object):
    template_base = None
    def __init__(self, joblist=[]):
        if not os.path.isdir(self.template_base):
            raise Exception('{} is not a valid path for PDB repository'.format(
                self.template_base))
        self.jobs = map(lambda k: Stage1(base = k,template_base = self.template_base), joblist)
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
