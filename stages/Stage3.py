import core
import collectChains

import pandas as pd
import os
from tqdm import tqdm
class Stage3(core.Stage):
    memberfiles = ['complexes.csv']
    memberrepos = ['alignments']
    outdirs = ['alignments_fiberdock']
    def __init__(self, base='wildtype'):
        super(Stage3,self).__init__(base = base)
        
        self.indir = 'alignments'
        self.outdir = 'alignments_fiberdock'
        assert(self.check_members())
    def check_members(self):
        tab = pd.read_csv(os.path.join(self.base,'complexes.csv'))
        print self.base 
        cc = collectChains.FiberDockRender(
            inpath = os.path.join(self.base, self.indir), 
            outpath = os.path.join(self.base, self.outdir))
        return set(tab.get('complex',pd.Series()).tolist()) == set(cc.files)
    def run(self):
        cc = collectChains.FiberDockRender(
            inpath = os.path.join(self.base, self.indir), 
            outpath = os.path.join(self.base, self.outdir))
        cc.collect()
class _checker(core.Stage):
    memberfiles = ['template_pair_note.csv.seg']
    memberrepos = []
    outdirs = []
    def scriptFix(self, tab):
        def func(row):
            if isinstance(row['codeA'], int) or isinstance(row['codeB'], int):
                return False
            return not (row['codeA'].isdigit() or row['codeB'].isdigit()) 
        return tab.loc[tab.apply(func, 1), :]
    def is_empty(self):
        tab = pd.read_csv(os.path.join(self.base, 'template_pair_note.csv.seg'))
        if tab.shape[0] == 0:
            return True
        tab = self.scriptFix(tab)
        return tab.shape[0] == 0
class Stage3_Sequential(object):
    def __init__(self, joblist):
        tmp = {k: _checker(k).is_empty() for k in joblist}
        self.empty_case = map(lambda k: k[0], filter(lambda k: k[1],
            tmp.items()))
        joblist = map(lambda k: k[0], filter(lambda k: not k[1],
            tmp.items()))
        self.jobs = map(lambda k: Stage3(base = k), joblist)
    def run(self):
        self.errorlist = []
        for jj in tqdm(self.jobs, ascii=True):
            try:
                jj.run()
            except Exception as e:
                self.errorlist.append(jj.base)
                continue
class Stage3_parallel(core.parallel):
    def __init__(self, joblist=[]):
        tmp = {k: _checker(k).is_empty() for k in joblist}
        self.empty_case = map(lambda k: k[0], filter(lambda k: k[1],
            tmp.items()))
        joblist = map(lambda k: k[0], filter(lambda k: not k[1],
            tmp.items()))
        self.jobs = map(lambda k: Stage3(base = k), joblist)

