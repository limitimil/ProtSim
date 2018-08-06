import os
import pandas as pd
import myModeller as mm
import datetime
from multiprocessing import Process


class Case(object):
    memberfiles = [
        'template_pair_note.csv.seg'
    ]
    memberrepos = [
        'fastas'
    ]
    outdirs = []
    def __init__(self, base='wildtype'):
        assert(os.path.isdir(base))
        self.base=os.path.abspath(base)
        assert(self._validity())
        self._prepare()
    def _validity(self):
        ret = True
        ret &= all(os.path.isfile(os.path.join(self.base,f) )\
            for f in self.memberfiles
        )
        ret &= all(os.path.isdir(os.path.join(self.base,f) )\
            for f in self.memberrepos
        )
        return ret
    def _prepare(self):
        for dd in self.outdirs:
            if os.path.isdir(os.path.join(self.base,dd)):
                continue
            else:
                os.mkdir(os.path.join(self.base,dd))
class Stage(Case):
    def errlog(self,msg ,filename='errlog', stdout=False):
        f = open(os.path.join(self.base, filename), 'a')
        f.write('#'*10 + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        + '#'*10 + '\n')
        f.write(msg + '\n')
        f.close()
        if stdout:
            print msg 
class Worker(Case):
    memberfiles = []
    memberrepos = []
    outdirs = []
    def run(self, cpu = 32):
        print 'this is worker, my id is {}'.format(self.base)
class parallel(object):
    def __init__(self, joblist=[]):
        self.jobs = map(lambda k: Worker(base = k), joblist)
    def progress(self, p):
        print 'Progress: %s/%s' %(p, len(self.jobs))
    def run(self, cpu=32):
        for i in xrange(0, len(self.jobs), cpu):
            self.progress(i)
            plist = []
            for s in self.jobs[i:i+cpu]:
                proc = Process(target = s.run)
                proc.start()
                plist.append(proc)
            for p in plist:
                p.join()
        print 'done'
