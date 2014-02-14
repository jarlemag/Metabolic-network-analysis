import subprocess
import re

def system(*args, **kwargs):
    kwargs.setdefault('stdout', subprocess.PIPE)
    proc = subprocess.Popen(args, **kwargs)
    out, err = proc.communicate()
    return out

#system('pep8.py')

print '''
def system(*args, **kwargs):
    kwargs.setdefault('stdout', subprocess.PIPE)
    proc = subprocess.Popen(args, **kwargs)
    out, err = proc.communicate()
    return out"
'''

print "files = system('git', 'status', '--porcelain')"
files = system('git', 'status', '--porcelain')

print "modified = re.compile('^[AM]+\s+(?P<name>.*\.py)', re.MULTILINE)"
modified = re.compile('^[AM]+\s+(?P<name>.*\.py)', re.MULTILINE)

print 'files:'
print files

print 'modified_files = modified.findall(files)'
matches= modified.findall(files)

print 'matches:'
print matches

#THIS IS A TEST.
