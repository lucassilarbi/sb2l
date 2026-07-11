#!/usr/bin/env python3
# python3 launcher for the bundled waf 1.7.5.
# The stock ./waf embeds a bz2 blob with null bytes that python3 refuses to
# compile. Here we extract waflib from that blob and apply the few py3.11 fixes
# waf 1.7.5 needs, then hand over to waf's entry point.
import os, sys, shutil, tarfile

VERSION = "1.7.5"
REVISION = "f6a5ffc696be5d2549129d9fd105b5a3"
C1, C2 = '#-', '#%'

here = os.path.dirname(os.path.abspath(__file__))
wafsrc = os.path.join(here, 'waf')
wafdir = os.path.join(here, 'waf3-%s-%s' % (VERSION, REVISION))


def b(x):
    return x.encode()


def unpack():
    # Extract the tar.bz2 archive embedded between the #==> / #<== markers.
    f = open(wafsrc, 'rb')
    while True:
        line = f.readline()
        if not line:
            sys.exit('waf3: no archive found in %s' % wafsrc)
        if line == b('#==>\n'):
            txt = f.readline()
            if f.readline() != b('#<==\n'):
                sys.exit('waf3: corrupt archive')
            break
    txt = txt[1:-1].replace(b(C1), b('\n')).replace(b(C2), b('\r'))
    shutil.rmtree(wafdir, ignore_errors=True)
    for x in ['Tools', 'extras']:
        os.makedirs(os.path.join(wafdir, 'waflib', x))
    cwd = os.getcwd()
    os.chdir(wafdir)
    open('t.bz2', 'wb').write(txt)
    t = tarfile.open('t.bz2')
    for x in t:
        t.extract(x)
    t.close()
    os.unlink('t.bz2')
    os.chdir(cwd)


def patch_py3():
    # py3.11 fixes for waf 1.7.5, applied to the extracted waflib (idempotent).
    wl = os.path.join(wafdir, 'waflib')

    def edit(rel, olds_news):
        p = os.path.join(wl, rel)
        s = open(p).read()
        for old, new in olds_news:
            s = s.replace(old, new)
        open(p, 'w').write(s)

    # 'rU' file mode was removed in py3.11.
    edit('ConfigSet.py', [("m='rU'", "m='r'")])
    edit('Context.py', [("m='rU'", "m='r'"), ("node.read('rU')", "node.read('r')")])
    # PEP 479: a generator may no longer raise StopIteration.
    edit('Node.py', [("raise StopIteration", "return")])
    # Make the locally-defined node_class picklable under py3.
    edit('Context.py', [
        ('self.node_class.__name__="Nod3"',
         'self.node_class.__name__="Nod3"\n\t\tself.node_class.__qualname__="Nod3"'
         '\n\t\twaflib.Node.Nod3=self.node_class'),
    ])


if not os.path.isdir(os.path.join(wafdir, 'waflib')):
    unpack()
    patch_py3()

sys.path.insert(0, wafdir)
from waflib import Scripting
Scripting.waf_entry_point(os.getcwd(), VERSION, wafdir)
