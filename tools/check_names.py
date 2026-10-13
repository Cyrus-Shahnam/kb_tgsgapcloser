#!/usr/bin/env python3
"""Name-chain check (lessons §2, §17). Run from the module root after every make."""
import glob, json, os, re, sys, yaml
errors = []
k = yaml.safe_load(open('kbase.yml'))['module-name'].strip()
spec = open(k + '.spec').read()
impl = open('lib/%s/%sImpl.py' % (k, k)).read()
names = {'kbase.yml': k,
         'KIDL': re.search(r'^\s*module\s+(\w+)', spec, re.M).group(1),
         'Impl class': re.search(r'^class\s+(\w+)', impl, re.M).group(1)}
for d in sorted(glob.glob('ui/narrative/methods/*/')):
    sj = json.load(open(os.path.join(d, 'spec.json')))
    smap = sj['behavior']['service-mapping']
    names['service-mapping ' + d] = smap['name']
    m = smap['method']
    if not re.search(r'^\s*funcdef\s+%s\b' % m, spec, re.M): errors.append('no funcdef ' + m)
    if '#BEGIN ' + m not in impl or '#END ' + m not in impl: errors.append('markers ' + m)
    sjp = {p['id'] for p in sj['parameters']}
    dy = set(yaml.safe_load(open(os.path.join(d, 'display.yaml')))['parameters'])
    im = {x['input_parameter'] for x in smap['input_mapping'] if 'input_parameter' in x}
    if sjp != dy: errors.append('display.yaml mismatch %s' % (sjp ^ dy))
    if not sjp <= im: errors.append('unwired params %s' % (sjp - im))
if len(set(names.values())) != 1: errors.append('names %s' % names)
print('\n'.join(errors) or 'all consistent: ' + k); sys.exit(1 if errors else 0)
