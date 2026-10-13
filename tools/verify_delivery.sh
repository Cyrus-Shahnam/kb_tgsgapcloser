#!/bin/bash
# Verify the delivered kb_tgsgapcloser files on the dev server (lessons §0, §15).
# Run from the module root (the folder containing kbase.yml):  bash tools/verify_delivery.sh
# Checks each file exists, is non-empty, has the expected first line and about
# the expected line count (+0..3 blank lines from pasting is fine), has no
# terminal heredoc wrapper, and parses. Then checks init leftovers and Makefile.
cd "$(dirname "$0")/.." || exit 1
[ -f kbase.yml ] || { echo "run from the module root (no kbase.yml here)"; exit 1; }
bad=0; ok=0
note() { echo "PROBLEM  $1"; bad=$((bad+1)); }
while IFS=$'\t' read -r path want first; do
  [ -z "$path" ] && continue
  if [ ! -f "$path" ]; then note "$path: missing"; continue; fi
  if [ ! -s "$path" ]; then note "$path: EMPTY (failed paste?)"; continue; fi
  got=$(wc -l < "$path"); line1=$(head -1 "$path" | tr -d '\r' | sed 's/[[:space:]]*$//')
  want1=$(printf '%s' "$first" | sed 's/[[:space:]]*$//')
  if [ "$line1" != "$want1" ]; then
    case "$line1" in cat*'<<'*) note "$path: starts with a terminal heredoc line -> sed -i '1d;/^EOF\$/,\$d' $path";;
      *) note "$path: first line is [$line1], expected [$want1] (swapped or wrong file?)";; esac
    continue
  fi
  if tail -1 "$path" | grep -qx 'EOF'; then note "$path: ends with EOF (heredoc wrapper pasted into file)"; continue; fi
  if [ "$got" -lt "$want" ] || [ "$got" -gt $((want+3)) ]; then
    note "$path: $got lines, expected $want (truncated paste or extra content?)"; continue; fi
  case "$path" in
    *.py) python3 -c "import ast,sys;ast.parse(open(sys.argv[1]).read())" "$path" 2>/dev/null || { note "$path: Python does not parse"; continue; } ;;
    *.json) python3 -c "import json,sys;json.load(open(sys.argv[1]))" "$path" 2>/dev/null || { note "$path: JSON does not parse"; continue; } ;;
    *.yml|*.yaml) if python3 -c "import yaml" 2>/dev/null; then
        python3 -c "import yaml,sys;yaml.safe_load(open(sys.argv[1]))" "$path" 2>/dev/null || { note "$path: YAML does not parse"; continue; }; fi ;;
    *.sh) bash -n "$path" 2>/dev/null || { note "$path: bash syntax error"; continue; } ;;
  esac
  ok=$((ok+1))
done <<'DATA'
.dockerignore	4	test_local
Dockerfile	56	FROM kbase/sdkpython:3.8.10
README.md	136	# kb_tgsgapcloser
env-tgsgapcloser.yml	8	name: tgsgapcloser
kb_tgsgapcloser.spec	57	/*
kbase.yml	17	module-name:
lib/kb_tgsgapcloser/kb_tgsgapcloserImpl.py	85	# -*- coding: utf-8 -*-
lib/kb_tgsgapcloser/tgsgapcloser_utils.py	678	# -*- coding: utf-8 -*-
scripts/entrypoint.sh	30	#!/bin/bash
test/data/make_gap_fixture.py	117	#!/usr/bin/env python3
test/kb_tgsgapcloser_server_test.py	322	# -*- coding: utf-8 -*-
tools/check_names.py	24	#!/usr/bin/env python3
tools/dry_run_local.py	100	#!/usr/bin/env python3
tools/fake_tgsgapcloser.py	86	#!/usr/bin/env python3
tools/plant_gaps.py	99	#!/usr/bin/env python3
tools/score_gapfill.py	137	#!/usr/bin/env python3
ui/narrative/methods/run_tgsgapcloser/display.yaml	148	#
ui/narrative/methods/run_tgsgapcloser/spec.json	372	{
DATA
[ -d ui/narrative/methods/run_kb_tgsgapcloser ] && note "init example app ui/narrative/methods/run_kb_tgsgapcloser still present: rm -rf it"
[ "$(grep -c '^ *funcdef' kb_tgsgapcloser.spec)" = "1" ] || note "kb_tgsgapcloser.spec should have exactly 1 funcdef (init's spec not replaced?)"
grep -qx '  - tgsgapcloser=1.2.1' env-tgsgapcloser.yml 2>/dev/null || note "env-tgsgapcloser.yml has no '  - tgsgapcloser=1.2.1' line (empty or lines joined? check: cat -A env-tgsgapcloser.yml)"
grep -nE 'exit1|echo"|getopt-T|&&for|for bin ' Dockerfile >/dev/null && note "Dockerfile has joined words (lost spaces from pasting): grep -nE 'exit1|echo\"|getopt-T|&&for' Dockerfile"
grep -n -- '--java' Makefile >/dev/null 2>&1 && note "Makefile still has --java flags (see guide step 3)"
for c in AssemblyUtilClient ReadsUtilsClient WorkspaceClient KBaseReportClient; do
  [ -f lib/installed_clients/$c.py ] || note "lib/installed_clients/$c.py missing (kb-sdk install ...)"; done
find . -type f -empty -not -path './test_local/*' -not -path './.git/*' -not -name '__init__.py' | sed 's/^/PROBLEM  empty file: /'
echo "----"; echo "$ok files OK, $bad problem(s)"; [ "$bad" -eq 0 ]
