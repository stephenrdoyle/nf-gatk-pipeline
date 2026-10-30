#!/usr/bin/env python3
"""
Static checks on the shell and Python code the pipeline runs.

Nextflow renders each process's script: block as a Groovy string before bash
sees it. Groovy interpolates ${...} and $name, and turns escapes such as \\n
and \\t into real characters, so a script can look right in the .nf file and
still be broken when it runs (an embedded Python '\\n' becoming a literal
newline, for example). This renders every script: and stub: block the way
Groovy would, with interpolated values replaced by a placeholder, then:
  - runs `bash -n` on the rendered script
  - compiles every embedded Python heredoc (<< 'PYEOF' ... PYEOF)
It also compiles bin/*.py and runs `bash -n` on run_pipeline.sh.

Usage: python3 tests/check_scripts.py [repo_dir]      (exit 1 on any failure)
"""
import ast
import glob
import os
import re
import subprocess
import sys
import tempfile

ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PLACEHOLDER = 'GROOVY_VALUE'
GROOVY_ESCAPES = {'\\': '\\', '$': '$', 't': '\t', 'n': '\n', 'r': '\r',
                  'b': '\b', 'f': '\f', '"': '"', "'": "'"}


def interpolations_to_placeholder(text):
    """Replace Groovy interpolations (unescaped ${...} and $name) with a placeholder."""
    out, i = [], 0
    while i < len(text):
        c = text[i]
        if c == '\\' and i + 1 < len(text):          # escaped char: keep for unescape step
            out.append(text[i:i + 2]); i += 2; continue
        if c == '$' and text[i + 1:i + 2] == '{':     # ${ ... } with nested braces
            depth, j = 0, i + 1
            while j < len(text):
                if text[j] == '{': depth += 1
                elif text[j] == '}':
                    depth -= 1
                    if depth == 0: break
                j += 1
            out.append(PLACEHOLDER); i = j + 1; continue
        m = re.match(r'\$[A-Za-z_][\w.]*', text[i:])    # $name / $name.property
        if c == '$' and m:
            out.append(PLACEHOLDER); i += len(m.group(0)); continue
        out.append(c); i += 1
    return ''.join(out)


def groovy_unescape(text):
    """Apply Groovy string escapes; backslash-newline is a line continuation."""
    out, i = [], 0
    while i < len(text):
        if text[i] == '\\' and i + 1 < len(text):
            nxt = text[i + 1]
            if nxt == '\n':
                i += 2; continue
            out.append(GROOVY_ESCAPES.get(nxt, '\\' + nxt)); i += 2; continue
        out.append(text[i]); i += 1
    return ''.join(out)


def script_blocks(src):
    """(kind, start_line, body) for each script:/stub: block's triple-quoted string."""
    # The script is the LAST triple-quoted string in the section: earlier ones
    # (e.g. `def x = cond ? """...""" : """..."""`) are Groovy values. Triple
    # quotes are paired in order; the section ends at the next stub:/script:
    # or at the process's closing brace in column 0.
    section_end = re.compile(r'^(?:[ \t]*(?:stub|script):|\})', re.M)
    for m in re.finditer(r'^\s*(script|stub):\s*$', src, flags=re.M):
        pos, last = m.end(), None
        while True:
            opening = src.find('"""', pos)
            if opening == -1:
                break
            if last and section_end.search(src[pos:opening]):
                break
            closing = src.find('"""', opening + 3)
            if closing == -1:
                break
            last, pos = (opening + 3, closing), closing + 3
        if last:
            start, end = last
            yield m.group(1), src.count('\n', 0, start) + 1, src[start:end]


def embedded_python(script):
    """(label, code) for Python run from a rendered bash script, as
    `python3 << 'TAG' ... TAG` heredocs or `python3 -c "..."` / `-c '...'`."""
    for m in re.finditer(r"python3?\s+(?:-\s+)?<<-?\s*(['\"]?)(\w+)\1\n(.*?)\n\s*\2\n", script, flags=re.S):
        yield f'heredoc {m.group(2)}', m.group(3)
    for m in re.finditer(r"python3?\s+-c\s+(['\"])", script):
        quote, i, code = m.group(1), m.end(), []
        while i < len(script) and script[i] != quote:
            # inside "...", bash removes the backslash before \ " $ `
            if quote == '"' and script[i] == '\\' and script[i + 1:i + 2] in ('\\', '"', '$', '`'):
                code.append(script[i + 1]); i += 2; continue
            code.append(script[i]); i += 1
        yield 'python3 -c', ''.join(code)


failures = []

for nf in sorted(glob.glob(os.path.join(ROOT, 'modules', '*.nf'))):
    name = os.path.relpath(nf, ROOT)
    src = open(nf).read()
    for kind, line, body in script_blocks(src):
        rendered = groovy_unescape(interpolations_to_placeholder(body))
        where = f'{name}:{line} ({kind})'

        with tempfile.NamedTemporaryFile('w', suffix='.sh', delete=False) as tmp:
            tmp.write(rendered)
        res = subprocess.run(['bash', '-n', tmp.name], capture_output=True, text=True)
        os.unlink(tmp.name)
        if res.returncode != 0:
            failures.append(f'{where}: bash -n: {res.stderr.strip().replace(tmp.name, "script")}')

        for label, py in embedded_python(rendered):
            try:
                ast.parse(py)
            except SyntaxError as e:
                failures.append(f'{where}: embedded Python ({label}), line {e.lineno}: {e.msg}: {(e.text or "").strip()!r}')

# Containers are assigned only in nextflow.config (withName ... container =).
# Singularity is enabled globally, so a container/conda directive left in a
# module would silently move that process off its farm module - changing its
# cache key and, for VALIDATE_BAM, re-running every HaplotypeCaller job.
for nf in sorted(glob.glob(os.path.join(ROOT, 'modules', '*.nf'))):
    for n, line in enumerate(open(nf), 1):
        if re.match(r'\s*(container|conda)\s', line):
            failures.append(f'{os.path.relpath(nf, ROOT)}:{n}: process declares its own '
                            f'{line.split()[0]}; assign containers in nextflow.config instead')

for py in sorted(glob.glob(os.path.join(ROOT, 'bin', '*.py'))):
    try:
        compile(open(py).read(), py, 'exec')
    except SyntaxError as e:
        failures.append(f'{os.path.relpath(py, ROOT)}:{e.lineno}: {e.msg}')

run_sh = os.path.join(ROOT, 'run_pipeline.sh')
if os.path.exists(run_sh):
    res = subprocess.run(['bash', '-n', run_sh], capture_output=True, text=True)
    if res.returncode != 0:
        failures.append(f'run_pipeline.sh: bash -n: {res.stderr.strip()}')

n_modules = len(glob.glob(os.path.join(ROOT, 'modules', '*.nf')))
if failures:
    print(f'FAILED: {len(failures)} problem(s)')
    for f in failures:
        print('  - ' + f)
    sys.exit(1)
print(f'OK: {n_modules} modules, bin/*.py and run_pipeline.sh pass')
