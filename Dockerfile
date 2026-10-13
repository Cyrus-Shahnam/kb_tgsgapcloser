FROM kbase/sdkpython:3.8.10
LABEL maintainer="ac.shahnam, cshahnam"

USER root

# -----------------------------------------------------------------------------
# micromamba. RUN curl (not ADD, which ignored the proxy on one host), with a
# fallback mirror.
RUN ( curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xvj -C /usr/local bin/micromamba ) \
    || ( curl -fsSL https://github.com/mamba-org/micromamba-releases/releases/latest/download/micromamba-linux-64.tar.bz2 \
        | tar -xvj -C /usr/local bin/micromamba )

# Isolated tool env. It is deliberately NOT on the image PATH: it would shadow
# the server Python. The Runner calls tools by absolute path and puts this bin
# dir on the PATH of the subprocess only (tgsgapcloser calls minimap2 by name).
ENV MAMBA_ROOT_PREFIX=/opt/micromamba
COPY env-tgsgapcloser.yml /tmp/env-tgsgapcloser.yml
# micromamba exits 0 with "Empty environment created" if the file has no
# readable package list (empty, or lines joined by a bad paste): fail loudly.
RUN micromamba create -y -p /opt/conda/envs/tgsgapcloser -f /tmp/env-tgsgapcloser.yml \
 && micromamba clean -a -y \
 && test -x /opt/conda/envs/tgsgapcloser/bin/tgsgapcloser \
 || { echo "ENV EMPTY OR INCOMPLETE - env-tgsgapcloser.yml as received:"; cat -A /tmp/env-tgsgapcloser.yml; exit 1; }

# Build-time check with the runtime PATH, scoped to this RUN. Every helper the
# driver script needs is checked; the driver uses bash, util-linux getopt, awk.
RUN export PATH=/opt/conda/envs/tgsgapcloser/bin:$PATH \
 && tgsgapcloser -h | grep 'Version : 1.2.1' \
 && for t in minimap2 racon bash awk getopt; do \
        command -v "$t" >/dev/null || { echo "MISSING helper: $t"; exit 1; }; done \
 && { getopt -T >/dev/null; test $? -eq 4 || { echo "getopt is not util-linux getopt"; exit 1; }; } \
 && for b in tgsgapcandidate tgsgapcloser tgsseqgen tgsseqsplit; do \
        test -x /opt/conda/envs/tgsgapcloser/bin/tgsgapcloserbin/$b \
        || { echo "MISSING binary: tgsgapcloserbin/$b"; exit 1; }; done \
 && minimap2 --version && racon --version \
 && echo TGS_ENV_OK

# Server-side Python packages the minimal base does not bundle. A separate,
# later layer, so adding to it does not force a conda re-solve.
RUN pip install --no-cache-dir jsonrpcbase pytest pytest-cov coverage
RUN python -c "import sys; print(sys.executable)" \
 && python -c "import pytest, jsonrpcbase; print('server python OK')"
# -----------------------------------------------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
