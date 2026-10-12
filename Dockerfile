FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n mbxdb-env -f /tmp/environment.yml \
    && micromamba clean -a -y

# conda-forge megahit 1.2.9 is built for AVX2 and dies (SIGILL) on older CPUs;
# install the upstream static generic x86-64 build (with build-time check) instead.
RUN mkdir -p /tmp/mh \
    && /opt/conda/envs/mbxdb-env/bin/python -c "import urllib.request; urllib.request.urlretrieve('https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz', '/tmp/mh/megahit.tar.gz')" \
    && tar xzf /tmp/mh/megahit.tar.gz -C /tmp/mh \
    && cp /tmp/mh/MEGAHIT-*-static/bin/* /opt/conda/envs/mbxdb-env/bin/ \
    && mkdir -p /opt/conda/envs/mbxdb-env/share \
    && cp -r /tmp/mh/MEGAHIT-*-static/share/megahit /opt/conda/envs/mbxdb-env/share/ \
    && for c in megahit_core megahit_core_popcnt megahit_core_no_hw_accel; do /opt/conda/envs/mbxdb-env/bin/$c dumpversion; done \
    && rm -rf /tmp/mh

ENV MAMBA_DEFAULT_ENV=mbxdb-env \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/mbxdb-env/bin:/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8
