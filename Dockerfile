FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml

RUN micromamba create -y -n mbxdb-env -f /tmp/environment.yml \
    && micromamba clean -a -y

ENV MAMBA_DEFAULT_ENV=mbxdb-env \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PATH=/opt/conda/envs/mbxdb-env/bin:/opt/conda/bin:$PATH \
    LC_ALL=C.UTF-8
