# syntax=docker/dockerfile:1.7
FROM quay.io/jupyter/scipy-notebook:2024-12-23

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ARG julia_channel=release

USER root

# Install Julia packages in a shared depot so they're not stuck in $HOME
ENV JULIA_DEPOT_PATH=/opt/julia \
    JULIA_PKGDIR=/opt/julia \
    JULIA_DEFAULT_CHANNEL=${julia_channel}

RUN apt-get update && \
    apt-get install -y --no-install-recommends curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# Tell Julia where conda libraries live so PyCall etc. can find them. Julia
# reads $(DEPOT_PATH[1])/config/startup.jl at startup; with
# JULIA_DEPOT_PATH=/opt/julia that is /opt/julia/config/startup.jl.
RUN mkdir -p "${JULIA_PKGDIR}/config" && \
    printf 'import Libdl\npush!(Libdl.DL_LOAD_PATH, "%s/lib")\n' "${CONDA_DIR}" \
        > "${JULIA_PKGDIR}/config/startup.jl" && \
    chown -R "${NB_USER}" "${JULIA_PKGDIR}" && \
    fix-permissions "${JULIA_PKGDIR}"

# Stage the package itself; src/ is intentionally the only payload from the repo
RUN mkdir /LatticeQM
COPY Project.toml /LatticeQM/Project.toml
COPY src /LatticeQM/src

RUN fix-permissions /LatticeQM && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

USER ${NB_UID}

# Install juliaup non-interactively, then resolve + precompile the LatticeQM env
RUN curl -fsSL https://install.julialang.org | \
        sh -s -- --yes --default-channel "${JULIA_DEFAULT_CHANNEL}"

ENV PATH="/home/${NB_USER}/.juliaup/bin:${PATH}"

RUN julia -e ' \
        import Pkg; \
        Pkg.develop(path="/LatticeQM"); \
        Pkg.add(["Plots", "ProgressMeter", "IJulia"]); \
        Pkg.precompile(); \
    '

USER root

RUN fix-permissions /home/${NB_USER} && \
    fix-permissions "${JULIA_PKGDIR}"

USER ${NB_UID}
