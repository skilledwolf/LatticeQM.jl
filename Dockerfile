FROM jupyter/scipy-notebook

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

ARG julia_channel

# install Julia packages in /opt/julia instead of ${HOME}
ENV JULIA_DEPOT_PATH=/opt/julia \
    JULIA_PKGDIR=/opt/julia \
    JULIA_DEFAULT_CHANNEL="${julia_channel:-release}"

RUN apt-get update && \
    apt-get install -y curl && \
    rm -rf /var/lib/apt/lists/*

# Show Julia where conda libraries are \
RUN mkdir /etc/julia && \
    echo "push!(Libdl.DL_LOAD_PATH, \"${CONDA_DIR}/lib\")" >> /etc/julia/juliarc.jl && \
    # Create JULIA_PKGDIR \
    mkdir -p "${JULIA_PKGDIR}" && \
    chown "${NB_USER}" "${JULIA_PKGDIR}" && \
    fix-permissions "${JULIA_PKGDIR}"


# Copy the julia package into the container
RUN mkdir /LatticeQM
COPY Project.toml /LatticeQM/Project.toml
COPY src /LatticeQM/src

RUN fix-permissions /LatticeQM && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"

RUN mkdir -p /opt/julia && \
    chown "${NB_USER}" /opt/julia && \
    fix-permissions /opt/julia

USER $NB_UID

RUN curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel "${JULIA_DEFAULT_CHANNEL}"

ENV PATH="/home/${NB_USER}/.juliaup/bin:${PATH}"

RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'import Pkg; Pkg.develop(path="/LatticeQM"); Pkg.add("Plots"); Pkg.add("ProgressMeter"); Pkg.add("IJulia"); Pkg.build(); Pkg.precompile();'

USER root

RUN fix-permissions /home/${NB_USER}
