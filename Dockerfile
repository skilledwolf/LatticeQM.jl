FROM "jupyter/scipy-notebook"

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

USER root

# install Julia packages in /opt/julia instead of ${HOME}
ENV JULIA_DEPOT_PATH=/opt/julia \
    JULIA_PKGDIR=/opt/julia \
    JULIA_VERSION="${julia_version:-1.8.2}"

RUN mkdir /opt/julia-${JULIA_VERSION} && \
    cd /tmp && \
    wget -q https://julialang-s3.julialang.org/bin/linux/x64/`echo ${JULIA_VERSION} | cut -d. -f 1,2`/julia-${JULIA_VERSION}-linux-x86_64.tar.gz && \
    tar xzf julia-${JULIA_VERSION}-linux-x86_64.tar.gz -C /opt/julia-${JULIA_VERSION} --strip-components=1 && \
    rm /tmp/julia-${JULIA_VERSION}-linux-x86_64.tar.gz

RUN ln -fs /opt/julia-*/bin/julia /usr/local/bin/julia

# Show Julia where conda libraries are \
RUN mkdir /etc/julia && \
    echo "push!(Libdl.DL_LOAD_PATH, \"${CONDA_DIR}/lib\")" >> /etc/julia/juliarc.jl && \
    # Create JULIA_PKGDIR \
    mkdir "${JULIA_PKGDIR}" && \
    chown "${NB_USER}" "${JULIA_PKGDIR}" && \
    fix-permissions "${JULIA_PKGDIR}"


# Copy the julia package into the container
RUN mkdir /LatticeQM
COPY Project.toml /LatticeQM/Project.toml
COPY src /LatticeQM/src

USER $NB_UID

RUN fix-permissions /LatticeQM && \
    fix-permissions "${CONDA_DIR}" && \
    fix-permissions "/home/${NB_USER}"


RUN julia -e 'import Pkg; Pkg.update()' && \
    julia -e 'import Pkg; Pkg.develop(path="/LatticeQM"); Pkg.add("Plots"); Pkg.add("ProgressMeter"); Pkg.add("IJulia"); Pkg.build(); Pkg.precompile();' && \
    fix-permissions /home/$NB_USER
