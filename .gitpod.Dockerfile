FROM gitpod/workspace-full

USER gitpod

# Install Julia
RUN sudo apt-get update && \
    sudo apt-get install -y wget curl && \
    curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel release 