FROM gitpod/workspace-full

USER gitpod

# Install Julia
RUN sudo apt-get update && \
    sudo apt-get install -y wget curl && \
    curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel release 

RUN sudo /usr/bin/pip3 install scipy

# Install Julia packages
RUN julia -e 'ENV["PYTHON"] = "/usr/bin/python3"; using Pkg; Pkg.add("PyCall")' && \
    julia -e 'using Pkg; Pkg.add(["Plots", "IJulia"])'