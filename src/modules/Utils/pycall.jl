
import PyCall # adds ~ 1.5 seconds to the load time 

function __init__()
    global brentq = PyCall.pyimport_conda("scipy.optimize", "scipy").brentq
    global cKDTree = PyCall.pyimport("scipy.spatial").cKDTree
end
