module DummySave

    # using HDF5

    export save, savedlm

    function save!() # file should be hdf5 output stream
        error("Should not be called directly!")
    end

    function save()
        error("Should not be called directly!")
    end

    function savedlm!() # file should be hdf5 output stream
        error("Should not be called directly!")
    end

    function savedlm()
        error("Should not be called directly!")
    end

    function save_wrapper(data, filename)
        mkpath(dirname(filename)) # make sure the path exists (create if need be)
        rm(filename, force=true) # make sure the file does not already exist (delete if need be)

        h5open(filename, "w") do file
            save!(file, data)
        end

        filename
    end


end
