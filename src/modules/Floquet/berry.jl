using ..Spectrum

Spectrum.getberry!(bands::Spectrum.BandData, h, ks, drive, M::Integer) = Spectrum.getberry_wf!(bands, Spectrum.wavefunctions(h, drive, M), ks)

