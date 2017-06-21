echo "Simple match of butane parameters in CHARMM-style FF."

../../charmm/match.py --prm cgenff.36.prm buta.psf mm-buta.npy out
../../charmm/write_prm.py out out/out.prm
cat out/out.prm
