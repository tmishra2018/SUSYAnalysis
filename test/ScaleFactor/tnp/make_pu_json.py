# Build egamma-tnp pileup correction (correctionlib JSON) from the analysis PU histograms:
# data = test/PUfiles/dataPU_<tag>.root ("pileup"), MC = PileUp_with_official_recipe/PUfiles/PUWeights_<tag>.root ("hPUTrue_MC")
import sys
import uproot
from egamma_tnp.utils.pileup import create_correction

tag, test_dir, out = sys.argv[1], sys.argv[2], sys.argv[3]
pu_mc = uproot.open(f"{test_dir}/PileUp_with_official_recipe/PUfiles/PUWeights_{tag}.root")["hPUTrue_MC"].values().astype(float)
create_correction(f"{test_dir}/PUfiles/dataPU_{tag}.root", pu_mc, outfile=out, normalize_pu_mc_array=True)
print("wrote", out)
