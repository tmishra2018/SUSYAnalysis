# Submit masked-bin limit jobs for the revised Run 2 cards (2026-11-02). usage: python3 submit.py T5WG|T6WG
import sys, os, glob
model = sys.argv[1]; mode = "unblind"
CFG = {  # tag in limits_Run2_unblind_<tag>: bins masked (1-18 mu-gamma, 19-36 e-gamma)
    "masked_20261102":            "13 18 34 35 36",
    "masked_egamma_20261102":     " ".join(map(str, range(1, 19))),
    "masked_mgamma_20261102":     " ".join(map(str, range(19, 37))),
    "masked_MET_lt200_20261102":  "4 5 6 7 8 9 13 14 15 16 17 18 22 23 24 25 26 27 31 32 33 34 35 36",
    "masked_MET_200_400_20261102":"1 2 3 7 8 9 10 11 12 16 17 18 19 20 21 25 26 27 28 29 30 34 35 36",
    "masked_MET_gt400_20261102":  "1 2 3 4 5 6 10 11 12 13 14 15 19 20 21 22 23 24 28 29 30 31 32 33",
}
here = os.path.dirname(os.path.abspath(__file__)); os.chdir(here)
split = f"cardsSplitted_{model}"
os.system(f"rm -rf {split} && mkdir -p {split}")
cards = sorted(glob.glob(f"/eos/uscms/store/user/tmishra/CombinedLimit/{model}/cards/cards_Run2_{mode}/counting_*.txt"))
chunks = [cards[i:i+2] for i in range(0, len(cards), 2)]
for j, c in enumerate(chunks):
    open(f"{split}/x{j}", "w").write("\n".join(c) + "\n")
print(f"{model}: {len(cards)} cards, {len(chunks)} jobs per config")
for tag, bins in CFG.items():
    if "--dry" not in sys.argv: os.makedirs(f"/eos/uscms/store/user/tmishra/CombinedLimit/{model}/limits/limits_Run2_{mode}_{tag}", exist_ok=True)
    logs = f"logs/{model}_{tag}"; os.makedirs(logs, exist_ok=True)
    sub = f"""universe              = vanilla
Executable            = runLimitMasked.sh
request_memory        = 4096
Should_Transfer_Files = YES
WhenToTransferOutput  = ON_EXIT_OR_EVICT
Transfer_Input_Files  = {split}
x509userproxy         = $ENV(X509_USER_PROXY)
+JobFlavour           = "nextweek"
Output    = {logs}/job$(job).stdout
Error     = {logs}/job$(job).stderr
Log       = {logs}/jobs.log
Arguments = "$(job) {model} {mode} {tag} '{bins}'"
Queue job from seq 0 {len(chunks)-1} |
"""
    f = f"sub_{model}_{tag}.jdl"; open(f, "w").write(sub)
    if "--dry" not in sys.argv:
        os.system(f"condor_submit {f}")
