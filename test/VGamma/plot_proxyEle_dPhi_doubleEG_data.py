import os
import ROOT
import math
import itertools


def get_mass_inv(i_pho, i_ele):
    photon_vec = ROOT.TLorentzVector()
    electron_vec = ROOT.TLorentzVector()
    photon_vec.SetPtEtaPhiE(chain_in.phoEt[i_pho], chain_in.phoEta[i_pho],
                            chain_in.phoPhi[i_pho], chain_in.phoE[i_pho])
    electron_vec.SetPtEtaPhiE(chain_in.elePt[i_ele], chain_in.eleEta[i_ele],
                              chain_in.elePhi[i_ele], chain_in.eleEn[i_ele])
    return (photon_vec + electron_vec).M()


ROOT.gROOT.Macro(os.path.expanduser('rootlogon.C'))
sw = ROOT.TStopwatch()
sw.Start()

chain_in = ROOT.TChain('ggNtuplizer/EventTree')
#chain_in.Add('skim_medPho_looseEle_doubleEG_data.root')
chain_in.Add('root://store/group/lpcsusystealth/ggNtuple_leppho/skim-DoubleEG_ReReco2016B.root')
n_entries = chain_in.GetEntries()
print 'Total number of events: ' + str(n_entries)

file_out = ROOT.TFile('plot_proxyEle_dPhi_doubleEG_data.root', 'recreate')
tree_out = ROOT.TTree('VarTree', '')
samples = ['good', 'proxy']
variables = ['dPhi', 'mt']
branches = {}
for i in samples:
    for j in variables:
        branch_name = i + '_ele_' + j
        branches[branch_name] = ROOT.std.vector('float')()
        tree_out.Branch(branch_name, branches[branch_name])

for j_entry in range(n_entries):
    i_entry = chain_in.LoadTree(j_entry)
    if i_entry < 0:
        break
    nb = chain_in.GetEntry(j_entry)
    if j_entry % 10000 == 0:
        print 'Reading event ' + str(j_entry)

    for i in samples:
        for j in variables:
            branches[i + '_ele_' + j].clear()

    photon_index = []
    for i in range(chain_in.nPho):
        if (chain_in.phoEt[i] > 35.0
                and abs(chain_in.phoEta[i]) < 2.5
                and chain_in.phoIDbit[i]>>1&1
                and not chain_in.phohasPixelSeed[i]):
            photon_index.append(i)

    good_electron_index = []
    proxy_electron_index = []
    for i in range(chain_in.nEle):
        if (chain_in.elePt[i] > 25.0
                and abs(chain_in.eleEta[i]) < 2.5
                and chain_in.eleD0[i] < 0.2
                and chain_in.eleDz[i] < 0.5):
            if chain_in.eleIDbit[i]>>2&1 and chain_in.elePFMiniIso[i] < 0.2:
                good_electron_index.append(i)
            # The elif below isn't quite the correct requirement, but it's
            # close enough for testing purposes.
            elif chain_in.eleIDbit[i]>>1&1 and chain_in.elePFMiniIso[i] < 0.4:
                proxy_electron_index.append(i)

    max_mass_inv = 0.0
    for i in photon_index:
        electron_indices = [good_electron_index, proxy_electron_index]
        for j in itertools.chain(*electron_indices):
            m_inv = get_mass_inv(i, j)
            if m_inv > max_mass_inv:
                max_mass_inv = m_inv

    passes_pho_trigger = chain_in.HLTPho>>14&1 == 1

    if (photon_index
            and (good_electron_index or proxy_electron_index)
            and max_mass_inv > 95.0
            and passes_pho_trigger
            and chain_in.pfMET > 20.0
            and chain_in.pfMET < 100.0):
        electron_lists = {'good_ele': good_electron_index,
                          'proxy_ele': proxy_electron_index}
        for name, index in electron_lists.items():
            for i in index:
                dPhi = ROOT.TVector2.Phi_mpi_pi(chain_in.pfMETPhi
                                                - chain_in.elePhi[i])
                mt = math.sqrt(2.0 * chain_in.pfMET * chain_in.elePt[i]
                               * (1.0 - math.cos(dPhi)))
                branches[name + '_dPhi'].push_back(dPhi)
                branches[name + '_mt'].push_back(mt)
        tree_out.Fill()

file_out.Write()
file_out.Close()

sw.Stop()
print 'Real time: ' + str(sw.RealTime() / 60.0) + ' minutes'
print 'CPU time:  ' + str(sw.CpuTime() / 60.0) + ' minutes'
