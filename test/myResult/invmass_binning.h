// =============================================================================
// invmass_binning.h
//
// Shared by pred_sig.C, pred_qcdBkg.C, pred_eleBkg_pT_nVtx.C,
// pred_jetBkg_pT_nVtx.C, pred_VGBkg.C, pred_rareBkg.C and plot_invmass.C so
// that the data, every background prediction and the plot all use the SAME
// invariant-mass binning and the SAME definition of the mass.
//
// Variable : m(l,gamma) of the selected lepton + photon pair of each event that
//            lands in a signal-region bin (Bin.findSignalBin(...) >= 0).
//              eg channel -> m(e,gamma)
//              mg channel -> m(mu,gamma)
//
// Mass definition: the trees only store (pT, eta, phi) of the lepton and the
// photon, so the mass is built in the massless approximation
//     m^2 = 2 pT_l pT_gamma ( cosh(eta_l - eta_gamma) - cos(phi_l - phi_gamma) )
// (lepton / photon masses are negligible at these energies).
//
// To change the binning, edit ONLY the three constants below and re-run all
// pred_*.C scripts (the plot script reads the binning back from this header).
// =============================================================================
#ifndef INVMASS_BINNING_H
#define INVMASS_BINNING_H

#include <cmath>
#include <algorithm>

// Uniform binning: nInvMassBins bins between invMassMin and invMassMax (GeV).
// Events above invMassMax are folded into the last bin (see foldInvMass()).
const int    nInvMassBins = 25;
const double invMassMin   = 0.0;
const double invMassMax   = 1000.0;

// Invariant mass (GeV) of the lepton + photon pair, massless approximation.
inline double calcLepPhoMass(double lepPt, double lepEta, double lepPhi,
                             double phoPt, double phoEta, double phoPhi)
{
	// cos() is 2pi-periodic, so dPhi does not need to be wrapped into [-pi,pi]
	double m2 = 2.0 * lepPt * phoPt * (std::cosh(lepEta - phoEta) - std::cos(lepPhi - phoPhi));
	return (m2 > 0.0) ? std::sqrt(m2) : 0.0;
}

// Fold the overflow into the last bin so that the sum over mass bins equals the
// total number of events in the signal region.
inline double foldInvMass(double m)
{
	return std::min(m, invMassMax - 1e-6);
}

#endif
