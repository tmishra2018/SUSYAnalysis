// apply_ISRWeight.C
// ===================
// Utility functions to apply ISR weights correctly.
//
//   Rec.1: Reads channel-specific weight file (mg or eg)
//   Rec.2: Uses gen bosonPt as the reweighting variable when available
//          Falls back to reco ISRJetPt if not
//   Rec.3: Uses process-specific weight file (ZGToLLG or WGToLNuG)
//   Rec.4: Rare backgrounds (TTG, WWG, WZG) get ISRWeight=1 with
//          a 100% systematic uncertainty flag, NOT the Zgamma weight
//   Rec.5: Returns shape-only weight by default; shape+norm available
//          via flag — normalization should be handled separately in
//          the background estimation fit
//
// USAGE in analysis_VGBkg.C or similar:
//
//   #include "apply_ISRWeight.C"
//
//   ISRWeightReader isr_mg_ZG, isr_mg_WG, isr_eg_ZG, isr_eg_WG;
//   isr_mg_ZG.init("mg","ZGToLLG","2018",false);
//   isr_mg_WG.init("mg","WGToLNuG","2018",false);
//   isr_eg_ZG.init("eg","ZGToLLG","2018",false);
//   isr_eg_WG.init("eg","WGToLNuG","2018",false);
//
//   // In event loop:
//   float w = isr_mg_ZG.getWeight(bosonPt_gen, mcType, "mg");

#ifndef APPLY_ISRWEIGHT_H
#define APPLY_ISRWEIGHT_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <map>

// ---- Process type classification ----------------------------------------
// Rec.4: categorise processes for correct weight assignment
enum ISRProcClass {
    ISR_ZGAMMA  = 0,  // ZGToLLG — use ZG weight file
    ISR_WGAMMA  = 1,  // WGToLNuG, WGJet* — use WG weight file
    ISR_RARE    = 2,  // TTGJets, WWG, WZG — ISRWeight=1, unc=100%
    ISR_UNKNOWN = 3
};

// Map mcType integers to ISRProcClass
// Adjust mcType values to match your MCType enum
ISRProcClass getISRProcClass(int mcType){
    switch(mcType){
        case 4:  return ISR_ZGAMMA;  // ZGInclusive
        case 0:  return ISR_WGAMMA;  // WGJetInclusive
        case 1:  return ISR_WGAMMA;  // WGJet40
        case 2:  return ISR_WGAMMA;  // WGJet130
        case 5:  return ISR_RARE;    // DYJets (no ISR correction needed)
        case 6:  return ISR_RARE;    // TTG
        case 7:  return ISR_RARE;    // WWG
        case 8:  return ISR_RARE;    // WZG
        case 9:  return ISR_RARE;    // TTJets
        default: return ISR_UNKNOWN;
    }
}

// ---- ISRWeightReader: loads a weight file and returns per-event weight ---
struct ISRWeightBin {
    double ptLow, ptHigh;
    double wShape;   // shape-only weight
    double wTotal;   // shape+norm weight
    double unc;      // statistical uncertainty
};

class ISRWeightReader {
public:
    std::string channel, procTag, yearTag;
    double normFactor;
    std::vector<ISRWeightBin> bins;
    bool loaded;

    ISRWeightReader() : normFactor(1.0), loaded(false) {}

    // Rec.1+3: init with channel and process tag
    void init(const std::string& ch, const std::string& proc,
              const std::string& year, bool applyNorm=false){
        channel   = ch;
        procTag   = proc;
        yearTag   = year;
        loaded    = false;
        bins.clear();
        normFactor = 1.0;

        std::string fname = Form("ISRweights_%s_%s_%s.txt",
                                  ch.c_str(), proc.c_str(), year.c_str());
        std::ifstream f(fname.c_str());
        if(!f.is_open()){
            std::cout << "[WARN] ISRWeightReader: cannot open " << fname << std::endl;
            std::cout << "       Will use ISRWeight=1 for all events." << std::endl;
            return;
        }

        std::string line;
        while(std::getline(f, line)){
            if(line.empty() || line[0]=='#') continue;
            if(line.find("normFactor") != std::string::npos){
                sscanf(line.c_str(), "normFactor  %lf", &normFactor);
                continue;
            }
            ISRWeightBin b;
            int n = sscanf(line.c_str(), "%lf %lf %lf %lf %lf",
                           &b.ptLow, &b.ptHigh, &b.wShape, &b.wTotal, &b.unc);
            if(n == 5) bins.push_back(b);
        }
        f.close();
        loaded = (bins.size() > 0);
        std::cout << "[ISRWeightReader] Loaded " << bins.size() << " bins from " << fname
                  << "  normFactor=" << normFactor << std::endl;
    }

    // Rec.2: ptVar = gen bosonPt if available, else reco ISRJetPt
    // Rec.4: rare backgrounds get weight=1 with unc flag
    // Rec.5: applyNorm controls whether normalization is included
    float getWeight(float ptVar, int mcTypeInt, bool applyNorm=false,
                    float* uncOut=nullptr) const {

        ISRProcClass cls = getISRProcClass(mcTypeInt);

        // Rec.4: rare backgrounds — return 1.0, set unc=1.0 (100%)
        if(cls == ISR_RARE || cls == ISR_UNKNOWN){
            if(uncOut) *uncOut = 1.0; // 100% uncertainty
            return 1.0;
        }

        if(!loaded){
            if(uncOut) *uncOut = 0.5;
            return 1.0;
        }

        if(ptVar > 799) ptVar = 799;
        for(const auto& b : bins){
            if(ptVar >= b.ptLow && ptVar < b.ptHigh){
                if(uncOut) *uncOut = b.unc;
                return applyNorm ? (float)b.wTotal : (float)b.wShape;
            }
        }
        // Overflow: use last bin
        if(!bins.empty()){
            const auto& b = bins.back();
            if(uncOut) *uncOut = b.unc;
            return applyNorm ? (float)b.wTotal : (float)b.wShape;
        }
        if(uncOut) *uncOut = 0.5;
        return 1.0;
    }

    float getNormFactor() const { return (float)normFactor; }
};


// ---- Convenience: global readers, initialised once ----------------------
// Instantiate these at the top of your analysis script then call
// initISRReaders() once before the event loop.

ISRWeightReader g_isr_mg_ZG, g_isr_mg_WG;
ISRWeightReader g_isr_eg_ZG, g_isr_eg_WG;

void initISRReaders(const std::string& yearTag, bool applyNorm=false){
    std::cout << "[ISR] Initialising ISR weight readers for year " << yearTag << std::endl;
    // Rec.1: separate per channel; Rec.3: separate per process
    g_isr_mg_ZG.init("mg", "ZGToLLG",  yearTag, applyNorm);
    g_isr_mg_WG.init("mg", "WGToLNuG", yearTag, applyNorm);
    g_isr_eg_ZG.init("eg", "ZGToLLG",  yearTag, applyNorm);
    g_isr_eg_WG.init("eg", "WGToLNuG", yearTag, applyNorm);
}

// Main accessor — call this in the event loop
// channel: "mg" or "eg"
// ptVar:   gen bosonPt if available (Rec.2), else reco ISRJetPt
// mcTypeInt: your MCType enum value
// applyNorm: Rec.5 — false by default (shape only); set true to include norm
float getISRWeight(const std::string& channel, float ptVar, int mcTypeInt,
                   bool applyNorm=false, float* uncOut=nullptr){

    ISRProcClass cls = getISRProcClass(mcTypeInt);

    // Rec.4: rare backgrounds
    if(cls == ISR_RARE || cls == ISR_UNKNOWN){
        if(uncOut) *uncOut = 1.0;
        return 1.0;
    }

    // Rec.1+3: select correct reader
    if(channel == "mg"){
        if(cls == ISR_ZGAMMA) return g_isr_mg_ZG.getWeight(ptVar, mcTypeInt, applyNorm, uncOut);
        if(cls == ISR_WGAMMA) return g_isr_mg_WG.getWeight(ptVar, mcTypeInt, applyNorm, uncOut);
    }
    else if(channel == "eg"){
        if(cls == ISR_ZGAMMA) return g_isr_eg_ZG.getWeight(ptVar, mcTypeInt, applyNorm, uncOut);
        if(cls == ISR_WGAMMA) return g_isr_eg_WG.getWeight(ptVar, mcTypeInt, applyNorm, uncOut);
    }

    if(uncOut) *uncOut = 0.5;
    return 1.0;
}

#endif // APPLY_ISRWEIGHT_H
