#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "TCanvas.h"
#include "TPie.h"
#include "TROOT.h"
#include "TColor.h"

using namespace std;

//------------------------------------------------------------
// Determine region tag from filename
//------------------------------------------------------------
string getRegionTag(const string& filename) {
    if (filename.find("eventcount") != string::npos) {
        return "SignalRegion";
    } else if (filename.find("BKG") != string::npos) {
        return "CR_MET-40-70";
    } else if (filename.find("VALID") != string::npos) {
        return "ValidationRegion";
    }
    return "UnknownRegion";
}

//------------------------------------------------------------
// Parse the log file and extract background contributions
//------------------------------------------------------------
map<string, map<string, float>> extractData(const string& filename) {

    ifstream file(filename);
    string line, channel;
    map<string, map<string, float>> data;   // channel → (background → value)

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return data;
    }

    // Two formats: CR/VR logs (Background/plot_bkg.C: "eg channel :", "elefakePho : 123") and the SR event
    // count (myResult/plot_eventct.C: "--- e+gamma channel ---", "e->gamma fake  : 123"); the latter labels
    // are mapped onto the former keys.
    const vector<pair<string, string>> labelToKey = {
        {"elefakePho", "elefakePho"}, {"e->gamma fake", "elefakePho"},
        {"jetfakePho", "jetfakePho"}, {"jet->gamma fake", "jetfakePho"},
        {"jetfakeLep", "jetfakeLep"}, {"fake lepton", "jetfakeLep"},
        {"VGamma", "VGamma"}, {"V+gamma", "VGamma"},
        {"Electroweak", "Electroweak"}};

    while (getline(file, line)) {

        if (line.find("eg channel") != string::npos || line.find("e+gamma channel") != string::npos) {
            channel = "eg";
        }
        else if (line.find("mg channel") != string::npos || line.find("mu+gamma channel") != string::npos) {
            channel = "mg";
        }
        if (channel.empty()) continue;

        string key;
        for (const auto& lk : labelToKey) {
            if (line.rfind(lk.first, 0) == 0) { key = lk.second; break; }   // label at the start of the line
        }
        size_t colon = line.find(':');
        if (key.empty() || colon == string::npos) continue;

        float value;
        stringstream ss(line.substr(colon + 1));
        ss >> value;

        if (!ss.fail()) {
            data[channel][key] = value;
        }
    }

    file.close();
    return data;
}

//------------------------------------------------------------
// Draw and save pie chart
//------------------------------------------------------------
void plotPieChart(map<string, float> values,
                  string year,
                  string channel,
                  string vfp_string,
                  string regionTag) {

    if (values.empty()) {
        cout << "No data available for channel: " << channel << endl;
        return;
    }

    vector<string> labels;
    vector<float> fractions;

    vector<int> colors = {
        kRed+1, kBlue+1, kGreen+2, kMagenta+2, kOrange+1,
        kCyan+2, kViolet+1, kSpring+4, kTeal+3, kAzure+2
    };

    float total = 0.0;
    for (const auto& entry : values) {
        labels.push_back(entry.first);
        fractions.push_back(entry.second);
        total += entry.second;
    }

    int n = fractions.size();
    float* valuesArray = fractions.data();

    TPie pie(
        (channel + "_pie_" + year).c_str(),
        (year + " " + vfp_string).c_str(),
        n,
        valuesArray
    );

    for (int i = 0; i < n; i++) {
        string label = labels[i] +
                       Form(" (%.1f%%)", 100.0 * fractions[i] / total);

        pie.SetEntryLabel(i, label.c_str());
        pie.SetEntryFillColor(i, colors[i % colors.size()]);
        pie.SetEntryRadiusOffset(i, -0.05);
    }

    pie.SetLabelFormat("%txt");
    pie.SetTextSize(0.03);
    pie.SetCircle(0.5, 0.5, 0.3);

    TCanvas canvas(
        ("c_" + channel + "_" + year).c_str(),
        (channel + " Channel").c_str(),
        800, 600
    );

    pie.Draw("");   // 2D pie

    string outputPath = "/uscms_data/d3/tmishra/Output/"
                      + channel + "_" + regionTag + "_" + year;

    if (!vfp_string.empty()) {
        outputPath += vfp_string;
    }

    outputPath += ".png";

    canvas.SaveAs(outputPath.c_str());
    cout << "Saved pie chart: " << outputPath << endl;
}

void plot_pie_chart(string filename, string year, string vfp_string) {

    gROOT->SetBatch(true);

    string regionTag = getRegionTag(filename);
    map<string, map<string, float>> data = extractData(filename);

    for (const auto& channelData : data) {
        plotPieChart(channelData.second,
                     year,
                     channelData.first,
                     vfp_string,
                     regionTag);
    }
}
