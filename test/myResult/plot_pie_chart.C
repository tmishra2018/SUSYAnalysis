#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "TCanvas.h"
#include "TPie.h"

using namespace std;

// Function to parse the log file and extract background contributions
map<string, map<string, float>> extractData(const string& filename) {
    ifstream file(filename);
    string line, channel;
    map<string, map<string, float>> data;  // Channel -> (Background Type -> Value)

    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return data;
    }

    while (getline(file, line)) {
        if (line.find("eg channel") != string::npos) {
            channel = "eg";
        } else if (line.find("mg channel") != string::npos) {
            channel = "mg";
        }

        if (!channel.empty() && 
            (line.find("elefakePho") != string::npos || 
             line.find("jetfakePho") != string::npos || 
             line.find("jetfakeLep") != string::npos || 
             line.find("VGamma") != string::npos || 
             line.find("TTJets") != string::npos ||
             line.find("TTGJets") != string::npos ||
             line.find("WWG") != string::npos ||
             line.find("WZG") != string::npos ||
             line.find("WW") != string::npos ||
             line.find("WZ") != string::npos)) {
		//line.find("Rare") != string::npos)) {

            string key;
            float value;
            stringstream ss(line);
            ss >> key;  // Read the label
            ss.ignore(256, ':'); // Ignore up to ':'
            ss >> value; // Read the number

            if (!ss.fail()) {
                data[channel][key] = value;
            }
        }
    }
    file.close();
    return data;
}

// Function to plot and save a 2D pie chart with percentages
void plotPieChart(map<string, float> values, string year, string channel, string vfp_string) {
    if (values.empty()) {
        cout << "No data available for channel: " << channel << endl;
        return;
    }

    vector<string> labels;
    vector<float> fractions;
	vector<int> colors = {kRed+1, kBlue+1, kGreen+2, kMagenta+2, kOrange+1, 
                      kCyan+2, kViolet+1, kSpring+4, kTeal+3, kAzure+2};

    float total = 0.0;
    for (const auto& entry : values) {
        labels.push_back(entry.first);
        fractions.push_back(entry.second);
        total += entry.second;
    }

    int n = fractions.size();
    float* valuesArray = fractions.data();

    TPie pie((channel + "_pie_" + year).c_str(), (year + vfp_string).c_str(), n, valuesArray);

    for (int i = 0; i < n; i++) {
        string percentLabel = labels[i] + Form(" (%.1f%%)", (fractions[i] / total) * 100);
        pie.SetEntryLabel(i, percentLabel.c_str());
        pie.SetEntryFillColor(i, colors[i % colors.size()]);
        pie.SetEntryRadiusOffset(i, -0.05); // Move labels slightly inward for readability
    }
    
    pie.SetLabelFormat("%txt");  // Show only text (labels with percentages)
    pie.SetTextSize(0.03);      // Adjust text size
    pie.SetCircle(0.5, 0.5, 0.3); // Adjust size and position of the pie chart
//    pie.SetAngularOffset(45);

    TCanvas canvas(("c_" + channel + "_" + year).c_str(), (channel + " Channel").c_str(), 800, 600);
    pie.Draw("");  // 2D Pie Chart (no "3d" option)

    // Construct output filename with vfp_string
    string outputPath = "/uscms_data/d3/tmishra/Output/" + channel + "_background_" + year;
    if (!vfp_string.empty()) {
        outputPath += "_" + vfp_string;
    }
    outputPath += ".pdf";

    canvas.SaveAs(outputPath.c_str());
    cout << "Saved pie chart: " << outputPath << endl;
}

// Function to read the log file and create pie charts
void plot_pie_chart(string filename, string year, string vfp_string) {
    gROOT->SetBatch(true);
    map<string, map<string, float>> data = extractData(filename);

    for (const auto& channelData : data) {
        plotPieChart(channelData.second, year, channelData.first, vfp_string);
    }
}

