double mS = 1.192642, mS1385 = 1.38283;
struct Resonance{
    TString name, family;
    double spin;
    int parity;
    double mass, width;
    double g_kn, g_kxi;
};
//Resonance Lambda1890 = ResonanceList["Lambda1890"];
map<TString, Resonance> ResonanceLists =
{
    // key, {key, family, spin, parity, mass,width, g_kn, g_kxi } 
    {"u_Lambda", {"Lambda", "Lambda", 1./2., +1, 1.115683, 0.0, -13.24, 3.52}},
    {"u_Sigma", {"Sigma", "Sigma", 1./2., +1, 1.192642, 0.0, 3.58, -13.25}},
    {"Lambda1890", {"Lambda1890", "Lambda", 3./2., 1, 1.890, 0.120, 0.84, -0.26}},
    {"Lambda2100", {"Lambda2100", "Lambda", 7./2., -1, 2.100, 0.200, 2.41, 2.95}},
    {"Sigma2030", {"Sigma2030", "Sigma", 7./2., +1, 2.030, 0.180, 0.82, -0.93}},
    {"Sigma2230", {"Sigma2230", "Sigma", 3./2., +1, 2.230, 0.345, 0.41, 0.34}},
    {"Sigma2250", {"Sigma2250", "Sigma", 7./2., -1, 2.290, 0.100, 0.59, 0.88}}
};
