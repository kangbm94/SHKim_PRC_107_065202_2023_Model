#pragma once

#include <TString.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct DiffCrossSection {
    // Dataset name assigned to each input file/row.
    // This replaces the duplicated vector<TString> sqrtS in the question.
    std::vector<TString> name;
    std::vector<double> sqrtS;
    std::vector<double> cth;
    std::vector<double> dsig;
    std::vector<double> dsigErr;

    void Clear()
    {
        name.clear();
        sqrtS.clear();
        cth.clear();
        dsig.clear();
        dsigErr.clear();
    }

    std::size_t Size() const { return sqrtS.size(); }
};

struct DiffCrossSectionSorted {
    TString name;  // Generated key: BaseName_0, BaseName_1, ...
    double sqrtS = 0.0;
    std::vector<double> cth;
    std::vector<double> dsig;
    std::vector<double> dsigErr;
};

// Filled by SortDiffCrossSection().
inline std::map<TString, DiffCrossSectionSorted> DiffDatasets;

namespace DiffCrossSectionReaderDetail {

inline void ValidateSizes(const DiffCrossSection& diffs)
{
    const std::size_t n = diffs.sqrtS.size();
    if (diffs.name.size() != n ||
        diffs.cth.size() != n ||
        diffs.dsig.size() != n ||
        diffs.dsigErr.size() != n) {
        throw std::runtime_error(
            "DiffCrossSection contains vectors with different sizes");
    }
}

inline bool IsBlankOrComment(const std::string& line)
{
    const std::size_t first = line.find_first_not_of(" \t\r\n");
    return first == std::string::npos || line[first] == '#';
}

struct Row {
    double sqrtS = 0.0;
    double cth = 0.0;
    double dsig = 0.0;
    double dsigErr = 0.0;
};

}  // namespace DiffCrossSectionReaderDetail

// Reads TSV files with columns
// #sqrts  #cth  #Xsection  #XsectionErr  [#cthe_L  #cthe_R]
// The optional angular-error columns are ignored because they are not members
// of DiffCrossSection.
inline void LoadDiffCrossSections(const std::vector<TString>& files,
                                  const std::vector<TString>& names,
                                  DiffCrossSection& diffs)
{
    if (files.size() != names.size()) {
        throw std::invalid_argument(
            "LoadDiffCrossSections: files and names must have the same size");
    }

    diffs.Clear();

    for (std::size_t ifile = 0; ifile < files.size(); ++ifile) {
        std::ifstream input(files[ifile].Data());
        if (!input.is_open()) {
            throw std::runtime_error(
                "Cannot open differential-cross-section file: " +
                std::string(files[ifile].Data()));
        }

        std::string line;
        std::size_t lineNumber = 0;
        while (std::getline(input, line)) {
            ++lineNumber;

            // Remove a UTF-8 byte-order mark if one is present.
            if (lineNumber == 1 && line.size() >= 3 &&
                static_cast<unsigned char>(line[0]) == 0xEF &&
                static_cast<unsigned char>(line[1]) == 0xBB &&
                static_cast<unsigned char>(line[2]) == 0xBF) {
                line.erase(0, 3);
            }

            if (DiffCrossSectionReaderDetail::IsBlankOrComment(line)) {
                continue;
            }

            double sqrtS = 0.0;
            double cth = 0.0;
            double dsig = 0.0;
            double dsigErr = 0.0;
            std::istringstream parser(line);
            if (!(parser >> sqrtS >> cth >> dsig >> dsigErr)) {
                std::ostringstream message;
                message << "Cannot parse " << files[ifile].Data()
                        << ":" << lineNumber
                        << ". Expected at least four numeric columns.";
                throw std::runtime_error(message.str());
            }

            if (!std::isfinite(sqrtS) || !std::isfinite(cth) ||
                !std::isfinite(dsig) || !std::isfinite(dsigErr)) {
                std::ostringstream message;
                message << "Non-finite value in " << files[ifile].Data()
                        << ":" << lineNumber;
                throw std::runtime_error(message.str());
            }
            if (cth < -1.0 || cth > 1.0 || dsigErr < 0.0) {
                std::ostringstream message;
                message << "Unphysical value in " << files[ifile].Data()
                        << ":" << lineNumber;
                throw std::runtime_error(message.str());
            }

            diffs.name.push_back(names[ifile]);
            diffs.sqrtS.push_back(sqrtS);
            diffs.cth.push_back(cth);
            diffs.dsig.push_back(dsig);
            diffs.dsigErr.push_back(dsigErr);
        }
    }
}

// Backward-compatible alias preserving the capitalization in the question.
inline void LoadDIffCrossSections(const std::vector<TString>& files,
                                  const std::vector<TString>& names,
                                  DiffCrossSection& diffs)
{
    LoadDiffCrossSections(files, names, diffs);
}

// Groups rows first by the supplied dataset name and then by sqrt(s).
// For each base name, energies are numbered in ascending order:
// Name_0, Name_1, ... . The angular points within each group are sorted by cth.
inline void SortDiffCrossSection(const DiffCrossSection& diffs,
                                 double sqrtSTolerance = 1.0e-6)
{
    using DiffCrossSectionReaderDetail::Row;

    DiffCrossSectionReaderDetail::ValidateSizes(diffs);
    if (sqrtSTolerance < 0.0) {
        throw std::invalid_argument(
            "SortDiffCrossSection: sqrtSTolerance must be non-negative");
    }

    DiffDatasets.clear();

    std::map<TString, std::vector<Row>> rowsByName;
    for (std::size_t i = 0; i < diffs.Size(); ++i) {
        rowsByName[diffs.name[i]].push_back(
            Row{diffs.sqrtS[i], diffs.cth[i], diffs.dsig[i], diffs.dsigErr[i]});
    }

    for (auto& nameAndRows : rowsByName) {
        const TString& baseName = nameAndRows.first;
        std::vector<Row>& rows = nameAndRows.second;

        std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) {
            if (a.sqrtS != b.sqrtS) return a.sqrtS < b.sqrtS;
            return a.cth < b.cth;
        });

        std::vector<std::vector<Row>> energyGroups;
        for (const Row& row : rows) {
            if (energyGroups.empty() ||
                std::abs(row.sqrtS - energyGroups.back().front().sqrtS) >
                    sqrtSTolerance) {
                energyGroups.push_back(std::vector<Row>{row});
            } else {
                energyGroups.back().push_back(row);
            }
        }

        for (std::size_t igroup = 0; igroup < energyGroups.size(); ++igroup) {
            std::vector<Row>& group = energyGroups[igroup];
            std::sort(group.begin(), group.end(), [](const Row& a, const Row& b) {
                return a.cth < b.cth;
            });

            const TString sortedName = TString::Format(
                "%s_%zu", baseName.Data(), igroup);

            DiffCrossSectionSorted sorted;
            sorted.name = sortedName;

            double sqrtSSum = 0.0;
            for (const Row& row : group) sqrtSSum += row.sqrtS;
            sorted.sqrtS = sqrtSSum / static_cast<double>(group.size());

            sorted.cth.reserve(group.size());
            sorted.dsig.reserve(group.size());
            sorted.dsigErr.reserve(group.size());
            for (const Row& row : group) {
                sorted.cth.push_back(row.cth);
                sorted.dsig.push_back(row.dsig);
                sorted.dsigErr.push_back(row.dsigErr);
            }

            DiffDatasets.emplace(sortedName, std::move(sorted));
        }
    }
}

struct ForwardDifferentialCrossSection {
    vector<TString> name;
    vector<double> beam_mom;
    vector<double> sqrts;
    vector<double> dsig;
    vector<double> dsigerr;
};

ForwardDifferentialCrossSection OldFwdData = {
    vector<TString>(),
    vector<double>(),
    vector<double>(),
    vector<double>(),
    vector<double>()
};
double CalculateSqrtS(double beam_mom)
{
  constexpr double mK = 0.493677;     // charged-kaon mass [GeV/c^2]
  constexpr double mp = 0.938272081;  // proton mass [GeV/c^2]

  const double beam_energy = std::sqrt(beam_mom * beam_mom + mK * mK);
  const double s = mK * mK + mp * mp + 2.0 * mp * beam_energy;
  return std::sqrt(s);
}
bool LoadFwdData(const TString &filename)
{
  std::ifstream input(filename.Data());
  if (!input.is_open()) {
    std::cerr << "Cannot open " << filename << std::endl;
    return false;
  }

  // Clear any previously loaded entries.
  OldFwdData.name.clear();
  OldFwdData.beam_mom.clear();
  OldFwdData.sqrts.clear();
  OldFwdData.dsig.clear();
  OldFwdData.dsigerr.clear();

  std::string line;
  int line_number = 0;

  while (std::getline(input, line)) {
    ++line_number;

    if (line.empty() || line[0] == '#')
      continue;

    std::istringstream row(line);
    std::string experiment;
    double beam_mom = 0.0;
    double dsig = 0.0;
    double dsigerr = 0.0;

    // operator>> treats tabs and spaces as delimiters.
    if (!(row >> experiment >> beam_mom >> dsig >> dsigerr)) {
      // Skip the column-header line.
      if (experiment == "Experiment")
        continue;

      std::cerr << "Malformed row at line " << line_number
                << ": " << line << std::endl;
      return false;
    }

    OldFwdData.name.emplace_back(experiment.c_str());
    OldFwdData.beam_mom.push_back(beam_mom);
    OldFwdData.sqrts.push_back(CalculateSqrtS(beam_mom));
    OldFwdData.dsig.push_back(dsig);
    OldFwdData.dsigerr.push_back(dsigerr);
  }

  return true;
}
double CMToLabCth(double sqrts, double cth_cm)
{
    constexpr double mK  = 0.493677;
    constexpr double mp  = 0.938272081;
    constexpr double mXi = 1.32171;

    const double s = sqrts * sqrts;

    const double EK_beam =
        (s - mp*mp - mK*mK) / (2.0*mp);
    const double pK_beam =
        std::sqrt(EK_beam*EK_beam - mK*mK);

    const double beta  = pK_beam / (EK_beam + mp);
    const double gamma = (EK_beam + mp) / sqrts;

    const double EK_cm =
        (s + mK*mK - mXi*mXi) / (2.0*sqrts);

    const double lambda =
        s*s + std::pow(mK, 4) + std::pow(mXi, 4)
        - 2.0*s*mK*mK
        - 2.0*s*mXi*mXi
        - 2.0*mK*mK*mXi*mXi;

    const double p_cm = std::sqrt(lambda) / (2.0*sqrts);
    const double sth_cm =
        std::sqrt(std::max(0.0, 1.0 - cth_cm*cth_cm));

    const double pz_lab =
        gamma*(p_cm*cth_cm + beta*EK_cm);
    const double pt_lab = p_cm*sth_cm;

    return pz_lab/std::sqrt(pz_lab*pz_lab + pt_lab*pt_lab);
}
double LabToCMCth(double sqrts, double cth_lab, double mK = 0.493677, double mp = 0.938272081, double mXi = 1.32171)
{

    const double s = sqrts*sqrts;

    const double EK_beam =
        (s - mp*mp - mK*mK)/(2.0*mp);
    const double pK_beam =
        std::sqrt(EK_beam*EK_beam - mK*mK);

    const double Etot_lab = EK_beam + mp;
    const double beta     = pK_beam/Etot_lab;
    const double gamma    = Etot_lab/sqrts;

    const double A =
        0.5*(s + mK*mK - mXi*mXi);

    /*
      From:
        Etot_lab * E_K+ - pK_beam*pK+*cos(theta_lab) = A

      Squaring gives a quadratic equation for pK+ in the lab.
    */
    const double a =
        Etot_lab*Etot_lab
        - pK_beam*pK_beam*cth_lab*cth_lab;

    const double b =
        A*pK_beam*cth_lab;

    const double discriminant =
        b*b
        - a*(Etot_lab*Etot_lab*mK*mK - A*A);

    if (discriminant < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double pK_lab =
        (b + std::sqrt(discriminant))/a;

    if (pK_lab < 0.0)
        return std::numeric_limits<double>::quiet_NaN();

    const double EK_lab =
        std::sqrt(pK_lab*pK_lab + mK*mK);

    const double sth_lab =
        std::sqrt(std::max(0.0, 1.0 - cth_lab*cth_lab));

    // Inverse Lorentz boost of the outgoing K+ momentum.
    const double pz_cm =
        gamma*(pK_lab*cth_lab - beta*EK_lab);
    const double pt_cm =
        pK_lab*sth_lab;

    return pz_cm/std::sqrt(pz_cm*pz_cm + pt_cm*pt_cm);
}