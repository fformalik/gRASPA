#pragma once
#include <iomanip>
//Shared by fxn_main.h (end of simulation) and write_data.h (WriteTMMCEvery, during simulation)//
static inline void Write_TMMC(size_t Cycle, Components& SystemComponents, size_t SystemIndex, std::string FileName = "TMMC_Statistics.data")
{
  std::ofstream textTMMCFile{};
  std::filesystem::path cwd = std::filesystem::current_path();

  std::string dirname="TMMC/System_" + std::to_string(SystemIndex) + "/";
  std::string fname  = dirname + "/" + FileName;

  std::filesystem::path directoryName = cwd /dirname;
  std::filesystem::path fileName = cwd /fname;
  std::filesystem::create_directories(directoryName);
  textTMMCFile = std::ofstream(fileName, std::ios::out);
  for(size_t i = 0; i < SystemComponents.NComponents.x; i++)
  {
    if(SystemComponents.Tmmc[i].DoTMMC)
    {
      textTMMCFile << "Component " << i << ": " << SystemComponents.MoleculeName[i] << " -> Updated " << SystemComponents.Tmmc[i].TMUpdateTimes << " times \n";
      textTMMCFile << "Min Macrostate : " << SystemComponents.Tmmc[i].MinMacrostate << '\n';
      textTMMCFile << "Max Macrostate : " << SystemComponents.Tmmc[i].MaxMacrostate << '\n';
      textTMMCFile << "Wang-Landau Factor : " << SystemComponents.Tmmc[i].WLFactor << '\n';
      textTMMCFile << "N NMol Bin CM[-1] CM[0] CM[1] WLBias ln_g TMBias lnpi Forward_lnpi Reverse_lnpi Histogram U_count U_mean_K U_var_K2" << '\n';
      for(size_t j = 0; j < SystemComponents.Tmmc[i].Histogram.size(); j++)
      {
        size_t N   = j / SystemComponents.Tmmc[i].nbinPerMacrostate;
        size_t bin = j % SystemComponents.Tmmc[i].nbinPerMacrostate;
        textTMMCFile << j << " " << N << " " << bin << " "; 
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].CMatrix[j].x << " " ; 
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].CMatrix[j].y << " " ;
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].CMatrix[j].z << " " ;
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].WLBias[j] << " " ; 
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].ln_g[j] << " "; 
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].TMBias[j] << " " ;
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].lnpi[j] << " "; 
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].forward_lnpi[j] << " " ;
        textTMMCFile /*<< std::setprecision (15)*/ << SystemComponents.Tmmc[i].reverse_lnpi[j] << " " ; 
        textTMMCFile << SystemComponents.Tmmc[i].Histogram[j] << " ";
        //energy statistics of the bin in K (internal energy unit * energy_to_kelvin), variance in K^2//
        {
          const double toK = 1.2027242847;
          const TMMC& T = SystemComponents.Tmmc[i];
          double cnt = T.EnergyCount.size() > j ? T.EnergyCount[j] : 0.0;
          double mu  = T.EnergyMean.size()  > j ? T.EnergyMean[j]  : 0.0;
          double var = (cnt > 1.0) ? T.EnergyM2[j] / (cnt - 1.0) : 0.0;
          textTMMCFile << std::setprecision(15) << cnt << " " << mu * toK << " " << var * toK * toK << std::setprecision(6) << '\n';
        }
      }
    }
  }
  textTMMCFile.close();
}
