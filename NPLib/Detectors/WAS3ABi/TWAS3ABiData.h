#ifndef __WAS3ABiDATA__
#define __WAS3ABiDATA__
/*****************************************************************************
 * Copyright (C) 2009-2016    this file is part of the NPTool Project        *
 *                                                                           *
 * For the licensing terms see $NPTOOL/Licence/NPTool_Licence                *
 * For the list of contributors see $NPTOOL/Licence/Contributors             *
 *****************************************************************************/

/*****************************************************************************
 * Original Author: Adrien MATTA  contact address: a.matta@surrey.ac.uk      *
 *                                                                           *
 * Creation Date  : November 2012                                            *
 * Last update    :                                                          *
 *---------------------------------------------------------------------------*
 * Decription:                                                               *
 *  This class hold the WAS3ABi Silicon array raw data (Made for TIG64 card)   *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * Comment:                                                                  *
 *                                                                           *
 *                                                                           *
 *****************************************************************************/

// STL
#include <vector>

// ROOT
#include "TNamed.h"

class TWAS3ABiData : public TNamed {
private:
  // WAS3ABi
  // Energy
  std::vector<UShort_t>   fWAS3ABi_StripFront_DetectorNbr;
  std::vector<UShort_t>   fWAS3ABi_StripFront_StripNbr;
  std::vector<Double_t>   fWAS3ABi_StripFront_Energy;
  std::vector<Double_t>   fWAS3ABi_StripFront_TimeCFD;
  std::vector<Double_t>   fWAS3ABi_StripFront_TimeLED;
  
  std::vector<UShort_t>   fWAS3ABi_StripBack_DetectorNbr;
  std::vector<UShort_t>   fWAS3ABi_StripBack_StripNbr;
  std::vector<Double_t>   fWAS3ABi_StripBack_Energy;
  std::vector<Double_t>   fWAS3ABi_StripBack_TimeCFD;
  std::vector<Double_t>   fWAS3ABi_StripBack_TimeLED;

public:
  TWAS3ABiData();
  virtual ~TWAS3ABiData();
  
  void Clear();
  void Clear(const Option_t*) {};
  void Dump() const;
  
  /////////////////////           SETTERS           ////////////////////////
  inline void SetFront_DetectorNbr(const UShort_t& DetNbr)  {fWAS3ABi_StripFront_DetectorNbr.push_back(DetNbr);}
  inline void SetFront_StripNbr(const UShort_t& StripNbr)   {fWAS3ABi_StripFront_StripNbr.push_back(StripNbr);}
  inline void SetFront_Energy(const Double_t& Energy)       {fWAS3ABi_StripFront_Energy.push_back(Energy);}
  inline void SetFront_TimeCFD(const Double_t& TimeCFD)     {fWAS3ABi_StripFront_TimeCFD.push_back(TimeCFD);}
  inline void SetFront_TimeLED(const Double_t& TimeLED)     {fWAS3ABi_StripFront_TimeLED.push_back(TimeLED);}

  inline void SetBack_DetectorNbr(const UShort_t& DetNbr)   {fWAS3ABi_StripBack_DetectorNbr.push_back(DetNbr);}
  inline void SetBack_StripNbr(const UShort_t& StripNbr)    {fWAS3ABi_StripBack_StripNbr.push_back(StripNbr);}
  inline void SetBack_Energy(const Double_t& Energy)        {fWAS3ABi_StripBack_Energy.push_back(Energy);}
  inline void SetBack_TimeCFD(const Double_t& TimeCFD)      {fWAS3ABi_StripBack_TimeCFD.push_back(TimeCFD);}
  inline void SetBack_TimeLED(const Double_t& TimeLED)      {fWAS3ABi_StripBack_TimeLED.push_back(TimeLED);}


  inline void SetFront(const UShort_t& DetNbr,const UShort_t& StripNbr,const Double_t& Energy,const Double_t& TimeCFD,const Double_t& TimeLED)	{
		SetFront_DetectorNbr(DetNbr);
		SetFront_StripNbr(StripNbr);
		SetFront_Energy(Energy);
		SetFront_TimeCFD(TimeCFD);
		SetFront_TimeLED(TimeLED);
	};
	inline void SetBack(const UShort_t &DetNbr,const UShort_t &StripNbr,const Double_t &Energy,const Double_t &TimeCFD,const Double_t &TimeLED)	{
		SetBack_DetectorNbr(DetNbr);
		SetBack_StripNbr(StripNbr);
		SetBack_Energy(Energy);
		SetBack_TimeCFD(TimeCFD);
		SetBack_TimeLED(TimeLED);
	};
  
  /////////////////////           GETTERS           ////////////////////////
  inline UShort_t GetFront_DetectorNbr(const unsigned int &i) const {return fWAS3ABi_StripFront_DetectorNbr[i];}//!
  inline UShort_t GetFront_StripNbr(const unsigned int &i)    const {return fWAS3ABi_StripFront_StripNbr[i];}//!
  inline Double_t GetFront_Energy(const unsigned int &i)      const {return fWAS3ABi_StripFront_Energy[i];}//!
  inline Double_t GetFront_TimeCFD(const unsigned int &i)     const {return fWAS3ABi_StripFront_TimeCFD[i];}//!
  inline Double_t GetFront_TimeLED(const unsigned int &i)     const {return fWAS3ABi_StripFront_TimeLED[i];}//!


  inline UShort_t GetBack_DetectorNbr(const unsigned int &i) const {return fWAS3ABi_StripBack_DetectorNbr[i];}//!
  inline UShort_t GetBack_StripNbr(const unsigned int &i)    const {return fWAS3ABi_StripBack_StripNbr[i];}//!
  inline Double_t GetBack_Energy(const unsigned int &i)      const {return fWAS3ABi_StripBack_Energy[i];}//!
  inline Double_t GetBack_TimeCFD(const unsigned int &i)     const {return fWAS3ABi_StripBack_TimeCFD[i];}//!
  inline Double_t GetBack_TimeLED(const unsigned int &i)     const {return fWAS3ABi_StripBack_TimeLED[i];}//!


  inline unsigned int GetMultiplicityFront() const {return fWAS3ABi_StripFront_DetectorNbr.size();}//!
  inline unsigned int GetMultiplicityBack()  const {return fWAS3ABi_StripBack_DetectorNbr.size();}//!
 
  ClassDef(TWAS3ABiData,1)  // WAS3ABiData structure
};

#endif
