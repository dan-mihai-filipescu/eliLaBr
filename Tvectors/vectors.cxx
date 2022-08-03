#include "TClass.h"
#include "vectors.h"
#include "Riostream.h"

namespace
  {
  struct Counter
    {
    std::string name;
    int count;
    Counter(const std::string& n) : name(n), count(0) {}
    ~Counter() { print(); }
    void print(const std::string& msg="")
      { cout << msg << " --- Counter: " << name << " " << count << endl; }
    };
  }

Counter TvectorsCount("Tvectors");

//-------------------------------------------------------------
ClassImp(Tvectors)
//-------------------------------------------------------------
Tvectors::Tvectors():AvgTimeContract(false)
  {
  nmax = NMAXPEAKS;
  PeaksPos.reserve(NMAXPEAKS);
  PeaksVal.reserve(NMAXPEAKS);
  npeaks = 0;
  TvectorsCount.count++;
  }

Tvectors::Tvectors(int Nmax):AvgTimeContract(false)
  {
  if (Nmax>NMAXPEAKS) Nmax = NMAXPEAKS;
  nmax = Nmax;
  PeaksPos.reserve(nmax);
  PeaksVal.reserve(nmax);
  npeaks = 0;
  TvectorsCount.count++;
  }

Tvectors::Tvectors(const Tvectors &vectors):AvgTimeContract(false)
  {
  nmax = vectors.nmax;
  Clear();
  PeaksPos.reserve(nmax);
  PeaksVal.reserve(nmax);
  npeaks = vectors.npeaks;
  if(npeaks>nmax) npeaks = nmax;
  for(int j=0; j<npeaks; j++)
    {
    PeaksPos.push_back(vectors.PeaksPos.at(j));
    PeaksVal.push_back(vectors.PeaksVal.at(j));
    }
  TvectorsCount.count++;
  }

Tvectors::~Tvectors()
  {
  TvectorsCount.count--;
  Clear();
  }

void Tvectors::Clear(Option_t *)
  {
  npeaks = 0;
  PeaksPos.erase(PeaksPos.begin(),PeaksPos.end());
  PeaksVal.erase(PeaksVal.begin(),PeaksVal.end());
  }

Tvectors& Tvectors::operator=(const Tvectors& vectors)
  {
  nmax = vectors.nmax;
  Clear();
  npeaks = vectors.npeaks;
  if(npeaks>nmax) npeaks = nmax;
  for(int j=0; j<npeaks; j++)
    {
    PeaksPos.push_back(vectors.PeaksPos.at(j));
    PeaksVal.push_back(vectors.PeaksVal.at(j));
    }
  return *this;
  }

void Tvectors::Add(float PkPos, float PkVal)
  {
  if(npeaks<nmax)
    {
    PeaksPos.push_back(PkPos);
    PeaksVal.push_back(PkVal);
    npeaks++;
    }
    else
    {
    printf("Vectors in Tvectors class are already at maximum number of elements (%d)!\n",nmax);
    printf("Cannot push_back any more!\n");
    }
  }

void Tvectors::Add(Tvectors vectors)
  {
  if(vectors.npeaks>0)
  for(int j=0; j<vectors.npeaks; j++)
    this->Add(vectors.PeaksPos.at(j),vectors.PeaksVal.at(j));
  }

void Tvectors::Set(int npks, std::vector<float> PksPos, std::vector<float> PksVal)
  {
  npeaks = npks;
  PeaksPos = PksPos;
  PeaksVal = PksVal;
  }

float Tvectors::GetPeakPos(int index)
  {
  if((0<=index)&&(index<=npeaks))
    return PeaksPos.at(index);
    else
    {
    printf("Cannot Get PeakPos!\n");
    printf("Index out of range!\n");
    return 0;
    }
  }

float Tvectors::GetPeakVal(int index)
  {
  if((0<=index)&&(index<=npeaks))
    return PeaksVal.at(index);
    else
    {
    printf("Cannot Get PeakVal!\n");
    printf("Index out of range!\n");
    return 0;
    }
  }

void Tvectors::QuickSort()
  {
  if(!(npeaks>0)) return;
  QuickSort(0,npeaks-1);
  }

void Tvectors::QuickSort(int left, int right)
  {
      int i = left, j = right;
      float PkPos, PkVal;
      float pivot = PeaksPos.at((left + right) / 2);

      /* partition */
      while (i <= j)
      {
            while (PeaksPos.at(i) < pivot)
                  i++;
            while (PeaksPos.at(j) > pivot)
                  j--;
            if (i <= j) {
                PkPos = PeaksPos.at(i); PkVal = PeaksVal.at(i);
                PeaksPos.at(i) = PeaksPos.at(j); PeaksVal.at(i) = PeaksVal.at(j);
                PeaksPos.at(j) = PkPos;   PeaksVal.at(j) = PkVal;
                  i++;
                  j--;
            }
      }

      /* recursion */
      if (left < j)
            QuickSort(left, j);
      if (i < right)
            QuickSort(i, right);
  }

void Tvectors::Contract(float usPackChan)
  {
  QuickSort();
  //QuickSort(0,npeaks-1);
  for(int j=0; j<(npeaks-1); j++)
    {
    if((PeaksPos.at(j+1)-PeaksPos.at(j))<=usPackChan)
      {
      if(AvgTimeContract) PeaksPos.at(j) = (PeaksPos.at(j)+PeaksPos.at(j+1)) / 2.;
      PeaksVal.at(j) = PeaksVal.at(j) + PeaksVal.at(j+1);
      for(int i=(j+1); i<(npeaks-1); i++)
        {
	PeaksPos.at(i) = PeaksPos.at(i+1);
	PeaksVal.at(i) = PeaksVal.at(i+1);
	}
      npeaks--;
      j--;
      PeaksPos.pop_back();
      PeaksVal.pop_back();
      }
    }
  }

//-------------------------------------------------------------
ClassImp(TSTLvector)
//-------------------------------------------------------------
TSTLvector::TSTLvector():AvgTimeContract(false)
  {
  fNMax = NMAXELEMENTS;
  fList.reserve(fNMax);
  fNelements = 0;
  }

TSTLvector::TSTLvector(Int_t Nmax):AvgTimeContract(false)
  {
  if (Nmax>NMAXELEMENTS) Nmax = NMAXELEMENTS;
  fNMax = Nmax;
  fList.reserve(fNMax);
  fNelements = 0;
  }

void TSTLvector::Clear(Option_t *)
  {
  fNelements = 0;
  fList.erase(fList.begin(),fList.end());
  }

void TSTLvector::ClearElements(Option_t *)
  {
  for(int j=0; j<fNelements; j++) fList.at(j).Clear();
  }

TSTLvector::~TSTLvector()
  {
  Clear();
  TvectorsCount.print();
  }

void TSTLvector::Add(Tvectors tvec)
  {
  if(fNelements<fNMax)
    {
    fList.push_back(tvec);
    fNelements++;
    }
    else
    {
    printf("Vector in TSTLvector class is already at maximum number of elements (%d)!\n",fNMax);
    printf("Cannot push_back any more!\n");
    }  
  }

Tvectors TSTLvector::Get(Int_t index)
  {
  Tvectors emptyTvectors;
  if((0<=index)&&(index<=fNelements))
    return fList.at(index);
    else
    {
    printf("Cannot Get fList element!\n");
    printf("Index out of range!\n");
    return emptyTvectors;
    }  
  }

void TSTLvector::Set(Int_t index, Tvectors VectorIn)
  {
  if((0<=index)&&(index<=fNelements))
    fList.at(index).Add(VectorIn);
    else
    {
    printf("Cannot Set fList element!\n");
    printf("Index out of range!\n");
    }    
  }

void TSTLvector::Contract(float usPackChan)
  {
  for(int j=0; j<fNelements; j++)
    {
    fList.at(j).AvgTimeContract = AvgTimeContract;
    fList.at(j).Contract(usPackChan);
    }
  }
