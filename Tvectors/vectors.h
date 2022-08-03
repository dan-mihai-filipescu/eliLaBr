#ifndef ROOT_TVECTORS
#define ROOT_TVECTORS

#include <vector>

#define NMAXPEAKS 50
#define NMAXELEMENTS 40

using namespace std;
using std::vector;

//-------------------------------------------------------------
class Tvectors
{
protected:
  int nmax;
  int npeaks;
  std::vector<float> PeaksPos;
  std::vector<float> PeaksVal;

public:
  bool AvgTimeContract;
  Tvectors();
  Tvectors(int Nmax);
  Tvectors(const Tvectors &);
  virtual ~Tvectors();
  
  void Add(float PkPos, float PkVal);
  void Add(Tvectors vectors);
  void Set(int npks, std::vector<float> PksPos, std::vector<float> PksVal);
  float GetPeakPos(int index);
  float GetPeakVal(int index);
  inline int GetNpeaks() {return npeaks;}
  void Clear(Option_t *option="");
  void QuickSort();
  void QuickSort(int left, int right);
  void Contract(float usPackChan);
  
  bool operator==(const Tvectors& c) const { return this==&c;}
  bool operator<(const Tvectors& c) const { return this<&c;}
  Tvectors& operator=(const Tvectors& c);
  //friend TBuffer &operator<<(TBuffer &b, const Tvectors *vectors);

  ClassDef(Tvectors,2) // the Tvectors class
};

inline TBuffer &operator>>(TBuffer &buf,Tvectors *&obj)
{
   obj = new Tvectors();
   obj->Streamer(buf);
   return buf;
}

inline TBuffer &operator<<(TBuffer &buf, const Tvectors *obj)
{
   ((Tvectors*)obj)->Streamer(buf);
   return buf;
}

//-------------------------------------------------------------
class TSTLvector
{
protected:
  Int_t             fNMax;
  Int_t             fNelements;
  vector <Tvectors> fList;

public:
  bool AvgTimeContract;
  TSTLvector();
  TSTLvector(Int_t Nmax);
  void Clear(Option_t *option="");
  void ClearElements(Option_t *option="");
  virtual ~TSTLvector();
  
  void Add(Tvectors tvec);
  Tvectors Get(Int_t index);
  inline Int_t GetNelements() {return fNelements;}
  void Set(Int_t index, Tvectors VectorIn);
  void Contract(float usPackChan);

  ClassDef(TSTLvector,2) // STL vector of Tvectors class
};

#endif


