// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME SelDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "AllPix_Hits_WriteToEntuple.h"
#include "AllPix_Frames_WriteToEntuple.h"
#include "allpix_dm.h"
#include "AllPixDigitAnimation.hh"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_FrameStruct(void *p);
   static void deleteArray_FrameStruct(void *p);
   static void destruct_FrameStruct(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FrameStruct*)
   {
      ::FrameStruct *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FrameStruct >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FrameStruct", ::FrameStruct::Class_Version(), "allpix_dm.h", 83,
                  typeid(::FrameStruct), DefineBehavior(ptr, ptr),
                  &::FrameStruct::Dictionary, isa_proxy, 4,
                  sizeof(::FrameStruct) );
      instance.SetDelete(&delete_FrameStruct);
      instance.SetDeleteArray(&deleteArray_FrameStruct);
      instance.SetDestructor(&destruct_FrameStruct);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FrameStruct*)
   {
      return GenerateInitInstanceLocal((::FrameStruct*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::FrameStruct*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FrameContainer(void *p = 0);
   static void *newArray_FrameContainer(Long_t size, void *p);
   static void delete_FrameContainer(void *p);
   static void deleteArray_FrameContainer(void *p);
   static void destruct_FrameContainer(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FrameContainer*)
   {
      ::FrameContainer *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FrameContainer >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FrameContainer", ::FrameContainer::Class_Version(), "allpix_dm.h", 29,
                  typeid(::FrameContainer), DefineBehavior(ptr, ptr),
                  &::FrameContainer::Dictionary, isa_proxy, 4,
                  sizeof(::FrameContainer) );
      instance.SetNew(&new_FrameContainer);
      instance.SetNewArray(&newArray_FrameContainer);
      instance.SetDelete(&delete_FrameContainer);
      instance.SetDeleteArray(&deleteArray_FrameContainer);
      instance.SetDestructor(&destruct_FrameContainer);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FrameContainer*)
   {
      return GenerateInitInstanceLocal((::FrameContainer*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::FrameContainer*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_WriteToNtuple(void *p);
   static void deleteArray_WriteToNtuple(void *p);
   static void destruct_WriteToNtuple(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::WriteToNtuple*)
   {
      ::WriteToNtuple *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::WriteToNtuple >(0);
      static ::ROOT::TGenericClassInfo 
         instance("WriteToNtuple", ::WriteToNtuple::Class_Version(), "AllPix_Frames_WriteToEntuple.h", 33,
                  typeid(::WriteToNtuple), DefineBehavior(ptr, ptr),
                  &::WriteToNtuple::Dictionary, isa_proxy, 4,
                  sizeof(::WriteToNtuple) );
      instance.SetDelete(&delete_WriteToNtuple);
      instance.SetDeleteArray(&deleteArray_WriteToNtuple);
      instance.SetDestructor(&destruct_WriteToNtuple);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::WriteToNtuple*)
   {
      return GenerateInitInstanceLocal((::WriteToNtuple*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::WriteToNtuple*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_SimpleHits(void *p = 0);
   static void *newArray_SimpleHits(Long_t size, void *p);
   static void delete_SimpleHits(void *p);
   static void deleteArray_SimpleHits(void *p);
   static void destruct_SimpleHits(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleHits*)
   {
      ::SimpleHits *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimpleHits >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimpleHits", ::SimpleHits::Class_Version(), "AllPix_Hits_WriteToEntuple.h", 33,
                  typeid(::SimpleHits), DefineBehavior(ptr, ptr),
                  &::SimpleHits::Dictionary, isa_proxy, 4,
                  sizeof(::SimpleHits) );
      instance.SetNew(&new_SimpleHits);
      instance.SetNewArray(&newArray_SimpleHits);
      instance.SetDelete(&delete_SimpleHits);
      instance.SetDeleteArray(&deleteArray_SimpleHits);
      instance.SetDestructor(&destruct_SimpleHits);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleHits*)
   {
      return GenerateInitInstanceLocal((::SimpleHits*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::SimpleHits*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_AllPixDigitAnimation(void *p);
   static void deleteArray_AllPixDigitAnimation(void *p);
   static void destruct_AllPixDigitAnimation(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AllPixDigitAnimation*)
   {
      ::AllPixDigitAnimation *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AllPixDigitAnimation >(0);
      static ::ROOT::TGenericClassInfo 
         instance("AllPixDigitAnimation", ::AllPixDigitAnimation::Class_Version(), "AllPixDigitAnimation.hh", 28,
                  typeid(::AllPixDigitAnimation), DefineBehavior(ptr, ptr),
                  &::AllPixDigitAnimation::Dictionary, isa_proxy, 4,
                  sizeof(::AllPixDigitAnimation) );
      instance.SetDelete(&delete_AllPixDigitAnimation);
      instance.SetDeleteArray(&deleteArray_AllPixDigitAnimation);
      instance.SetDestructor(&destruct_AllPixDigitAnimation);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AllPixDigitAnimation*)
   {
      return GenerateInitInstanceLocal((::AllPixDigitAnimation*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::AllPixDigitAnimation*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr FrameStruct::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FrameStruct::Class_Name()
{
   return "FrameStruct";
}

//______________________________________________________________________________
const char *FrameStruct::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FrameStruct*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FrameStruct::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FrameStruct*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FrameStruct::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FrameStruct*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FrameStruct::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FrameStruct*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FrameContainer::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FrameContainer::Class_Name()
{
   return "FrameContainer";
}

//______________________________________________________________________________
const char *FrameContainer::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FrameContainer*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FrameContainer::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FrameContainer*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FrameContainer::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FrameContainer*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FrameContainer::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FrameContainer*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr WriteToNtuple::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *WriteToNtuple::Class_Name()
{
   return "WriteToNtuple";
}

//______________________________________________________________________________
const char *WriteToNtuple::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WriteToNtuple*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int WriteToNtuple::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::WriteToNtuple*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *WriteToNtuple::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WriteToNtuple*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *WriteToNtuple::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::WriteToNtuple*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr SimpleHits::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimpleHits::Class_Name()
{
   return "SimpleHits";
}

//______________________________________________________________________________
const char *SimpleHits::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimpleHits*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimpleHits::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimpleHits*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimpleHits::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimpleHits*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimpleHits::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimpleHits*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr AllPixDigitAnimation::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *AllPixDigitAnimation::Class_Name()
{
   return "AllPixDigitAnimation";
}

//______________________________________________________________________________
const char *AllPixDigitAnimation::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AllPixDigitAnimation*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int AllPixDigitAnimation::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AllPixDigitAnimation*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AllPixDigitAnimation::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AllPixDigitAnimation*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AllPixDigitAnimation::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AllPixDigitAnimation*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void FrameStruct::Streamer(TBuffer &R__b)
{
   // Stream an object of class FrameStruct.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FrameStruct::Class(),this);
   } else {
      R__b.WriteClassBuffer(FrameStruct::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_FrameStruct(void *p) {
      delete ((::FrameStruct*)p);
   }
   static void deleteArray_FrameStruct(void *p) {
      delete [] ((::FrameStruct*)p);
   }
   static void destruct_FrameStruct(void *p) {
      typedef ::FrameStruct current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FrameStruct

//______________________________________________________________________________
void FrameContainer::Streamer(TBuffer &R__b)
{
   // Stream an object of class FrameContainer.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(FrameContainer::Class(),this);
   } else {
      R__b.WriteClassBuffer(FrameContainer::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FrameContainer(void *p) {
      return  p ? new(p) ::FrameContainer : new ::FrameContainer;
   }
   static void *newArray_FrameContainer(Long_t nElements, void *p) {
      return p ? new(p) ::FrameContainer[nElements] : new ::FrameContainer[nElements];
   }
   // Wrapper around operator delete
   static void delete_FrameContainer(void *p) {
      delete ((::FrameContainer*)p);
   }
   static void deleteArray_FrameContainer(void *p) {
      delete [] ((::FrameContainer*)p);
   }
   static void destruct_FrameContainer(void *p) {
      typedef ::FrameContainer current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::FrameContainer

//______________________________________________________________________________
void WriteToNtuple::Streamer(TBuffer &R__b)
{
   // Stream an object of class WriteToNtuple.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(WriteToNtuple::Class(),this);
   } else {
      R__b.WriteClassBuffer(WriteToNtuple::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_WriteToNtuple(void *p) {
      delete ((::WriteToNtuple*)p);
   }
   static void deleteArray_WriteToNtuple(void *p) {
      delete [] ((::WriteToNtuple*)p);
   }
   static void destruct_WriteToNtuple(void *p) {
      typedef ::WriteToNtuple current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::WriteToNtuple

//______________________________________________________________________________
void SimpleHits::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimpleHits.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(SimpleHits::Class(),this);
   } else {
      R__b.WriteClassBuffer(SimpleHits::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimpleHits(void *p) {
      return  p ? new(p) ::SimpleHits : new ::SimpleHits;
   }
   static void *newArray_SimpleHits(Long_t nElements, void *p) {
      return p ? new(p) ::SimpleHits[nElements] : new ::SimpleHits[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimpleHits(void *p) {
      delete ((::SimpleHits*)p);
   }
   static void deleteArray_SimpleHits(void *p) {
      delete [] ((::SimpleHits*)p);
   }
   static void destruct_SimpleHits(void *p) {
      typedef ::SimpleHits current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::SimpleHits

//______________________________________________________________________________
void AllPixDigitAnimation::Streamer(TBuffer &R__b)
{
   // Stream an object of class AllPixDigitAnimation.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AllPixDigitAnimation::Class(),this);
   } else {
      R__b.WriteClassBuffer(AllPixDigitAnimation::Class(),this);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_AllPixDigitAnimation(void *p) {
      delete ((::AllPixDigitAnimation*)p);
   }
   static void deleteArray_AllPixDigitAnimation(void *p) {
      delete [] ((::AllPixDigitAnimation*)p);
   }
   static void destruct_AllPixDigitAnimation(void *p) {
      typedef ::AllPixDigitAnimation current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AllPixDigitAnimation

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 210,
                  typeid(vector<string>), DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 210,
                  typeid(vector<float>), DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 210,
                  typeid(vector<double>), DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = 0);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 210,
                  typeid(vector<TVector3>), DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const vector<TVector3>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TVector3>*)0x0)->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete ((vector<TVector3>*)p);
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] ((vector<TVector3>*)p);
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *maplEintcOintgR_Dictionary();
   static void maplEintcOintgR_TClassManip(TClass*);
   static void *new_maplEintcOintgR(void *p = 0);
   static void *newArray_maplEintcOintgR(Long_t size, void *p);
   static void delete_maplEintcOintgR(void *p);
   static void deleteArray_maplEintcOintgR(void *p);
   static void destruct_maplEintcOintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,int>*)
   {
      map<int,int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,int>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,int>", -2, "map", 96,
                  typeid(map<int,int>), DefineBehavior(ptr, ptr),
                  &maplEintcOintgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,int>) );
      instance.SetNew(&new_maplEintcOintgR);
      instance.SetNewArray(&newArray_maplEintcOintgR);
      instance.SetDelete(&delete_maplEintcOintgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOintgR);
      instance.SetDestructor(&destruct_maplEintcOintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,int> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<int,int>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,int>*)0x0)->GetClass();
      maplEintcOintgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOintgR(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) map<int,int> : new map<int,int>;
   }
   static void *newArray_maplEintcOintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) map<int,int>[nElements] : new map<int,int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOintgR(void *p) {
      delete ((map<int,int>*)p);
   }
   static void deleteArray_maplEintcOintgR(void *p) {
      delete [] ((map<int,int>*)p);
   }
   static void destruct_maplEintcOintgR(void *p) {
      typedef map<int,int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,int>

namespace {
  void TriggerDictionaryInitialization_SelDict_Impl() {
    static const char* headers[] = {
"AllPix_Hits_WriteToEntuple.h",
"AllPix_Frames_WriteToEntuple.h",
"allpix_dm.h",
"AllPixDigitAnimation.hh",
0
    };
    static const char* includePaths[] = {
"./include",
"/home/idarraga/analysis/root-6.04.02/include",
"/home/idarraga/allpix/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$allpix_dm.h")))  FrameStruct;
class __attribute__((annotate("$clingAutoload$allpix_dm.h")))  FrameContainer;
class __attribute__((annotate("$clingAutoload$AllPix_Frames_WriteToEntuple.h")))  WriteToNtuple;
class __attribute__((annotate("$clingAutoload$AllPix_Hits_WriteToEntuple.h")))  SimpleHits;
class __attribute__((annotate(R"ATTRDUMP(A geometry and event track visualization class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$AllPixDigitAnimation.hh")))  AllPixDigitAnimation;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "AllPix_Hits_WriteToEntuple.h"
#include "AllPix_Frames_WriteToEntuple.h"
#include "allpix_dm.h"
#include "AllPixDigitAnimation.hh"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AllPixDigitAnimation", payloadCode, "@",
"FrameContainer", payloadCode, "@",
"FrameStruct", payloadCode, "@",
"SimpleHits", payloadCode, "@",
"WriteToNtuple", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("SelDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_SelDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_SelDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_SelDict() {
  TriggerDictionaryInitialization_SelDict_Impl();
}
