///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * HeapTracker.cc
 *
 *  Created on: Apr 22, 2014
 *      Author: tsharpe
 */
#include "feudal/HeapUse.h"
#include "feudal/HashSet.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <X11/Xlib.h>
// MakeDepend: archived
// MakeDepend: library X11
namespace
{
class HeapUseHasher
{
public:
    typedef HeapUse argument_type;
    size_t operator()( HeapUse const& heapUse ) const
    { return heapUse.getAddr(); }
};
typedef HashSet<HeapUse,HeapUseHasher> HeapUseSet;

class XWin
{
public:
    XWin()
    { mDisplay = XOpenDisplay(getenv("DISPLAY"));
      if ( !mDisplay )
      { std::cout << "Couldn't open display." << std::endl; exit(1); }
      mBlack = BlackPixel(mDisplay,DefaultScreen(mDisplay));
      mWhite = WhitePixel(mDisplay,DefaultScreen(mDisplay));
      mWin = XCreateSimpleWindow(mDisplay,DefaultRootWindow(mDisplay),
                                  0,0,1024,1024,5,mBlack,mBlack);
      XSelectInput(mDisplay,mWin,StructureNotifyMask);
      XMapWindow(mDisplay,mWin);
      mGC = XCreateGC(mDisplay,mWin,0,nullptr);
      XEvent e;
      do XNextEvent(mDisplay,&e);
      while ( e.type != MapNotify );
      XSelectInput(mDisplay,mWin,0l); }

    void paintPixel( int row, int col, bool on )
    { XSetForeground(mDisplay,mGC,on?mWhite:mBlack);
      XDrawPoint(mDisplay,mWin,mGC,col,1023-row);
      XFlush(mDisplay); }

private:
    Display* mDisplay;
    int mBlack;
    int mWhite;
    Window mWin;
    GC mGC;
};

// per Mb use stats for a Tb of address space
class InUseTracker
{
    static size_t const ONE_MEGA = 1024ul * 1024ul;
    static size_t const ONE_MEGA_MASK = ONE_MEGA - 1ul;
    static size_t const ONE_TERA = ONE_MEGA * ONE_MEGA;
    static size_t const ONE_TERA_MASK = ONE_TERA - 1ul;

public:
    InUseTracker( XWin& xWin ) : mXWin(xWin) {}

    void alloc( HeapUse const& heapUse )
    { size_t beg = heapUse.getAddr() & ONE_TERA_MASK;
      unsigned* pCount = &mInUse[beg>>20];
      size_t end = beg + heapUse.getSize();
      if ( end & ONE_TERA ) end = ONE_TERA;
      size_t use = (~beg & ONE_MEGA_MASK) + 1;
      while ( beg != end )
      { use = std::min(use,end-beg);
        if ( !*pCount ) pixel(beg>>20,true);
        *pCount++ += use;
        beg += use;
        use = ONE_MEGA; } }

    void dealloc( HeapUse const& heapUse )
    { size_t beg = heapUse.getAddr() & ONE_TERA_MASK;
      unsigned* pCount = &mInUse[beg>>20];
      size_t end = beg + heapUse.getSize();
      if ( end & ONE_TERA ) end = ONE_TERA;
      size_t use = (~beg & ONE_MEGA_MASK) + 1;
      while ( beg != end )
      { use = std::min(use,end-beg);
        if ( !(*pCount++ -= use) ) pixel(beg>>20,false);
        beg += use;
        use = ONE_MEGA; } }


private:
    void pixel( size_t offset, bool on )
    { mXWin.paintPixel(offset/1024,offset&0x3ff,on); }

    std::array<unsigned,ONE_MEGA> mInUse;
    XWin& mXWin;
};

}

int main( int arc, char** argv )
{
    HeapUseSet allocSet(10000000ul);
    XWin xWin;
    InUseTracker inUseTrackerLow(xWin);
    InUseTracker inUseTrackerHigh(xWin);
    HeapUse heapUse;

    std::cout << std::hex;
    std::ifstream in("/wga/scr4/tsharpe/memuse");
    while ( in >> heapUse )
    {
        if ( heapUse.getSize() )
        {
            allocSet.insertUniqueValueNoLocking(heapUse);
            size_t whichTb = heapUse.getAddr() >> 40;
            if ( !whichTb ) inUseTrackerLow.alloc(heapUse);
            else if ( whichTb == 0x7f ) inUseTrackerHigh.alloc(heapUse);
            else
                std::cout << "Ignoring alloc at "
                          << heapUse.getAddr() << std::endl;
        }
        else
        {
            HeapUse const* pHeapUse = allocSet.lookup(heapUse);
            if ( !pHeapUse )
            {
                std::cout << "Couldn't find alloc at "
                          << heapUse.getAddr() << std::endl;
                continue;
            }
            size_t whichTb = heapUse.getAddr() >> 40;
            if ( !whichTb ) inUseTrackerLow.dealloc(*pHeapUse);
            else if ( whichTb == 0x7f ) inUseTrackerHigh.dealloc(*pHeapUse);
            else
                std::cout << "Ignoring dealloc at "
                          << heapUse.getAddr() << std::endl;
            allocSet.removeNoLocking(*pHeapUse);
        }
    }
    if ( !in.eof() )
        std::cout << "Failed to read input to EOF.";
    for ( auto const& hhs : allocSet )
        for ( HeapUse const& hu : hhs )
            std::cout << hu.getAddr() << '\t'
                        << hu.getSize() << '\n';
}
