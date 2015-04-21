///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2015) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Define hardcoded directories for 10X experiments.

#include "CoreTools.h"

inline void SetTenxDirs( const int N, String& dir, String& odir, String& tdir )
{    dir = "/wga/scr4/jaffe/GapToy", tdir = "/wga/scr4/vendor/tenx";
     if ( N == 1 )
     {    dir += "/51400.newchem/a.final";
          odir = dir;
          tdir += "/2014-12-17";    }
     if ( N == 2 )
     {    dir += "/52009.HCC1143+BL/a.final";
          odir = dir + "/tenx_normal";
          tdir += "/2015-02-11/HCC1143_BL";    }
     if ( N == 3 )
     {    dir += "/52009.HCC1143+BL/a.final";
          odir = dir + "/tenx_tumor";
          tdir += "/2015-02-13";    }
     if ( N == 4 )
     {    dir += "/52398.XDP2/a.final";
          odir = dir + "/M";
          tdir += "/2015-02-11/S0117_XDP";    }
     if ( N == 5 )
     {    dir += "/52398.XDP2/a.final";
          odir = dir + "/F";
          tdir += "/2015-02-11/S0118_XDP";    }
     if ( N == 6 )
     {    dir += "/52398.XDP2/a.final";
          odir = dir + "/S";
          tdir += "/2015-02-11/S0119_XDP";    }
     {    Mkdir777(dir);
          Mkdir777(odir);    }
     tdir += "/parsed";
     Mkdir777(tdir);    }
