///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <dirent.h>
#include "MainTools.h"
#include "system/file/Directory.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "util/Fastq.h"
#include "paths/long/magic/NanoData.h"
#include "paths/long/magic/SimpleDataFrame.h"


namespace {
    using fastq::ReadFastq;
    void ReadOneFastq(String const& filename,
            basevector& bases_out, qualvector& quals_out,
            const char phred_offset) {
        vecbasevector vbases;
        vecqualvector vquals;
        size_t nread = ReadFastq( filename, vbases, vquals, phred_offset );
        if ( nread != 1 ) FatalErr("Expected one read in Fastq " + filename
                    + ", but saw " + ToString(nread) );

        bases_out = vbases[0];
        quals_out = vquals[0];
    }
};

int main( int argc, char* argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String_OrDefault(CONSENSUS_PREFIX,  "cons_");
     CommandArgument_String_OrDefault(TEMPLATE_PREFIX,   "temp_");
     CommandArgument_String_OrDefault(COMPLEMENT_PREFIX, "comp_");
     CommandArgument_Int_OrDefault_Doc(VER,2,"format version, default=2 for R7.2+");
     CommandArgument_String(FASTQ_DIR);
     CommandArgument_String(OUT_HEAD);
     EndCommandArguments;

     ForceAssert(VER==1 || VER==2);

     String const fastq = ".fastq";
     String const model = ".model";
     String const events = ".events";
     String const align = ".align";
     String const hairpin = ".hairpin";

     vec<String> consensae = AllFiles( FASTQ_DIR );
     // delete anything not a consensus read
     consensae.erase(
             std::remove_if( consensae.begin(), consensae.end(),
                     [CONSENSUS_PREFIX,TEMPLATE_PREFIX,COMPLEMENT_PREFIX,fastq](String const& s)
                     { auto base = Basename(s);
                             return (
                             !s.EndsWith(fastq) ||
                             !base.StartsWith(CONSENSUS_PREFIX) ||
                             base.StartsWith(TEMPLATE_PREFIX) ||
                             base.StartsWith(COMPLEMENT_PREFIX) ); } ),
                     consensae.end() );

     const char phred_offset = 33;

     NanoData nano;

     size_t count = 0;
     for ( auto const& consensus : consensae ) {
         String const base = Basename(consensus).After(CONSENSUS_PREFIX).Before(fastq);
         String const dir  = FASTQ_DIR + "/";   // AllFiles() returns files relative to ./ for some reason

         cout << "importing read " << count++ << " for base = " << base << endl;

         basevector bases_out;
         qualvector quals_out;

         try {
             ReadOneFastq(dir + CONSENSUS_PREFIX + base + fastq, bases_out, quals_out, phred_offset);
             nano.AddPartial(NanoData::CONS, bases_out, quals_out );
             ReadOneFastq(dir + TEMPLATE_PREFIX + base + fastq, bases_out, quals_out, phred_offset);
             nano.AddPartial(NanoData::TEMP, bases_out, quals_out );
             ReadOneFastq(dir + COMPLEMENT_PREFIX + base + fastq, bases_out, quals_out, phred_offset);
             nano.AddPartial(NanoData::COMP, bases_out, quals_out );
             nano.AddRead();
         } catch ( std::runtime_error const& e ) {
             FatalErr(e.what());
         }

         String modelfile = dir + COMPLEMENT_PREFIX + base + model;
         if ( VER == 1 ) nano.AddModel(SimpleDataFrame( modelfile, "s,i,f,f,f,f,f" ), NanoData::COMP);
         else nano.AddModel(SimpleDataFrame(modelfile, "s,f,f,f,f,f" ), NanoData::COMP);
         modelfile = dir + TEMPLATE_PREFIX + base + model;
         if ( VER == 1 ) nano.AddModel(SimpleDataFrame( modelfile, "s,i,f,f,f,f,f" ), NanoData::TEMP);
         else nano.AddModel(SimpleDataFrame(modelfile, "s,f,f,f,f,f" ), NanoData::TEMP);

         String eventsfile = dir + COMPLEMENT_PREFIX + base + events;
         nano.AddEvents(SimpleDataFrame( eventsfile, "f,f,f,f,s,f,i,f,s,f,f,f,f,f"), NanoData::COMP);
         eventsfile = dir + TEMPLATE_PREFIX + base + events;
         nano.AddEvents(SimpleDataFrame( eventsfile, "f,f,f,f,s,f,i,f,s,f,f,f,f,f"), NanoData::TEMP);

         String alignfile = dir + CONSENSUS_PREFIX + base + align;
         nano.AddAligns(SimpleDataFrame( alignfile, "i,i,s" ));

         String hairpinfile = dir + base + hairpin;
         nano.AddHairpin(SimpleDataFrame( hairpinfile, "i,i"));
     }

     BinaryWriter::writeFile(OUT_HEAD+".nano", nano);

     NanoData nano2;
     BinaryReader::readFile(OUT_HEAD+".nano", &nano2);

     auto ev = nano2.GetEvents(NanoData::COMP);
     cout << "read back in " << ev.size() << " events" << endl;
     if ( ev.size() ) {
         cout << "last one is:" << endl;
         ev.back().Dump();
     }
     auto mo = nano2.GetModels(NanoData::COMP);
     cout << "read back in " << mo.size() << " models" << endl;
     if ( mo.size() ) {
         cout << "last one is: " << endl;
         mo.back().Dump();
     }
}
