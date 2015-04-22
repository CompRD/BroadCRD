# Sat May 15 2010

Hi,

You can run the "example.py" out of the box.  I have checked in .so file that
should allow you to run the example script.

To learn what's going on, run example.py, look at example.py, look at
scaffold.i, and make.sh.  

Here is how you can build "scaffold.a", needed to link _scaffolds.so:

Step 0: $ svn update
Step 1: rerun ./configure
Step 2: setenv / export $FPIC == 'yes'
Step 3: $ rm -rf obj
Step 4: $ make -j 8 WriteScaffoldGraph

When the compiler finishes building the objects, make will write out an archive
file, link the binary, and then remove the archive.  Doap!  It's the archive
we're after!  Replay the 'ar' command on your screen.  Mine looks like this:

$ ar -qc obj/lib_WriteScaffoldGraph_temp.a obj/Alignment.o obj/AnnotatedContig.o obj/Basevector.o obj/BlockAlign.o obj/Charvector.o obj/Equiv.o obj/FastIfstream.o obj/Fastavector.o obj/IndexedAlignmentPlusVector.o obj/Intvector.o obj/MemberOf.o obj/Misc.o obj/PackAlign.o obj/PairsManager.o obj/ParseRange.o obj/ParseSet.o obj/PrintAlignment.o obj/Qualvector.o obj/ReadPairing.o obj/ScoreAlignment.o obj/SeqInterval.o obj/Superb.o obj/TokenizeString.o obj/Vec.o obj/VecAlignmentPlus.o obj/dna/Bases.o obj/feudal/BaseVec.o obj/feudal/BinaryStream.o obj/feudal/BitVec.o obj/feudal/CharString.o obj/feudal/FeudalControlBlock.o obj/feudal/FeudalFileReader.o obj/feudal/FeudalFileWriter.o obj/feudal/FeudalTools.o obj/feudal/Generic.o obj/feudal/Mempool.o obj/feudal/Oob.o obj/graph/Digraph.o obj/layout/ContigActualloc.o obj/lookup/LookAlign.o obj/math/Arith.o obj/math/Calculator.o obj/math/Functions.o obj/math/HoInterval.o obj/math/Matrix.o obj/math/Permutation.o obj/pairwise_aligners/Mutmer.o obj/paths/BuildScaffoldGraph.o obj/paths/Sepdev.o obj/paths/WriteScaffoldGraph.o obj/paths/reporting/CLinkBundle.o obj/paths/reporting/COffset.o obj/paths/reporting/CSuperLinks.o obj/random/Random.o obj/system/Assert.o obj/system/Crash.o obj/system/Exit.o obj/system/LockedData.o obj/system/MemTracker.o obj/system/ParsedArgs.o obj/system/ProcBuf.o obj/system/RunTime.o obj/system/SysConf.o obj/system/System.o obj/system/ThreadsafeIO.o obj/system/TraceVal.o obj/system/Types.o obj/system/UseGDB.o obj/system/WorklistUtils.o obj/util/Kronos.o

This will make an archive file at 'obj/lib_WriteScaffoldGraph_temp.a'

Step 5: Replay the archive command to get back the archive (see above).
Step 6: copy and rename the archive: 
    $ cp obj/lib_WriteScaffoldGraph_temp.a swig/scaffold.a
Step 7: cd in the 'swig' directoryw
Step 8: run the 'make.sh' script:
    $ ./make.sh
Step 9: run the ./example.py script:
    $ ./example.py

Feedback welcome!


