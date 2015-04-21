# ~/Arachne/jumping/findLibrary.py *.xml ${3}${1} $1.names
# JoinClones READS=$2 NAMES=$1.names OUT_PREFIX=$1 > $1.joinclones.out
#
MakeLookupTable SOURCE= $1.joined.fasta OUT_HEAD= $1.joined LOOKUP_ONLY=True IGNORE_SHORT=True
QueryLookupTable L= $1.joined.lookup  SEQS=stuffer.fasta MF=1000000 K=12 KB=0 MM=20 PARSEABLE=True > $1.joined.to.stuffer.qltout
RemoveStuffer QLTOUT= $1.joined.to.stuffer.qltout READS= $1.joined.fasta OUT= $1.joined.split.fasta  > $1.joined.removestuffer.out
QueryLookupTable L=../mouse35.lookup SEQS=$1.joined.split.fasta MF=10000 K=12 KB=0 MM=20 PARSEABLE=True VISUAL=True > $1.joined.split.vs.mouse.qltout 
EvaluateAlignmentPairs READS= $1.joined.split.fasta QLTOUT= $1.joined.split.vs.mouse.qltout > $1.joined.split.eval

MakeLookupTable SOURCE= $1.partial.fasta OUT_HEAD= $1.partial LOOKUP_ONLY=True IGNORE_SHORT=True
QueryLookupTable L= $1.partial.lookup  SEQS=stuffer.fasta MF=1000000 K=12 KB=0 MM=20 PARSEABLE=True > $1.partial.to.stuffer.qltout
RemoveStuffer QLTOUT= $1.partial.to.stuffer.qltout READS= $1.partial.fasta OUT= $1.partial.split.fasta  > $1.partial.removestuffer.out
QueryLookupTable L=../mouse35.lookup SEQS=$1.partial.split.fasta MF=10000 K=12 KB=0 MM=20 PARSEABLE=True VISUAL=True > $1.partial.split.vs.mouse.qltout 
EvaluateAlignmentPairs READS= $1.partial.split.fasta QLTOUT= $1.partial.split.vs.mouse.qltout > $1.partial.split.eval

MakeLookupTable SOURCE= $1.unjoined.fasta OUT_HEAD= $1.unjoined LOOKUP_ONLY=True IGNORE_SHORT=True
QueryLookupTable L= $1.unjoined.lookup  SEQS=stuffer.fasta MF=1000000 K=12 KB=0 MM=20 PARSEABLE=True > $1.unjoined.to.stuffer.qltout
RemoveStuffer QLTOUT= $1.unjoined.to.stuffer.qltout READS= $1.unjoined.fasta OUT= $1.unjoined.split.fasta PAIRS=True > $1.unjoined.removestuffer.out
QueryLookupTable L=../mouse35.lookup SEQS=$1.unjoined.split.fasta MF=10000 K=12 KB=0 MM=20 PARSEABLE=True VISUAL=True > $1.unjoined.split.vs.mouse.qltout 
EvaluateAlignmentPairs READS= $1.unjoined.split.fasta QLTOUT= $1.unjoined.split.vs.mouse.qltout > $1.unjoined.split.eval
