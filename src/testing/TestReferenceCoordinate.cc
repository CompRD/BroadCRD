#include "MainTools.h"
#include "dna/ReferenceCoordinate.h"

int main(int argc, char **argv) {
    RunTime();

    BeginCommandArguments;
    CommandArgument_String(REF);
    CommandArgument_StringSet(RANGES);
    EndCommandArguments;

    ReferenceCoordinateMapper rcm(REF);
    vec<ReferenceRange> ranges = rcm.ParseRanges(RANGES);

    for (unsigned int rangeindex = 0; rangeindex < ranges.size(); rangeindex++) {
        cout << rangeindex << " " << ranges[rangeindex].AsGenomicRange() << " " << ranges[rangeindex].AsLinearRange() << endl;
    }

    return 0;
}
