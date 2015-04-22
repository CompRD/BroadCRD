#include "lookup/LookAlign.h"
class SerialQltout {
    public:
        SerialQltout() {};
        SerialQltout(String filename, bool sort = false) { 
            Open(filename, sort);
        }

        ~SerialQltout() { Close(); }

        void Open(String filename, bool sort = false) {
            if (sort) {
                // side-effect: sqltout defined
                SortByReference(filename);
                ifqltout = new ifstream(sqltout->c_str());
            } else {
                sqltout = NULL;
                ifqltout = new ifstream(filename.c_str());
            }     
        }

        void SortByReference(String filename) {
            String tmpDir = temp_file::tmpDir();
            sqltout = new temp_file(tmpDir+"/Coverage_qltout_XXXXXX");

            SystemSucceed("grep \"^QUERY\" " + filename + " | sort -n -k7 -k8 > " + (*sqltout));
        }

        bool Eof(void) {
            return ifqltout->fail();
        }

        bool Next(look_align &a, unsigned int *length = NULL) {
            String qltline;
            string temp;

            do {
                getline(*ifqltout, temp);
                qltline = temp;
            } while (!qltline.Contains("QUERY") && !ifqltout->fail());

            if (ifqltout->fail()) {
                return 0;
            }

            a.ReadParseable(qltline);

            if (length != NULL) {
                (*length) = qltline.size();
            }

            return 1;
        }

        bool NextSet(vec<look_align> &aligns) {
            look_align lafirst;
            aligns.clear();

            if (Next(lafirst)) {
                aligns.push_back(lafirst);

                look_align temp;
                unsigned int length = 0;
                while (Next(temp, &length)) {
                    if (temp.query_id != lafirst.query_id) {
                        ulonglong pos = ifqltout->tellg();
                        ifqltout->seekg(pos - length - 1);
                        break;
                    }

                    aligns.push_back(temp);
                }

                return 1;
            }

            return 0;
        }

        void Close(void) {
            if (ifqltout != NULL) {
                ifqltout->close();
                delete ifqltout;
                ifqltout = NULL;
            }

            if (sqltout != NULL) {
                delete sqltout;
                sqltout = NULL;
            }
        }

    private:
        temp_file *sqltout;
        ifstream *ifqltout;
};
