///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "amino/Amino.h"

AminoAcid AminoAcid::gTable[]
{
   {'A',"Ala","Alanine"},       // 0
   {'R',"Arg","Arginine"},      // 1
   {'N',"Asn","Asparagine"},    // 2
   {'D',"Asp","Aspartic acid"}, // 3
   {'C',"Cys","Cysteine"},      // 4
   {'E',"Glu","Glutamic acid"}, // 5
   {'Q',"Gln","Glutamine"},     // 6
   {'G',"Gly","Glycine"},       // 7
   {'H',"His","Histidine"},     // 8
   {'I',"Ile","Isoleucine"},    // 9
   {'L',"Leu","Leucine"},       //10
   {'K',"Lys","Lysine"},        //11
   {'M',"Met","Methionine"},    //12
   {'F',"Phe","Phenylalanine"}, //13
   {'P',"Pro","Proline"},       //14
   {'O',"Pyl","Pyrrolysine"},   //15
   {'U',"Sec","Selenocysteine"},//16
   {'S',"Ser","Serine"},        //17
   {'Z',"Stp","Stop"},          //18
   {'T',"Thr","Threonine"},     //19
   {'W',"Trp","Tryptophan"},    //20
   {'Y',"Tyr","Tyrosine"},      //21
   {'V',"Val","Valine"},        //22
   {'X',"Unk","Unknown"}        //23
};

unsigned AminoAcid::gCodonToAA[]
{ 11, 2,11, 2,  //AAN
  19,19,19,19,  //ACN
   1,17, 1,17,  //AGN
   9, 9,12, 9,  //AUN
   6, 8, 6, 8,  //CAN
  14,14,14,14,  //CCN
   1, 1, 1, 1,  //CGN
  10,10,10,10,  //CUN
   5, 3, 5, 3,  //GAN
   0, 0, 0, 0,  //GCN
   7, 7, 7, 7,  //GGN
  22,22,22,22,  //GUN
  18,21,18,21,  //UAN
  17,17,17,17,  //UCN
  18, 4,20, 4,  //UGN
  10,13,10,13 };//UUN

unsigned AminoAcid::gSymbolToAA[]
{ 23, 0,23, 4, 3, 5,13, 7,   // ?A?CDEFG
   8, 9,23,11,10,12, 2,15,   // HI?KLMNO
  14, 6, 1,17,19,16,22,20,   // PQRSTUVW
  23,21,18,23,23,23,23,23 }; // ?YZ?????
