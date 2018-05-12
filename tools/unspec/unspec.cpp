/* unspec, a program for extracting scans from data files created by spec
 * command-line driven, many useful options
 *
 * Copyright (C) 1994-2017 Petr Mikulik
 *
 * e-mail: mikulik@physics.muni.cz
 * web:    http://www.sci.muni.cz/~mikulik/
 *
 */
#define VERSION "18. 1. 2017"
/*
 * The latest version of the program can be found on the above-mentioned
 * web page.
 */

/*
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */


/*
   Compilation:
     GNU C/C++ under Linux, unixes, et al:
	g++ -o unspec unspec.cpp -O2 -s
     GNU C/C++ under Windows:
	g++ -o unspec.exe unspec.cpp -O2 -s
     GNU C/C++ under IBM OS/2 Warp:
	gcc -o unspec.exe unspec.cpp -O2 -s -lstdcpp
     Clang++ compiler:
	clang++ -o unspec unspec.cpp -O2
     etc. --- see the enclosed Makefile.
   Installation:
     Copy the executable into any directory on your PATH or make an alias
     for it.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#if 1 // new C++ norm
    #include <iostream>
    #include <fstream>
    #include <iomanip>
    using namespace std;
#else // old C++ norm
    #include <iostream.h>
    #include <fstream.h>
    #include <iomanip.h>
#endif


#define PSD // support for the "@A" spectrum of a PSD/MCA

#define MAX_OP_LINES 100 // max. number of #O and #P lines can be 99
#define MAX_OP_MOTORSonLINE 99 // no more than 99 motors on #O or #P line
// tech note: larger numbers would need change of (100*n+k) to (1000*n+k) earlier in the code

#define SCAN_LOWER -999999
#define SCAN_UPPER  999999
#define WRITE_COMMENTS "SDTPL"


#if defined(__unix) || defined(__VMS) || defined(__hpux) || defined(__APPLE__)
  #define stricmp strcasecmp   // at least Linux/gcc2.7.2, HP CC, DEC CXX
  #define strnicmp strncasecmp
#endif

#ifdef DEBUG
  #define OKK cout <<__LINE__ <<" in " <<__FUNCTION__ <<"(): ";
  #define OK cout <<__LINE__ <<" in " <<__FUNCTION__ <<"() is OK\n";
  #define OK_ENTER cout <<__LINE__ <<" in " <<__FUNCTION__ <<"() ENTERED\n";
  #define OK_EXIT cout <<__LINE__ <<" in " <<__FUNCTION__ <<"() EXITING\n";
  #define outOK  out <<"LINE " <<__LINE__<<" is OK.\n";
#else
  #define OK_ENTER
  #define OK_EXIT
#endif


/* Global constants set from the command line options.
 */
int scan_from = SCAN_LOWER, scan_to = SCAN_UPPER;
int use_scans_range = 0;
char *scan_substring_filter = 0;
int write_grid = 1;
  // option '-#'
char *WriteComments = strdup(WRITE_COMMENTS);
  // option '-c'
int single_file = 0;
  // option `-1`: output goes to a single file
int no_scannb_in_outfile = 0;
  // option `-2`: no scan number in the output file
int out2stdout = 0;
  // option `-`: output goes to stdout
int scan_digits = 2;
  // options '-3', '-4'
int reindex_scans = 0, new_scan_index = 1;
  // option '-i'
int blank_after_pts = 0;
  // option '-b'
int blank_between_scans = 1;
  // option '-B'
int write_header = 0;
  // option '-H'
#ifdef PSD
int write_psd = 0;
  // options '-Pn' and '-Py': enable/disable writing out PSD/MCA spectra
int psd_roi_from = -1, psd_roi_to = -1;
  // option '-Pr': write only selected ROI; numbering from 1..max_channels
int psd_every = 1;
  // option '-Pe': write only every psd_every-th channel
int psd_average = -1;
  // option '-Pa': average channel intensities over psd_average channels
int psd_write_rows = 0;
  // option '-PR': write PSD spectra in rows (default in columns)
int psd_write_aA = 0;
  // option '-P@A': write PSD spectra in columns on new lines starting @A
int psd_ang_scale = 0;
double psd_centre = 0, psd_chans_1deg = 1;
  // recalculation PSD channels -> angles if psd_ang_scale==1, for '-PR'
  // options '-P0', '-P1': channel of angle=0 and nb of channels/1 deg
char *psd_mca_file_prefix = 0;
//int psd_which_mcadev = -1;
int psd_which_mcadev = -1;
  // which @MCADEV goes to output?
  // -1   ... the first one (default)
  // >= 0 ... explicitly choosen to match @MCADEV 0, @MCADEV 1, etc.
#endif


// Global variables.
ofstream *fOut = 0;  // Output stream, if file.
ostream *Out = 0;    // Output stream, can be file or cout.
#define out *Out     // Make me happy.

char S[65530]; // This is the input line here -- with its max length!
int MaxScanNb = -2, ScanNb = -1, DataFileScanNb, PointNb;


// Variables for autoguessing motor names if option "-w =":
int want_autoguess_names = 0;	// is motor names autoguessing ("-w =,...") requested?
int autoguess_names = 0;    	// how many motor names were autoguessed
char autoguess_name[MAX_OP_LINES][16]; // usually no more than 3 names are autoguessed
int autoguess_which_column_namewhere[9]; // there is no scan/mesh with more than 8 motors moving together

#ifdef PSD
// Keep here the current PSD/MCA spectra.
class t_psd {
  public:
    char rawline[65535];   // this is max length of psd data line

    int channels;	   // counts 0..4095, while user gives 1..
    int *intens;           // integer intens from rawline, 0..4095
    double *aintens;       // averaged intensities, max 0..4095

    int curr_which_mcadev; // which @MCADEV to read in the current scan?
			   // User wants psd_which_mcadev, but it must be found
			   // in the scan header unless default.

    t_psd();		   // constructor
    void do_intensities(); // fill intens+aintens from rawline
  };

t_psd psd;


t_psd::t_psd()
{
intens = 0;
aintens = 0;
}


/* Transform from words to numerical values required for averaging or
   when calculating total roi intensity.
   Thus this parses words from rawline to an array of numbers, and does
   averaging if requested.
*/
void t_psd::do_intensities ()
{
if (!intens) intens = new int [4096];
if (!aintens) aintens = new double [4096];
channels = 0; // counts 0..4095, while users gives 1..
int i;
char *c;
for (c=rawline; *c; ) { // Convert words to integers.
  while (*c && isspace(*c)) c++; // skip spaces
  if (!*c) break;
  if (channels>4095) break; // out of range
  sscanf(c,"%i",&intens[channels]);
  channels++;
  while (*c && !isspace(*c)) c++; // skip the number
  }
if (psd_average > 1) { // Average the intensity into aintens.
  /* Averaging over ((pts-1)/2) points left and (pts/2) points right,
     e.g. for psd_average=4 averaging each i-1..i+2 and =5  i-2..i+2.
  */
  int sum = 0;
  int L = -(psd_average-1)/2, R = psd_average/2; // last left and right range points
  L--; R--;
  int pts = R+1;
  for (i=0; i<=R; i++) sum += intens[i];
  for (i=0; i<channels; i++) {
    if (++R<channels) { sum += intens[R]; pts++; }
    if (L>=0) { sum -= intens[L]; pts--; }
    L++;
    aintens[i] = (double)sum / pts;
    }
  }
}

#endif


/*
void char2double ( char *optarg, double& X )
{
char *c; X=strtod(optarg,&c);
if (c[0]) { cerr <<"Wrong number " <<optarg <<endl; exit(1); }
}


void char2int ( char *optarg, int& X )
{
char *c; double d; d=strtod(optarg,&c);
if (c[0]) { cerr <<"Wrong number " <<optarg <<endl; exit(1); }
X=int(d+0.5);
}
*/


void ArgsRequired ( int howmany, int optarg, int argc, char *argv[] )
{
if (optarg+howmany>=argc) {
  cerr <<"Hey, you! Option " <<argv[optarg] <<" requires " <<howmany <<" argument";
  if (howmany>1) cerr <<'s';
  cerr <<".\n";
  exit(1); }
}


/* Get rid of the first word of the string.
*/
char* stripword1 ( char *S )
{
char *c = S; for ( ; *c && !isspace(*c); c++) ;
	     for ( ; *c && isspace(*c); c++ ) ;
char *s = S; for (; (*s++ = *c++); ) ;
     *s = 0; // command on separate line - shutup clang++ warning
return S;
}


/* Return number of words in the string.
*/
int words ( char *S )
{
int w = 0;
for (char *c = S; *c; w++) {
  while (*c && isspace(*c)) c++;
  if (!*c) break;
  while (*c && !isspace(*c)) c++;
  }
return w;
}


/* Copy n-th word from the S to w, supposing there are at least 1..n words
   in S, otherwise 0 is returned. There must be enough memory allocated in w.
*/
void copy_word ( char *to, char *S, int n )
{
int w = 0; // word counter
char *c = S;
*to = 0;
if (n<=0) return;
while (1) { // find n-th word
  while (*c && isspace(*c)) c++; // skip spaces
  if (!*c) return; // less than n words
  if (++w==n) break; // n-th word found
  while (*c && !isspace(*c)) c++;
  if (!*c) return; // less than n words
  }
char *from = c;
while (*c && !isspace(*c)) c++;
n = c - from;
memcpy(to,from,n); *(to+n) = 0;
}


/* WORDINDEX in REXX: returns the position of the first character in the n-th
   blank-delimited word in string. n must be a positive whole number. If
   there are fewer than n words in the string, 0 is returned.
   Examples: WORDINDEX('Now is the time',3) -> 8,
	     WORDINDEX('Now is the time',6) -> 0.
   I.e., words are count from 1.   Here, in C: returns pointer to that char.
*/
char* wordindex ( char *S, int n )
{
if (!S || !*S || n<1) return 0;
int w = 0;
while (*S) {
  while (*S && isspace(*S)) S++;  // skip spaces
  if (!*S) return 0;
  w++;
  if (w==n) return S; // n-th word found
  while (*S && !isspace(*S)) S++; // skip current word
}
return 0;
}


/* Like wordindex, but returns the position of the end (the last char) of
   the n-th word.
*/
char* ewordindex ( char *S, int n )
{
S = wordindex(S,n);
if (!S || !*S) return 0;
while (*S && !isspace(*S)) S++; // skip current word
return --S;
}


// Structure scans_range.

class t_scans_range {
  private:
    int n;          // number of elements
    int alloc;      // allocated elements
    int *from, *to; // note that to==from for flag 0
    int max_scan_nb;
  private:
    void need_one_more();
    int set_max_scan_nb();
  public:
    t_scans_range () { alloc=n=0; from = to = 0; max_scan_nb = -1; }
    void init ( const char *S );
    int contains ( int x );
    int get_n () { return n; }
    int get_from ( int i ) { return from[i]; }
    int get_to ( int i ) { return to[i]; }
    int get_max_scan_nb () { if (max_scan_nb<0) set_max_scan_nb(); return max_scan_nb; }
};


void t_scans_range::need_one_more ()
{
if (n!=alloc) return;
alloc += 8; // allocated array by chunks of 8
from = (int*) realloc( from, alloc*sizeof(int) );
to = (int*) realloc( to, alloc*sizeof(int) );
}


/* Parse string like "1,2,3,4:6,10,20:30,11" into the arrays.
*/
void t_scans_range::init ( const char *S )
{
if (!S || !*S) return;
char *tmp = strdup(S);
char *curr = tmp, *next;
int is_range;
// notes: S is begin of next, s is begin of current
n = 0;
while (curr) {
  need_one_more();
  is_range = 0;
  next = curr+1;
  while (*next && *next!=',') {
    if (*next==':') is_range = 1;
    next++;
  }
  if (*next) *next++ = 0;
    else next = 0;
  if (is_range) {
      sscanf(curr,"%i:%i",&from[n],&to[n]);
      if (to[n] < from[n]) { int t = to[n]; to[n] = from[n]; from[n] = t; }
      }
    else {
      sscanf(curr,"%i",&from[n]);
      to[n] = from[n];
      }
  n++;
  curr = next;
}
free(tmp);
max_scan_nb = -1;
}


/* Returns 1 if the struct contains number x, 0 otherwise.
*/
int t_scans_range::contains ( int x )
{
for (int i=0; i<n; i++)
  if (from[i]<=x && x<=to[i]) return 1;
return 0;
}


/* Sets the max_scan_nb, and returns it.
*/
int t_scans_range::set_max_scan_nb ()
{
for (int i=0; i<n; i++)
  if (to[i] > max_scan_nb) max_scan_nb = to[i];
return max_scan_nb;
}


#ifdef PSD

// Structure scans_range.

class t_psd_integr : public t_scans_range {
  public:
    t_psd_integr() { t_scans_range(); }
};

#endif /* PSD */


/* Usage / help.
*/
void Usage ()
{
cerr <<"## (c) Petr Mikulik  (http://www.sci.muni.cz/~mikulik/)  1994--2017\n\n";
cerr <<"Usage:\n";
cerr <<
"  unspec input_spec_file [ - | out_files_name [out_files_extension] ] [options]\n"
"\n"
"In general, this program is used to split the input spec file into a series of\n"
"output files each containing one scan only (unless option -1 or -2 is given or\n"
"stdout goes to stdout).\n"
"\n"
"If the output file is -, e.g. 'unspec xxx.spec - <opts>'), then the extracted\n"
"scans go all to screen (stdout). If the output is <out_name>, e.g. using\n"
"'unspec xxx.spec abc <opts>', then the output files are abc01, abc02, abc03,\n"
"etc. If the output file is specified by <out_name> <out_ext>, e.g. using\n"
"'unspec xxx.spec abc dat <opts>', then the output files are abc01.dat,\n"
"abc02.dat, abc03.dat, etc. Using option -3 changes the above to abc001, abc002,\n"
"...; abc001.dat, abc002.dat,..., option '-4' to abc0001, abc0002,...\n"
"With option -1 the output file name is <out_name><scan_from_nb><out_ext>,\n"
"while with option -2 the output file name does not contain the scan number,\n"
"thus being identical to the result of\n"
"'unspec input.spec - -1 >my_output.dat'.\n"
"\n";
cerr <<
"\nOptions:\n"
"  -s  <scan>        ... extract only this scan (default: all scans)\n"
"  -r  <from>  <to>  ... extract this range of scans (default: all scans)\n"
"  -f  <from>        ... extract all scans since this one\n"
"  -t  <to>          ... extract all scans until this one\n"
"  -R  <ranges>      ... extract given scans, syntax like: 1,5,8:12,30\n"
"  -T  ascan         ... only scans containing this string, e.g. ascan or a2scan;\n"
"                        use grep or text editor to look at lines with #S\n"
"  -1                ... all scans go to one file (i.e. append scans)\n"
"  -2                ... all scans go to one file, no scan number in its name\n"
"  -3                ... 3-digit scan number file indexing (default is 2)\n"
"  -4                ... 4-digit scan number file indexing (default is 2)\n"
"  -b  <n>           ... in each scan, write a blank line after each n-th point;\n"
"                        useful mainly when the scan has been produced by 'mesh'\n"
"  -S                ... spec file is not read after the last scan given via the\n"
"                        options -s, -r, -R and -t. This switch lets unspec read\n"
"                        the file until the end. Could be useful for a curious\n"
"                        case with two same scan numbers or scan indexing error\n"
"  -B                ... no blank line between scans (useful for writing out\n"
"                        matrices, use together '-B -[1|2] -#')\n"
"  -H                ... copy header of the input spec file to the output file\n"
"  -#                ... do not write out any spec comment (# lines) of scans\n"
"  -c  COMMENTSCHARS ... write these scan #comments (default "<<WRITE_COMMENTS<<")\n"
"                        Using lowercase letters strips the #comment name\n"
"                        (useful when importing data to a spreadsheet). Use\n"
"                        empty string, i.e. '-c \"\"), to output all comments\n"
"  -i                ... reindex scans (ignore scan number in #S, start from 1)\n"
#ifdef PSD
"  -Pn               ... don't write PSD/MCA spectra (lines with @A)\n"
"  -Py               ... write PSD/MCA spectra (lines with @A)\n"
"  -PMCADEV <n>      ... write the first PSD/MCA data found (-1, default) or\n"
"                        only those belonging to @MCADEV (0,1,...) as indicated\n"
"                        in scan headers (if multiple MCA devices connected)\n"
"  -Pr  <from>  <to> ... write out only ROI for PSD spectra (from>=1); writing\n"
"                        ROI can be suppressed by to<from\n"
"  -Pe  <n>          ... write out every n-th channel (default 1 = all channels)\n"
"  -Pa  <n>          ... average (smoothen) intensity over n PSD channels\n"
"  -PI  <ranges>     ... print columns of intensities integrated over the given\n"
"                        PSD channels, or output only specified PSD channels;\n"
"                        e.g. '-PI 100:400,255,250:260'\n"
"  -P@A...           ... write PSD spectra in columns on new lines starting @A\n"
"  -PR               ... write out PSD spectra in rows with the channel number\n"
"                        or with angular positions, see options '-P0' and '-P1'\n"
"  -P0  <x>          ... channel of the primary beam, i.e. angle=0; see '-PR'\n"
"  -P1  <x>          ... number of channels per 1 degree; see '-PR'\n"
"  -PM  <fileprefix> ... PSD/MCA spectra can be saved in the main spec file or\n"
"                        in separate files, see option 'to' in 'mcasetup' in\n"
"                        spec. Unspec supposes the former by default, while the\n"
"                        latter choice needs to know that <fileprefix>; then\n"
"                        spectra recorded during scans will be found and read.\n"
"                        Note that 'mcasave' files cannot be inserted into the\n"
"                        spec file, since they are not mentioned therein (hint:\n"
"                        see Solutions below what to do in that case)\n"
#endif
"  -v                ... be verbose (do not print text on screen)\n"
"  -w <columns>      ... which columns go to output (default: whole scan line)\n"
"                        <columns> specification can list columns, ranges of\n"
"                        columns, and motor or detector names (as found in #L\n"
"                        and #O lines of the spec file, ignoring letter case and\n"
"                        spaces). Since #L are column titles, unspec can reduce\n"
"                        number of columns on output, and extend it by values #P\n"
"                        of motors fixed during that scan.\n"
"                          Example:  -w 1,2:4,-1,-2,Theta,ome,scintil\n"
"                        selects column 1 to 5, the last and last but one, and\n"
"                        those of Theta and Ome motors and Scintil detector.\n"
"                        Further, the following special names are understood:\n"
"                          '#' print scan number\n"
"                          '@' print point number in the scan\n"
"                          '=' auto find + print all motors changed in the scan\n"
"                              (most standard scans and meshes are understood)\n"
"                          Example:  -w =,det\n"
"  -L                ... display software license\n"
"\n";
cerr <<
"Examples, where  'gesi'  is an experimental data file produced by spec:\n"
"  unspec gesi x dat -3      ... writes files x000.dat, x001.dat, x002.dat,...\n"
"  unspec gesi - -s 10 -#    ... writes scan 10 without comments to stdout\n"
"  unspec gesi p -r 20 22    ... makes p20, p21 and p22 with scans 20..22\n"
"  unspec gesi p -T \"ascan  th\" ... extract all ascans with th motor\n"
"  unspec gesi p q -s 6 -c l ... makes p06.q, its header describes data columns\n"
"  unspec gesi z -1 -w th,det -3 -f 90 ... makes single file z090 with all scans\n"
"                                from 90 and columns of motors 'th' and 'det'\n"
"  unspec gesi z.dat -2 -s 5 ... makes single file z.dat with scan number 5\n"
"  unspec gesi s dat -w =,det... write detector and the moving motor positions\n"
"\n"
"Solutions:\n"
" (*) Use the command \"grep #S mydata.spec\" to list all scans in a spec file.\n"
#ifdef PSD
" (*) Convert PSD/MCA files from 'mcaacq' / 'mcasave' commands to a single table\n"
"     1. These mca files gesi_mca/gesi_0sss_0nnn.mca were made during a scan\n"
"        written in .spec file:\n"
"        unspec gesi.spec mcascan dat -s 13 [OPTS] -PM gesi_mca/gesi\n"
"        where OPTS can be like  -w phi,delta,filter,mon  -PR\n"
"     2. Use only .mca files -- however, there is no information about other\n"
"        detectors like sec, mon, filter, so you cannot do proper normalization:\n"
"        rm -f out.dat\n"
"        for i in gesi_0*.mca; do unspec $i - [OPTS] >>out.dat; done\n"
"        where OPTS can be like  -w phi,delta -PR  (filter, mon not available)\n"
" (*) Make a single spec file with all mca spectra saved separately in mca/:\n"
"     unspec b1.spec b1-all.spec -2 -H -PM mca/b1 -P@A\n"
#endif
"\n"
"Examples for using unspec with gnuplot 4.0 (or later, see www.gnuplot.info):\n"
"For plotting a single scan of a running experiment, type a command like\n"
"  gnuplot> plot \"<unspec abc.spec - -v -# -s 10\" u 1:7 title \"scan 10\"\n"
"For plotting a map by spec's mesh (number of inner motor intervals 100), do\n"
"  gnuplot> !unspec abc.spec map.dat -2 -# -s 19 -w ome,tth,det -b 101\n"
"(mesh is scan nb 19) or when the map is measured by a series of scans 19..69\n"
"  gnuplot> !unspec abc.spec map.dat -2 -# -r 19 69 -w ome,tth,corr_det\n"
"and then\n"
"  gnuplot> set pm3d map; set autoscale fix; set log cb; splot 'map.dat'\n"
"\n"
"Note that you can print this documentation by commands (supposing bash shell):\n"
"  unspec 2>unspec.txt; unspec -L 2>>unspec.txt; a2ps -B -2 -o out.ps unspec.txt\n"
"\n";
exit(1);
}


/* Show the program license.
*/
void License ()
{
cerr <<
"\n"
"   Copyright (C) 1994--2017 Petr Mikulik\n\n"
"   This program is free software; you can redistribute it and/or modify\n"
"   it under the terms of the GNU General Public License as published by\n"
"   the Free Software Foundation; either version 2, or (at your option)\n"
"   any later version.\n"
"\n"
"   This program is distributed in the hope that it will be useful,\n"
"   but WITHOUT ANY WARRANTY; without even the implied warranty of\n"
"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"
"   GNU General Public License for more details.\n"
"\n"
"   You should have received a copy of the GNU General Public License\n"
"   along with this program; if not, write to the Free Software\n"
"   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.\n"
"\n\n"
"   The latest distribution of the program can be downloaded from the\n"
"   author's homepage http://www.sci.muni.cz/~mikulik/. Your comments or\n"
"   suggestions are welcome.\n"
"\n";
exit(0);
}


/* The following two variables are used for the processing of the input data
   (from a scan). Like in awk: NF=number of fields, dollar[i]=i-th column
   (dynamically allocated string), i=1..NF; dollar[0] is the whole line.
   dollar_title[i] is the description (title) of the column, read after #L.
*/
int NF = 0;
#ifndef PSD
const int MAX_NF = 256;
#else
const int MAX_NF = 4500;
#endif
  // there are no more than 30 columns in spec anyway unless using PSD / MCA
char* dollar[MAX_NF];
int DTNF = 0;
char* dollar_title[MAX_NF];

/* Parse columns of the input line to dollar array and NF.
*/
void fill_dollars ( char *S )
{
char *b, *p = S;
NF = 0;
if (dollar[0]) free(dollar[0]); // $0 contains the whole input line
dollar[0] = strdup(S);
for (NF=0; *p; ) {
  NF++;
  if (NF>=MAX_NF) {
    cerr <<"FATAL ERROR: please recompile unspec with larger MAX_NF.\n"; exit(1); }
  if (dollar[NF]) free(dollar[NF]);
  b = p; // begin of the current word
  while (*p && !isspace(*p)) p++; // skip until space or eoln
  dollar[NF] = (char*)malloc(p-b+1);
  memcpy(dollar[NF],b,p-b); dollar[NF][p-b] = 0;
  while (*p && isspace(*p)) p++; // skip until no space or eoln
  }
}


/* Parse column titles of the input line to dollar_title array.
   S starts with "#L " and keywords (which can contain 0 or 1 space) are
   separated by two spaces.
*/
void fill_dollar_titles ( char *S )
{
char *b, *p = S+3;
for (DTNF=0; *p; ) {
  DTNF++;
  if (DTNF>=MAX_NF) {
    cerr <<"FATAL ERROR: please recompile unspec with larger MAX_NF.\n"; exit(1); }
  if (dollar_title[DTNF]) free(dollar_title[DTNF]);
  b = p; // begin of the current word
  while (*p && (!isspace(*p) || !isspace(*(p+1)))) p++; // skip until 2 spaces or eoln
  dollar_title[DTNF] = (char*)malloc(p-b+1);
  char *q = b, *x = dollar_title[DTNF];
  for (; q != p; q++) if (!isspace(*q)) *x++ = *q; // tolower(*q);
  *x = 0;
  while (*p && isspace(*p)) p++; // skip until no space or eoln
  }
}


/* Names of motors parsed from #O0, #O1 etc. lines: supports up to
   #O99 (see MAX_OP_LINES), up to 99 (see MAX_OP_MOTORSonLINE) motors per line,
   and a motor name up to 16 chars.
   MotorValuesPn are values of the motor positions fixed at the current
   scan.
*/
char MotorNamesOn[MAX_OP_LINES][MAX_OP_MOTORSonLINE][16];
char MotorValuesPn[MAX_OP_LINES][MAX_OP_MOTORSonLINE][32];


/* Show known motor names (used for debugging).
*/
void ShowMotorNames ()
{
int n, k, u;
printf("\n*** List of all motors ***\n");
for (n=0; n<=99; n++) {
  u=0;
  for (k=0; k<=99; k++) {
    if (MotorNamesOn[n][k][0]) {
      printf("[%i,%i]: '%s'\t", n,k,MotorNamesOn[n][k]);
      u = 1;
      }
    }
    if (u) printf("\n");
  }
printf("\n");
}


/* Parse line like #O1, #O2, ... with motor names.
*/
void SetMotorNamesO_n ( int n, char *S )
{
memset( MotorNamesOn[n], 0, sizeof(MotorNamesOn[n]) );
char *b, *p = S;
if (S[0]=='#') p += 4; // #On not yet parsed out
for (int word=0; *p; word++) {
  while (isspace(*p)) p++;
  b = p; // begin of the current word
  while (*p && (!isspace(*p) || !isspace(*(p+1)))) p++; // skip until 2 spaces or eoln
  // MotorNamesOn[n][word] = (char*)malloc(p-b+1);
  // copy without spaces and lowercased
  char *a = MotorNamesOn[n][word];
  for (; b != p; b++) if (!isspace(*b)) *a++ = *b; // tolower(*b);
  *a = 0;
  while (*p && isspace(*p)) p++; // skip until no space or eoln
  }
}


/* The following number and arrays remember which columns were selected for
   output.
   Array of integers which_column contains positive number (number of
   columns), negative numbers (extraction from the end of the line) or
   0, in which case the symbolic name was given. Note that which_columns_str
   may not be only column in the current scan, but also a value of a motor
   fixed during this scan.
     When the column is given by a symbolic name (see Motor... variables)
   then which_column_namewhere knows where from to get the value:
   zero means non-existing motor/detector, positive number is the column of
   the dollar array, and negative addresses MotorValuesPn[n][k] by
   -(<n>*100+k).
*/
int nb_which_columns = 0;
int *which_column = 0;
char **which_column_str = 0;
int *which_column_namewhere = 0;
int which_column_namewhere_updatereq = 1; // setting which_column_namewhere required

/* Parse a string in the form:
	 "1,2,4,8:11,-1,-4"
   to which_column and which_column_str arrays and nb_which_columns.
*/
void parse_which_columns ( char *S )
{
if (which_column) delete [] which_column;
if (!S || !*S) { which_column=0; nb_which_columns=0; return; }
which_column = new int [ 2*MAX_NF ];
which_column_str = new char* [ 2*MAX_NF ];
which_column_namewhere = new int [ 2*MAX_NF ];
want_autoguess_names = 0;
autoguess_names = 0;
char c, *s = strdup(S);
char *p = s, *e;
int from, to;
while (*p) {
  e = strchr(p,','); // end of token
  if (e) *e = 0;
  if (strchr(p,':')) { // range of columns (two integers)
      if ( sscanf(p, "%i:%i%c", &from, &to, &c) != 2 ) {
	cout <<"ERROR: not valid range: " <<p <<'\n'; exit(1); }
      for ( ; from<=to; from++) {
	if (nb_which_columns >= 2*MAX_NF) {
	  cerr <<"ERROR: do you really want so many columns for output?\n"; exit(1); }
	which_column[ nb_which_columns++ ] = from;
	}
      }
    else { // single column specified either by a column number or motor/detector name
      if ( sscanf(p, "%i%c", &from,&c) != 1 ) {
	// It's a string = name, thus copy it lowercased without spaces.
	if (*p=='=' && !*(p+1)) want_autoguess_names = 1;
	which_column_str[nb_which_columns] = new char [strlen(p)+1];
	char *q = p, *b = which_column_str[nb_which_columns];
	for (; *q; q++) if (!isspace(*q)) *b++ = tolower(*q);
	*b = 0;
	from = 0; // indicator of using which_column_str
	}
      if (nb_which_columns >= 2*MAX_NF) {
	cerr <<"ERROR: do you really want so many columns for output?\n"; exit(1); }
      which_column[ nb_which_columns++ ] = from;
      }
  if (!e) break;
  p = e+1;
  }
}


/* Update which_column_namewhere: positive for dollar array, negative to
   MotorValuesPn[n][k].
*/
void Update_which_column_namewhere ()
{
int i, n, k;
// Firstly, exact match is required.
for (i=0; i<nb_which_columns; i++) {
  which_column_namewhere[i] = 0; // set it to 0, i.e. not given
  if (which_column[i]) continue; // non-zero means a column given by string
  int found = 0;
  for (k=1; !found && k<=DTNF; k++)
    if (!stricmp(which_column_str[i],dollar_title[k])) { found=1; break; };
  if (found) {
    which_column_namewhere[i] = k;
    continue;
    }
  found = 0;
  for (n=0; !found && n<MAX_OP_LINES; n++)
    for (k=0; !found && k<MAX_OP_MOTORSonLINE; k++)
      if (!stricmp(which_column_str[i],MotorNamesOn[n][k])) { found=1; break; };
  if (found) which_column_namewhere[i] = -(100*n+k);
    else which_column_namewhere[i] = 0;
  }
// Secondly, partial match is sufficient.
for (i=0; i<nb_which_columns; i++) {
  if (which_column_namewhere[i]) continue; // column already found
  if (which_column[i]) continue; // non-zero means a column given by string
  int found = 0;
  int len = strlen(which_column_str[i]);
  for (k=1; !found && k<=DTNF; k++)
    if (!strnicmp(which_column_str[i],dollar_title[k],len)) { found=1; break; };
  if (found) {
    which_column_namewhere[i] = k;
    continue;
    }
  found = 0;
  for (n=0; !found && n<MAX_OP_LINES; n++)
    for (k=0; !found && k<MAX_OP_MOTORSonLINE; k++)
      if (!strnicmp(which_column_str[i],MotorNamesOn[n][k],len)) { found=1; break; };
  if (found) which_column_namewhere[i] = -(100*n+k);
    else which_column_namewhere[i] = 0;
  }
// Finally, fill an extra field for autoguessed motors.
if (autoguess_names)
for (i=0; i<nb_which_columns; i++) {
  if (which_column_namewhere[i]) continue; // column already found
  if (which_column[i]) continue; // non-zero means a column given by string
  if (strcmp(which_column_str[i],"=")) continue; // continue if not "="
for (int a=0; a<autoguess_names; a++) {
  which_column_namewhere[i] = 30000; // flag to search autoguessed column elsewhere
  int found = 0;
  int len = strlen(autoguess_name[a]);
  for (k=1; !found && k<=DTNF; k++)
    if (!strnicmp(autoguess_name[a],dollar_title[k],len)) { found=1; break; };
  if (found) {
    autoguess_which_column_namewhere[a] = k;
    continue;
    }
  found = 0;
  for (n=0; !found && n<MAX_OP_LINES; n++)
    for (k=0; !found && k<MAX_OP_MOTORSonLINE; k++)
      if (!strnicmp(autoguess_name[a],MotorNamesOn[n][k],len)) { found=1; break; };
  if (found) autoguess_which_column_namewhere[a] = -(100*n+k);
    else autoguess_which_column_namewhere[a] = 0;
  }
} // for all autoguess_names
// Motor not found.
which_column_namewhere_updatereq = 0;
}


/* There was an automatic guess of motor names requested. Here the scan command
line after "#S <n>" is passed, e.g. "ascan th 0 7 10 0.1", and according to the
scan type motor names are extracted.
Scan types:
	ct 2					(2 words)
	{h|k|l}scan 0 7 10 2			(5 words)
	{a|d|step}scan Motor 0 7 10 2		(6 words)
	{a|d}2scan Motor1 0 7 Motor2 0 14 10 2		   (9 words)
	hklscan -0.1 0.1 -2.05 -1.95 0.04 0.06 100 0.5	   (9 words)
	[hkl]mesh Motor1 1.9 2.1 10 Motor2 -2.1 -1.9 10 2  (10 words)
	{a|d}3scan Motor1 0 7 Motor2 0 14 Motor3 0 21 10 2 (12 words)

*/
void parse_motors_from_scan_command ( char *S )
{
autoguess_names = 0;
while (*S && isspace(*S)) S++; // ignore leading spaces
switch (words(S)) {
  case 2: break;
  case 5: // {h|k|l}scan 0 7 10 2
	  autoguess_names = 1;
	  autoguess_name[0][0] = S[0];
	  autoguess_name[0][1] = 0;
	  break;
  case 6: // {a|d|step}scan Motor 0 7 10 2
	  autoguess_names = 1;
	  copy_word( autoguess_name[0], S, 2);
	  break;
  case 9: if (!strncmp(S,"hklscan",7)) {
	    // hklscan -0.1 0.1 -2.05 -1.95 0.04 0.06 100 0.5
	    autoguess_names = 3;
	    strcpy(autoguess_name[0],"h");
	    strcpy(autoguess_name[1],"k");
	    strcpy(autoguess_name[2],"l");
	    }
	  else {
	    // {a|d}2scan Motor1 0 7 Motor2 0 14 10 2
	    autoguess_names = 2;
	    copy_word( autoguess_name[0], S, 2);
	    copy_word( autoguess_name[1], S, 5);
	    }
	  break;
  case 10: // [hkl]mesh Motor1 1.9 2.1 10 Motor2 -2.1 -1.9 10 2
	  autoguess_names = 2;
	  copy_word( autoguess_name[0], S, 2);
	  copy_word( autoguess_name[1], S, 6);
	  break;
  case 12: // {a|d}3scan Motor1 0 7 Motor2 0 14 Motor3 0 21 10 2
	  autoguess_names = 3;
	  copy_word( autoguess_name[0], S, 2);
	  copy_word( autoguess_name[1], S, 5);
	  copy_word( autoguess_name[2], S, 8);
	  break;
  default: cerr <<" <no rule for motor autoguess>";
  }
}


/* Fill the OutCols string with motor values.
*/
void
fill_OutCols ( char *OutCols )
{
int i;
char tmp128[128];
/* Only selected columns are to be printed on output, thus
   parse the line to dollar array and NF. */
fill_dollars(S);
int printtab = 0;
for (i=0; i<nb_which_columns; i++)
  if (abs(which_column[i]) <= NF) {
    if (printtab)
      strcat( OutCols, "\t" );
    else
      printtab = -1;
    if (which_column[i] > 0) // column given by column nb
      strcat( OutCols, dollar[which_column[i]] );
    else
      if (which_column[i] < 0) // column wrt last column
	strcat( OutCols, dollar[NF+which_column[i]+1] );
      else { // column given by its title
	if (which_column_namewhere_updatereq)
	  Update_which_column_namewhere();
	int where = which_column_namewhere[i];
	int autoguess_where = 0;
	if (/* autoguess_names && */ where==30000) {
	  autoguess_where = 1;
	  where = autoguess_which_column_namewhere[ 0 ];
	}
	autoguess_next_motor: ;
	if (where > 0)
	  strcat( OutCols, dollar[ where ] );
	else
	  if (where < 0) {
	    int n = where / -100 - 1,
	    k = -where % 100;
	    if (*MotorValuesPn[n][k])
	      strcat( OutCols, MotorValuesPn[n][k] );
	    else // motor exists, but value is unknown in this scan
	      strcat( OutCols, "-1e38" );
	  } else
	      if (!strcmp(which_column_str[i],"#")) { // output scan number
		sprintf(tmp128,"%i",ScanNb);
		strcat( OutCols, tmp128 );
	      } else
		  if (!strcmp(which_column_str[i],"@")) { // output point number
		      sprintf(tmp128,"%i",PointNb );
		      strcat( OutCols, tmp128 );
		  } else // non-existing motor: which_column_namewhere[i]==0
		      printtab = 0;
	  if (autoguess_where) {
	    where = autoguess_which_column_namewhere[ autoguess_where++ ];
	    if (autoguess_where <= autoguess_names) {
	      if (printtab)
		strcat( OutCols, "\t" );
	      printtab = 1;
	      goto autoguess_next_motor;
	    }
	  }
      }
      if (printtab<0) printtab = 1;
  }
}


/* Print a table with all motor names and positions.
*/
void print_all_motor_positions ()
{
for (int m1=0; m1<MAX_OP_LINES; m1++) {
  cout <<"  " <<m1 <<":";
  for (int m2=0; m2<MAX_OP_LINES; m2++)
    cout <<"  " <<MotorNamesOn[m1][m2] <<": " <<MotorValuesPn[m1][m2];
  cout <<"\n";
}
}


#ifdef PSD
/*
Look whether a PSD line (MCA data) follows:
  (*) if "@A" follows: if write_psd then read the line and return 1,
      otherwise (i.e. !write_psd) skip it and return 0.
  (*) else "@A" does not follow: if !write_psd then return 0, otherwise
      (i.e. write_psd) see whether an <prefix>_scannb_ptnb.mca exists,
      and parse it to read "@A" line and return 0 or 1 accordingly.
If found, then data are read, but parsed only if they belong to the required
device (header line #@MCADEV 0).
Return values:
    0...PSD line not found,
    1...PSD line read and parsed,
    2...PSD line found and ignored (not parsed) because it was requested
	to read data from another device.
*/
int read_PSD_line ( ifstream& inp, int  write_it )
{
int len;
static char *S = 0;
const int sizeofS = 65536;
if (!S) S = new char [sizeofS];
if (inp.peek()=='@') { // There can be a PSD/MCA spectrum on the next line.
  unsigned long pos = inp.tellg();
  inp.getline(S,sizeofS-1);
  if (S[len=strlen(S)-1]=='\r') S[len] = 0;
  if (strncmp("@A",S,2)) {
    inp.seekg(pos); // Go back, that line is not a PSD spectrum.
    return 0; // mca data not available
    }
  // A PSD line follows.
  int stop = 0;
  int len;
  // psd.rawline[0] = 0; ... this must be done earlier because of potential multiple @MCADEV
  if (!write_it) { // Skip the @A .. \ .. \ ... part.
    while (1) {
      len = strlen(S) - 1;
      stop = S[len] != '\\'; // The last PSD line has no trailing backslash.
      if (stop) break;
      inp.getline(S,sizeofS-1);
      if (S[len=strlen(S)-1]=='\r') S[len] = 0;
      }
    return 2; // mca data found but ignored
    }
  // Get the @A line into psd.rawline.
  char *e, *c = S+2;
  while (1) {
    while (*c && isspace(*c)) c++; // strip early spaces
    len = strlen(c) - 1;
    stop = c[len] != '\\'; // The last PSD line has no trailing backslash.
    if (!stop) c[len] = 0;
    for ( e=c+len; e>c && isspace(*e); e-- ) ; // skip trailing spaces
    *++e = 0;
    if (psd.rawline[0]) { // separate lines by a single space
      if (c>S && *(c-1)==' ') c--;
	else strcat( psd.rawline, " " );
      }
    strcat( psd.rawline, c );
    if (stop) break;
    inp.getline(S,sizeofS-1);
    if (S[len=strlen(S)-1]=='\r') S[len] = 0;
    c = S;
    }
  // cout <<"   ... Got " <<psd.rawline <<"\n";
  return 1; // mca data found and parsed
  }
// psd.rawline[0] = 0; ... this is done earlier
if (!psd_mca_file_prefix) return 0; // mca data not found
// try to open / see whether file <prefix>_scannb_ptnb.mca exists
sprintf(S,"%s_%03i_%03i.mca", psd_mca_file_prefix,ScanNb,PointNb-1);
ifstream inpmca;
inpmca.open(S,ios::in);
if (!inpmca) return 0; // mca data not found
cerr <<" (" <<S <<")";
while (inpmca) {
  inpmca >> S;
  if (!strcmp(S,"@A")) break;
  inpmca.ignore(65530,'\n');
  }
if (!inpmca) return 0; // mca data not found
if (!write_it) return 2; // mca data found but no need to read them
// Get the @A line into psd.rawline (same code as a little bit above).
char notfirst=0, *e, *c;
int stop = 0;
inpmca.getline(S,sizeofS-1);
if (S[len=strlen(S)-1]=='\r') S[len] = 0;
while (1) {
  c = S;
  while (*c && isspace(*c)) c++; // strip early spaces
  len = strlen(c) - 1;
  stop = c[len] != '\\'; // The last PSD line has no trailing backslash.
  if (!stop) c[len] = 0;
  for ( e=c+len; e>c && isspace(*e); e-- ) ; // skip trailing spaces
  *++e = 0;
  if (notfirst) {
    if (c>S && *(c-1)==' ') c--;
      else strcat( psd.rawline, " " );
    }
  strcat( psd.rawline, c );
  if (stop || !inpmca) break;
  inpmca.getline(S,sizeofS-1);
  if (S[len=strlen(S)-1]=='\r') S[len] = 0;
  notfirst = 1;
  }
return 1; // mca data found and parsed
}


// Read or ignored all sequent MCA spectra.
// Thus call read_PSD_line(inp,write_it) until it returns zero.
int read_PSD_line ( ifstream& inp )
{
static int counter=0;
int local = 0;
int write_if_found = write_psd && (psd_which_mcadev == -1);
int psd_lines_found = -1;
int data_written = 0;
int res = 1;
psd.rawline[0] = 0; // empty MCA data (zero length string)
while (res) {
    counter++; local++;
    if (psd.curr_which_mcadev >= 0) // it was explicitly set on command line which mcadev to write
      write_if_found = write_psd && (psd_lines_found+1 == psd.curr_which_mcadev);
    res = read_PSD_line(inp, write_if_found);
    if (res > 0 && write_if_found) {
	write_if_found = 0;
    	data_written = 1;
    }
    if (res > 0)
	psd_lines_found++;
  }
// cout <<"  rawline: " <<psd.rawline <<endl;
return data_written;
}

#endif // PSD


/*
 * Main.
 */

int main ( int argc, char *argv[] )
{
// Write header only in non-verbose mode.
int verbose = 0;
int optarg = 2;
while (optarg<argc && !verbose)
   if (!strcmp(argv[optarg++],"-v")) verbose = 1;

if (!verbose || argc<=2) {
  cerr <<"\n## Program unspec (version " <<VERSION <<") ##\n";
  }

if (argc==2 && !strcmp(argv[1],"-L")) License();
if (argc<=2) Usage();

// Local variables:
int MaxScanNb_last = -65535;
int pt_in_scan = -1; // index going from blank_after_pts to 0;
int last_scan = 0; // last requested scan (from options -s, -r, -t, -R)
int stop_last_scan = 1; // enable stopping after reading the file after max. requested scans
t_scans_range scans_range;
int skip_scan = 0; // 1 if the current scan is to be skipped
#ifdef PSD
t_psd_integr psd_integr;
#endif

// Initialization.
memset( &dollar[0], 0, sizeof(dollar) );
memset( &dollar_title[0], 0, sizeof(dollar_title) );

optarg = 2;
if (argv[optarg][0]=='-') {
  if (!strcmp(argv[optarg],"-L")) License();
  if (argv[optarg][1]) Usage();
  }
char *FileNameBase = argv[optarg++];
char *FileNameExt = 0;
if (!strcmp(FileNameBase,"-")) { // write output to stdout
      out2stdout = 1;
      FileNameExt = 0;
      }
  else if (optarg<argc && argv[optarg][0] != '-') // file name extension given, not an option
	    FileNameExt = argv[optarg++];

// Parse command line options.
for (; optarg<argc; optarg++) {
  if (!strcmp(argv[optarg],"-s")) {
    ArgsRequired(1,optarg, argc,argv);
    scan_to = scan_from = (int) atol(argv[++optarg]);
    last_scan = scan_to;
    }
  else
  if (!strcmp(argv[optarg],"-f")) {
    ArgsRequired(1,optarg, argc,argv);
    scan_from = (int) atol(argv[++optarg]);
    last_scan = 0;
    }
  else
  if (!strcmp(argv[optarg],"-t")) {
    ArgsRequired(1,optarg, argc,argv);
    scan_to = (int) atol(argv[++optarg]);
    last_scan = scan_to;
    }
  else
  if (!strcmp(argv[optarg],"-r")) {
    ArgsRequired(2,optarg, argc,argv);
    scan_from = (int) atol(argv[++optarg]);
    scan_to = (int) atol(argv[++optarg]);
    if (scan_from>scan_to) {
      cerr <<"ERROR: Scan from > scan to, please revoke your idea.\n";
      exit(1); }
    last_scan = scan_to;
    }
  else
  if (!strcmp(argv[optarg],"-R")) {
    ArgsRequired(1,optarg, argc,argv);
    scans_range.init( argv[++optarg] );
    use_scans_range = 1;
    last_scan = scans_range.get_max_scan_nb();
    }
  else
  if (!strcmp(argv[optarg],"-T")) {
    if (scan_substring_filter) free(scan_substring_filter);
    scan_substring_filter = strdup(argv[++optarg]);
    }
  else
  if (!strcmp(argv[optarg],"-c")) {
    ArgsRequired(1,optarg, argc,argv);
    WriteComments = argv[++optarg]; }
  else
  if (!strcmp(argv[optarg],"-#"))
    write_grid = 0;
  else
  if (!strcmp(argv[optarg],"-H"))
    write_header = 1;
  else
  if (!strcmp(argv[optarg],"-1")) {
    single_file = 1;
    no_scannb_in_outfile = 0;
    }
  else
  if (!strcmp(argv[optarg],"-2")) {
    single_file = 1;
    no_scannb_in_outfile = 1;
    }
  else
  if (!strcmp(argv[optarg],"-3"))
    scan_digits = 3;
  else
  if (!strcmp(argv[optarg],"-4"))
    scan_digits = 4;
  else
  if (!strcmp(argv[optarg],"-i"))
    reindex_scans = 1;
  else
  if (!strcmp(argv[optarg],"-v"))
    verbose = 1;
  else
  if (!strcmp(argv[optarg],"-w")) {
    ArgsRequired(1,optarg, argc,argv);
    parse_which_columns( argv[++optarg] );
    }
  else
  if (!strcmp(argv[optarg],"-L"))
    License();
  else
  if (!strcmp(argv[optarg],"-b")) {
    ArgsRequired(1,optarg, argc,argv);
    blank_after_pts = (int) atol(argv[++optarg]);
    if (blank_after_pts < 0) blank_after_pts = 0;
    }
  else
  if (!strcmp(argv[optarg],"-B"))
    blank_between_scans = 0;
  else
#ifdef PSD
  if (!strcmp(argv[optarg],"-Pn")) {
    write_psd = 0;
    }
  else
  if (!strcmp(argv[optarg],"-Py")) {
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-PMCADEV")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_which_mcadev = atol(argv[++optarg]);
    if (psd_which_mcadev < -1) psd_which_mcadev = -1; // reset to default
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-Pr")) {
    ArgsRequired(2,optarg, argc,argv);
    psd_roi_from = (int) atol(argv[++optarg]);
    psd_roi_to = (int) atol(argv[++optarg]);
#if 0
    if (psd_roi_from>psd_roi_to) {
      cerr <<"ERROR: PSD ROI from > to, please revoke your idea.\n";
      exit(1); }
#endif
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-PR")) {
    psd_write_rows = 1;
    psd_write_aA = 0;
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-P@A")) {
    psd_write_aA = 1;
    psd_write_rows = 0;
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-Pe")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_every = (int) atol(argv[++optarg]);
    if (psd_every < 1) psd_every = 1;
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-Pa")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_average = (int) atol(argv[++optarg]);
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-P0")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_centre = atof(argv[++optarg]);
    psd_ang_scale = 1;
    }
  else
  if (!strcmp(argv[optarg],"-P1")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_chans_1deg = atof(argv[++optarg]);
    psd_ang_scale = 1;
    }
  else
  if (!strcmp(argv[optarg],"-PM")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_mca_file_prefix = argv[++optarg];
    write_psd = 1;
    }
  else
  if (!strcmp(argv[optarg],"-PI")) {
    ArgsRequired(1,optarg, argc,argv);
    psd_integr.init( argv[++optarg] );
    write_psd = 1;
    }
  else
#endif
  if (!strcmp(argv[optarg],"-S"))
    stop_last_scan = 0;
  else {
    cerr <<"ERROR: Wrong option " <<argv[optarg] <<endl;
    exit(1);
    }
  }

char *WriteComments_Upper = strdup(WriteComments);
for (char *p = WriteComments_Upper; *p; p++) *p = toupper(*p);

// Input and output files declaration.
ifstream inp;
if (!verbose) {
  cerr <<"Input file " <<argv[1] <<"; extracting ";
  if (use_scans_range) cerr <<"selected scans";
  else if (scan_from==SCAN_LOWER && scan_to==SCAN_UPPER) cerr <<"all scans";
  else if (scan_from==SCAN_LOWER) cerr <<"scans upto " <<scan_to;
  else if (scan_to==scan_from) cerr <<"only scan number " <<scan_from;
  else if (scan_to==SCAN_UPPER) cerr <<"scans from " <<scan_from;
	 else cerr <<"scans from " <<scan_from <<" to " <<scan_to;
  }
  if (scan_substring_filter) cerr <<" containing " <<scan_substring_filter;
  cerr <<endl;

#if defined(__GNUC__) && (__GNUC__ < 3)
  #define ios_binary ios::bin
#else
  #define ios_binary ios::binary
#endif
inp.open( argv[1], ios::in|ios_binary );
if (!inp) { cerr <<"\nERROR: Cannot open file " <<argv[1] <<".\n"; return 1; };

int i;
char OutFile[127];
char *header = 0;
#ifndef PSD
char OutCols[4096];
#else
char OutCols[2*4096];
#endif
 /* OutCols server for storing columns of motors/counters before going to
    output: they are used once for usual data and several times if output
    PSD spectra in rows.
      More space is allocate for PSD, since integrated intensity goes there
    too (intensities from the option '-PI'.
  */
#ifdef PSD
int psd_curr_mcadev = 0; // enumerates appearance of (sequence of) @MCADEV in the current scan
#endif

/* Copy header of the input file to the output file: copy everything until
the first blank line of the word "#S".
*/
if (write_header) {
  int len;
  unsigned long pos;
  unsigned int header_alloc = 1024;
  header = (char*) malloc(header_alloc);
  header[0] = 0;
  do {
    pos = inp.tellg();
    inp.getline(OutCols,sizeof(OutCols)-1);
    if (!strncmp("#S ",OutCols,3)) {
      inp.seekg(pos);
      break; }
    if (OutCols[len=strlen(OutCols)-1]=='\r') OutCols[len] = 0;
    if (!header || header_alloc < strlen(header)+len+1) {
      header_alloc += 256;
      header = (char*) realloc(header,header_alloc);
      }
    strcat(OutCols,"\n");
    strcat(header,OutCols);
    }
  while (OutCols[0]);
  }

if (out2stdout) {
  Out = (ostream*)(&cout);
  single_file = 3;
  if (header && header[0]) out <<header;
  }

/* Here begins the main loop going through all scans in the spec file.
*/
BeginLoop: ;
while (!(S[0]=='#' && S[1]=='S' && S[2]==0)) {
  inp >> S;
  if (!inp) goto TotalEnd;
  if (S[1]=='O' && S[2]) {
      // Line #O1, #O2, ... with motor names. Thus remember them.
      int n = S[2]-'0'; // #O0..#O9
      if (isdigit(S[3])) n = n*10 + (S[3]-'0'); // #O10..#O99
      inp.getline(S,sizeof(S)-1);
      SetMotorNamesO_n( n, S );
      }
    else
      if (!(S[0]=='#' && S[1]=='S')) inp.ignore(32535,'\n');
  }
if (MaxScanNb < ScanNb) MaxScanNb = ScanNb;
inp >> DataFileScanNb;
ScanNb = DataFileScanNb;
if (stop_last_scan && last_scan > 0 && ScanNb>last_scan) goto TotalEnd; // finish: don't read further scans
if (blank_after_pts) pt_in_scan = blank_after_pts+1;
inp.getline(S,sizeof(S)-1);
skip_scan = ((!use_scans_range && (ScanNb<scan_from || ScanNb>scan_to)) || // ignore this scan because of range
	     (use_scans_range && !scans_range.contains(ScanNb)));
if (!skip_scan) // ignore this scan because of substring
    if (scan_substring_filter && !strstr(S, scan_substring_filter)) skip_scan = 1;

  if (skip_scan) {
      // cerr <<"ignoring " <<ScanNb <<"  S=|" <<S <<"|\n";
      }
  else {                                  // read and write this scan
      if (single_file<2) {
	if (reindex_scans)
	    ScanNb = new_scan_index++;
	  else
	    if (ScanNb <= MaxScanNb && MaxScanNb > MaxScanNb_last) {
	      cerr <<"\n\n*** SCAN NUMBERING ERROR IN SPEC FILE: scan " <<ScanNb
		   <<" is before " <<MaxScanNb <<" (try option -i) ***\n\n";
	      MaxScanNb_last = MaxScanNb;
	      }
	// open the output file
	if (no_scannb_in_outfile) {
	    if (FileNameExt)
		sprintf(OutFile, "%s%s", FileNameBase, FileNameExt);
	      else
		sprintf(OutFile, "%s", FileNameBase);
	    }
	else {
	    const char FormatDigits2[] = "%s%02i";
	    const char FormatDigits3[] = "%s%03i";
	    const char FormatDigits4[] = "%s%04i";
	    char *Format1 = (char*)&FormatDigits2[0];
	    switch (scan_digits) {
		case 3: Format1 = (char*)&FormatDigits3[0];
			break;
		case 4: Format1 = (char*)&FormatDigits4[0];
			break;
	    }
	    sprintf(OutFile, Format1, FileNameBase, ScanNb);
	    if (FileNameExt) {
		strcat(OutFile, ".");
		strcat(OutFile, FileNameExt);
	        }
	    }

	if (!verbose) {
	  #if 1
	  cerr <<"[" <<DataFileScanNb <<" => " <<OutFile <<"] ";
	  #else // the old method
	  cerr <<"  Scan ";
	  if (DataFileScanNb<10) cerr <<" ";
	  cerr <<DataFileScanNb <<" goes to file " <<OutFile <<endl;
	  #endif
	  }
	if (!fOut) fOut = new ofstream;
	fOut->open(OutFile, ios_binary);
	if (!fOut->is_open()) { cerr <<"\nERROR: Cannot open file " <<OutFile <<".\n"; return 1; };
	Out = (ostream*)fOut;
	if (single_file) single_file = 2; // do not open the output file anymore
	if (header && header[0]) out <<header;
	}
      else
	if (!out2stdout && !verbose)
	  #if 1
	  cerr <<"[" <<OutFile <<" += " <<ScanNb <<"] ";
	  #else
	  cerr <<"  Scan " <<ScanNb <<"   added to   " <<OutFile <<endl;
	  #endif

    if (S[i=strlen(S)-1]=='\r') S[i] = 0;
    if (single_file > 1) {
      static int notfirst = 0;
      if (notfirst && blank_between_scans) out <<'\n';
      notfirst = 1; }
    if (write_grid && (!WriteComments[0] || strchr(WriteComments,'S')))
      out <<"#S " <<ScanNb <<' ' <<S <<'\n';
    if (want_autoguess_names) parse_motors_from_scan_command(S);
    which_column_namewhere_updatereq = 1; // must update motor values according to their names
    PointNb = 0;
#ifdef PSD
    // initialize PSD variables at the beginning of each scan
    psd.curr_which_mcadev = (psd_which_mcadev == -1 ? -1 : -2); // -1: default, -2: need to find @MCADEV
    psd_curr_mcadev = 0; // enumerates appearance of (sequence of) @MCADEV in the current scan
#endif
    do { // Main loop over lines in a scan.
      inp.getline(S,sizeof(S)-1);
      if (S[i=strlen(S)-1]=='\r') S[i] = 0;
#ifdef DEBUG
      cout <<"\n********************* " <<__FUNCTION__ <<"(): " <<__LINE__
	   <<" **********************\n********************* " <<S <<" ***\n";
#endif
      if (!*S) continue; // blank line
//    if (!strcmp(S,"in miread")) continue; // remove some crazy stuff at ID19 (appeared long time ago)
      if (S[0]!='#') { // Read and parse the data line.
	  PointNb++;
	  if (blank_after_pts) {
	    if (--pt_in_scan==0) {
		pt_in_scan = blank_after_pts;
		out <<"\n";
		}
	    }
	  *OutCols = 0; // Output of motors/counters: printed at the end or before each PSD channel.
	  if (!nb_which_columns) { // whole line is printed
	      strcat(OutCols,S);
	      }
	    else {
	      /* Only selected columns are to be printed on output, thus
		 parse the line to dollar array and NF.
	      */
#if 1
	      /* The code in routine fill_OutCols(OutCols) is identical to that
	         below as the routine was added on 25. 4. 2004 for filling of
		 OutCols for single .mca files.
	      */
	      fill_OutCols(OutCols);
#else
	      /* ... and thus (see the above remark), this code can be removed
	         any later.
	      */
	      fill_dollars(S);
	      int printtab = 0;
	      char tmp128[128];
	      for (i=0; i<nb_which_columns; i++)
		if (abs(which_column[i]) <= NF) {
		  if (printtab) strcat( OutCols, "\t" ); else printtab = -1;
		  if (which_column[i] > 0) // column given by column nb
		    strcat( OutCols, dollar[which_column[i]] );
		  else if (which_column[i] < 0) // column wrt last column
		    strcat( OutCols, dollar[NF+which_column[i]+1] );
		  else { // column given by its title
		    if (which_column_namewhere_updatereq) Update_which_column_namewhere();
		    int where = which_column_namewhere[i];
		    int autoguess_where = 0;
		    if (/* autoguess_names && */ where==30000) {
		      autoguess_where = 1;
		      where = autoguess_which_column_namewhere[ 0 ];
		      }
		    autoguess_next_motor: ;
		    if (where > 0)
			strcat( OutCols, dollar[ where ] );
		    else if (where < 0) {
		      int n = where / -100 - 1,
			  k = -where % 100;
		      if (*MotorValuesPn[n][k]) strcat( OutCols, MotorValuesPn[n][k] );
			else strcat( OutCols, "-1e38" ); // motor exists, but value is unknown in this scan
		      }
		    else if (!strcmp(which_column_str[i],"#")) { // output scan number
		      sprintf(tmp128,"%i",ScanNb);
		      strcat( OutCols, tmp128 );
		      }
		    else if (!strcmp(which_column_str[i],"@")) { // output point number
		      sprintf(tmp128,"%i",PointNb );
		      strcat( OutCols, tmp128 );
		      }
		    else // non-existing motor: which_column_namewhere[i]==0
		      { printtab = 0; }; // out <<"1e38";
		    if (autoguess_where) {
		      where = autoguess_which_column_namewhere[ autoguess_where++ ];
		      if (autoguess_where <= autoguess_names) {
			if (printtab) strcat( OutCols, "\t" ); printtab = 1;
			goto autoguess_next_motor;
			}
		      }
		    }
		  if (printtab<0) printtab = 1;
		  }
#endif
	      }
#ifndef PSD
	      out <<OutCols <<'\n';
#else
	      if (!write_psd) out <<OutCols <<'\n';
		// else OutCols is now properly set for use in PSD
#endif
	  }
	else { // Read and parse a comment line.
	  OutCols[0] = 0; // No OutCols (forget the old value when reading a new comment line).
	  if (S[1]=='O' && S[2]) {
	    /* Line #O1, #O2, ... with motor names. Thus remember them.
	       Well, these #O happen only in the preamble, but no problem
	       to keep the code here: maybe spec file was restarted.
	    */
	    int n = S[2]-'0'; // #O0..#O9
	    if (isdigit(S[3])) n = n*10 + (S[3]-'0'); // #O10..#O99
	    SetMotorNamesO_n( n, S );
	    }
	  if (S[1]=='P' && S[2]) {
	    // Set motor values.
	    int n = S[2]-'0'; // #P0..#P9
	    int offs = 3;
	    if (isdigit(S[3])) { // #P10..#P99
		n = n*10 + (S[3]-'0');
		offs++;
	    }
	    memset( MotorValuesPn[n], 0, sizeof(MotorValuesPn[n]) );
	    int word = 0;
	    for (char *p, *b = S+offs; *b; word++) { // for all words
	      while (*b && isspace(*b)) b++;
	      if (!*b) break;
	      p = &MotorValuesPn[n][word][0];
	      while (*b && !isspace(*b)) *p++ = *b++;
	      *p = 0;
	      }
	    }
	  if (nb_which_columns && S[1]=='L') {
	    /* Line with column description (motor and detector names) has to be
	       updated when only some of the columns are selected for output.
	    */
	    fill_dollar_titles(S);
	    S[3] = 0;
	    int print2sp = 0; // flag whether to print 2 spaces (columns separation)
	    for (i=0; i<nb_which_columns; i++)
	      if (abs(which_column[i]) <= DTNF) {
		  #if 1 // separate fields by two spaces, as spec does
		    if (print2sp) strcat(S,"  "); else print2sp = -1;
		  #else // separate fields by tabulator
		    if (print2sp) strcat(S,"\t"); print2sp = 1;
		  #endif
		  if (which_column[i] > 0)
		    strcat(S,dollar_title[which_column[i]]);
		  else if (which_column[i] < 0)
		    strcat(S,dollar_title[DTNF+which_column[i]+1]);
		  else {
		    if (which_column_namewhere_updatereq) Update_which_column_namewhere();
		    int where = which_column_namewhere[i];
		    int autoguess_where = 0;
		    if (/* autoguess_names && */ where==30000) {
		      autoguess_where = 1;
		      where = autoguess_which_column_namewhere[ 0 ];
		      }
		    L_autoguess_next_motor: ;
		    if (where > 0) // such motor/detector exists
			strcat(S,dollar_title[where]);
		    else if (where < 0) { // such motor/detector exists
			int n = -where / 100 - 1,
			    k = -where % 100;
			strcat(S,MotorNamesOn[n][k]);
		      }
		    else if (!strcmp(which_column_str[i],"#")) { // output scan number
		      strcat(S,"ScanNb");
		      }
		    else if (!strcmp(which_column_str[i],"@")) { // output point number
		      strcat(S,"PointNb");
		      }
		    else print2sp = 0;
		    if (autoguess_where) {
		      where = autoguess_which_column_namewhere[ autoguess_where++ ];
		      if (autoguess_where <= autoguess_names) {
			if (print2sp) strcat(S,"  "); print2sp = 1;
			goto L_autoguess_next_motor;
			}
		      }
		    }
		  if (print2sp<0) print2sp = 1;
		  }
#ifdef PSD
	    if (write_psd) {
	      for (i=0; i<psd_integr.get_n(); i++) {
		if (psd_integr.get_from(i)==psd_integr.get_to(i))
		    sprintf(S+strlen(S), "  MCA%i", psd_integr.get_from(i));
		  else
		    sprintf(S+strlen(S), "  MCA%i:%i", psd_integr.get_from(i),psd_integr.get_to(i));
		}
	      if (!(psd_roi_from > 0 && psd_roi_from > psd_roi_to)) { // no roi written
		 if (psd_write_rows) {
		    #if 0
		    strcat(S, psd_ang_scale ? "  MCAang" : "  MCAchan");
		    #else
		    strcat(S, "  MCAchan");
		    #endif
		    strcat(S, "  MCAint");
		    }
		  else {
		    if (!psd_write_aA) strcat(S, "  MCAroi");
		    }
		 }
	       }
#endif
	    }
	  if (write_grid) { // write allowed comments from the scan header
	    if (!WriteComments[0] || strchr(WriteComments,S[1])) out <<S <<'\n';
	      else
	    if (strchr(WriteComments_Upper,S[1])) out <<stripword1(S) <<'\n';
	    }
#ifdef DEBUG
	  cout <<"--- " <<__LINE__ <<" parsing of a comment line done.\n";
#if 0
	  cout <<"--- and current motor positions are:\n";
	  cout <<"  " <<S <<endl;
	  cout <<"  nb_which_columns=" <<nb_which_columns <<" OutCols=|" <<OutCols <<"|\n";
	  print_all_motor_positions();
#endif
#endif
#ifdef PSD
	  if (write_psd && psd_which_mcadev >= 0 && !strncmp(S+1, "@MCADEV ", 8)) { // need to parse #@MCADEV 1
	    if (strlen(S) >= 10) {
		int n = S[9]-'0'; // #@MCADEV 0..9
		if (isdigit(S[10])) n = n*10 + (S[10]-'0'); // #@MCADEV 10..99
		if (n == psd_which_mcadev)
		    psd.curr_which_mcadev = psd_curr_mcadev; // enumerates @MCADEV in the current scan
		psd_curr_mcadev++;
	    }

	  }
#endif
	}
#ifdef PSD
      /* There was a motor/counters line (then OutCols is set to a non-zero
      string), or there was a comment. If there was not comment, then
      look whether there follows a PSD line:
      if yes and !write_psd, then skip it; if !follows then try to read it
      from an .mca file. If successfully read, then output it as requested.
      If there was a comment, then a PSD line can immediately follow.
      */
      if ( (S[0]=='#' && inp.peek()=='@' && read_PSD_line(inp)) ||
	   (S[0]!='#' && read_PSD_line(inp)) ) {
	if (nb_which_columns > 0 && !*OutCols) {
	  /* Obviously reading not a .spec file, but a single .mca file. In the
	     former case, OutCols with motor values are filled when reading the
	     data line. In the latter, there are no data lines, only once @A line,
	     so we have to let fill OutCols here if it is (still) empty.
	     Note: .mca file contains only motor values, but not detectors.
	  */
	  fill_OutCols(OutCols);
	}

	if (psd_average > 1 || psd_integr.get_n() > 0) // rawline => (a)intensities
	  psd.do_intensities();
	for (i=0; i<psd_integr.get_n(); i++) { // intensity integrated over PSD channels
	  int from, to;
	  from = psd_integr.get_from(i) - 1; // count from and to 0..
	  to = psd_integr.get_to(i) - 1;
	  if (from<0) from = 0;
	  if (to>=psd.channels) to = psd.channels-1;
	  if (from > to)
	      sprintf(OutCols+strlen(OutCols),"\t-1");
	    else if (psd_average <= 1) {
	      int res = 0, k = from;
	      while (k<=to) res += psd.intens[k++];
	      sprintf(OutCols+strlen(OutCols),"\t%i",res);
	      }
	    else {
	      double res = 0;
	      int k = from;
	      while (k<=to) res += psd.aintens[k++];
	      sprintf(OutCols+strlen(OutCols),"\t%g",res);
	      }
	  }
	// now write the roi
	if (psd_roi_from > 0 && psd_roi_from > psd_roi_to) { // no roi written
	  if (*OutCols) out <<OutCols;
	  }
	else
	if (psd_roi_from <= 0 && psd_every <= 1 && !psd_write_rows && psd_average <= 1) {
	  // No ROI nor every nor averaging given, thus just copy the raw line.
	  if (*OutCols) out <<OutCols << (psd_write_aA ? "" : "\t");
	  if (psd_write_aA) out <<"\n@A ";
	  out <<psd.rawline;
	  }
	else
	if (psd_average <= 1) {
	  // Output in words (psd_every-th word from psd_from to psd_to).
	  // Note: psd_roi_from-th word counts from 1.
	  char *from = wordindex(psd.rawline,psd_roi_from > 0 ? psd_roi_from : 1);
	  if (from) {
	    if (psd_every <= 1 && !psd_write_rows) {
		// Write the whole substring from..to.
		if (psd_roi_from > 0) {
		  char *to = ewordindex(from,psd_roi_to-psd_roi_from+1); // end of word
		  if (to && *to) *(to+1) = 0;
		  }
		if (!psd_write_aA) out <<OutCols <<'\t' <<from;
		  else out <<OutCols <<"\n@A " <<from;
		}
	      else { // from:every:to
		char etmp, *e, *curr = from;
		int last = (psd_roi_from > 0) ? psd_roi_to : psd.channels+32000;
		int chan = (psd_roi_from > 0) ? psd_roi_from : 1;
		if (!psd_write_aA) out <<OutCols <<'\t';
		  else out <<OutCols <<"\n@A ";
		do {
		  e = ewordindex(curr,1); // firstly, print current word
		  etmp = *++e;
		  *e = 0;
		  if (curr!=from) {
		    if (!psd_write_rows) out <<' ';
		      else out <<OutCols <<'\t';
		    }
		  if (psd_write_rows) {
		    if (!psd_ang_scale) out <<chan <<'\t';
		      else out <<((chan-psd_centre)/psd_chans_1deg) <<'\t';
		    }
		  out <<curr;
		  if (psd_write_rows) {
		    out <<'\n';
		    chan += psd_every; }
		  *e = etmp;
		  curr = wordindex(e,psd_every); // now find next word
		  }
		while (curr && chan<=last);
		}
	    }
	  }
	else {
	  // Output the averaged intensities as numbers.
	  if (psd_average < 2) cerr <<"PSD_AVG>0?! FATAL.\n";
	  // Start from user-specified channel or from the 'psd_every' region
	  if (!psd_write_rows) {
	    if (!psd_write_aA) out <<OutCols <<'\t';
	      else out <<OutCols <<"\n@A ";
	    }
	  int chan = (psd_roi_from > 0) ? psd_roi_from-1 : (psd_every-1) / 2;
	  int to = (psd_roi_from > 0) ? psd_roi_to-1 : psd.channels-1;
	  while (chan <= to) {
	    if (psd_write_rows) {
		out <<OutCols <<'\t';
		if (!psd_ang_scale) out <<'\t' <<chan+1 <<'\t';
		  else out <<((chan+1-psd_centre)/psd_chans_1deg) <<'\t';
		}
	      else { if (chan!=0) out <<' '; }
	    if (psd_average > 1) out <<psd.aintens[chan];
	      else out <<psd.intens[chan];
	    if (psd_write_rows) out <<'\n';
	    chan += psd_every;
	    }
	  }
        out <<'\n';
	} // end of outputting a PSD line
      else { // PSD line was not read
	if (OutCols[0] && write_psd) out <<OutCols <<'\n';
	}
#endif
      }
    while (S[0] != 0);

    if (!single_file) fOut->close();
    }
goto BeginLoop;

TotalEnd:;
if (single_file==2) fOut->close();

if (!verbose) {
  cerr <<endl;
  if (scan_from==SCAN_LOWER) cerr <<"Uff, extraction completed. Have fun!\n"; else
  if (scan_from==scan_to) cerr <<"No problem to extract a single scan. Enjoy it!\n"; else
    cerr <<"Grrr, enjoy the extraction.\n";
  }
return 0;
}

// eof unspec.cpp
