#include "gfit.h"

#define MAXSTR		256

using namespace std;

/* 	Data File Format


	1. lines begining with # will be regarded as comments
	2. RESIDUE ### [0|1]
	3. 0.05 1.22 [0.12]
	4. &

	e.g.
	
	RESIDUE 121 1
	0.1	0.02	0.001
	0.2	0.03	0.002
	0.4	0.05	0.001
	0.6	0.12	0.002
	0.75	0.18	0.002
	0.9	0.20	0.001
	1.0	0.21	0.001
	1.2	0.22	0.001
	1.3	0.22	0.002
	&
	RESIDUE 122 0
	...
*/


/* Check DATAFILE */
void chk_data (const char *file, int &NR, int &NX)
{
char line[MAXSTR];
vector <string> c;//parsed columns
set <int> R;
int l,resid,flag,npts,nfit;

NR=NX=npts=nfit=0;

ifstream data_file(file);
if (data_file.fail()) {
	stringstream error;
	error <<"cannot open data file, "<<file;
	terminate(error.str());
	}

while ((! data_file.eof())) {
	data_file.getline(line, sizeof(line));
	if(strlen(line) && line[0] != '#' ) // skip comments
		{// effective line
		l++;
		parse(line,c);//parse line into c (vector of string)
		if(c[0]=="RESIDUE" && c.size()==3) {
			sscanf(c[1].c_str(),"%d",&resid);
			sscanf(c[2].c_str(),"%d",&flag);
			NR++;
			if (flag!=0) nfit++;
			if (R.find(resid) == R.end()) // new residue
				R.insert(resid); 
			else terminate("Redundant Residue Exists");
			npts=0;
			}
		else if(c[0]=="&") {
			if (npts > NX) NX = npts;
			}
		else npts++;
		} // effective line
	} // EOF

string a(INDENT,' ');
cout <<a<<left<<setw(40-INDENT);
cout <<"name of data file"<<": "<<file<<endl;
cout <<a<<left<<setw(40-INDENT);
cout <<"number of lines"<<": "<<l<<endl;
cout <<a<<left<<setw(40-INDENT);
cout <<"number of residues to be fitted"<<": "<<nfit;
cout <<" out of "<<NR<<endl;
cout <<a<<left<<setw(40-INDENT);
cout <<"number of maximum data points "<<": "<<NX<<endl;
}

/* Read DATAFILE  */
void read_data (const char *name,CSP &A)
{
char line[MAXSTR];
vector <string> c;//parsed columns
double x,y,dy;
int i,r,resid,flag;

r=0;//residue index
ifstream data_file(name);
while (!data_file.eof()) {
data_file.getline(line, sizeof(line));
if(strlen(line) && line[0] != '#') {// effective line

	parse(line,c);//parse line into c (vector of string)

	if (c[0]=="RESIDUE" && c.size()==3) {
		sscanf(c[1].c_str(),"%d",&resid);
		sscanf(c[2].c_str(),"%d",&flag);
		A.resid[r]=resid;
		if (flag!=0) A.flag[r] = true;
		i=0;
		}

	else if(c[0]=="&") 
		r++; // inc. residue index. r=[0...NR-1]

	else { // CSP data
		sscanf(c[0].c_str(),"%lf",&x);
		sscanf(c[1].c_str(),"%lf",&y);
		sscanf(c[2].c_str(),"%lf",&dy);
		gsl_matrix_set(A.x,r,i,x);
		gsl_matrix_set(A.y,r,i,y);
		gsl_matrix_set(A.dy,r,i,dy);
		i++;
		}
	}// effective line
}// EOF
data_file.close();
}

void read_parameter(const char *name, CSP &A)
{
string left(INDENT,' ');
char line[MAXSTR];
vector <string> c;//parsed columns
int i,l;
bool understood;

ifstream parameter_file(name);
if (parameter_file.fail()) {
	string error("cannot open parameter file");
	terminate(error);
	}

l=0;
while (! parameter_file.eof()) {

	parameter_file.getline(line, sizeof(line));
	l++;

	if(!strlen(line) || line[0] == '#') continue; /* skip comments */
	
	parse(line,c);//parse line into c (vector of string)

	understood = false;

	/* fit [-critx2 #] [-mc #] [-mctrim #] */
        if(c[0]=="fit") {
                understood = true;
                for (i=1;i<c.size();i++) {
			if (c[i]=="-mc") {
				sscanf(c[i+1].c_str(),"%d",&A.MC);
				cout << left;
				cout << "number of MC simulations "<< A.MC;
				cout << endl;
				i++;
				}
			else if (c[i]=="-mctrim") {
				sscanf(c[i+1].c_str(),"%d",&A.MC_trim);
				cout << left;
				cout << "number of MC simulations trimmed ";
				cout << "2*"<< A.MC_trim;
				cout << endl;
				i++;
				}
			else if (c[i]=="-critx2") {
				sscanf(c[i+1].c_str(),"%lf",&A.critx2);
				cout << left;
				cout << "acceptance for Chi-square criteria ";
				cout << A.critx2;
				cout << endl;
				i++;
				}
			else {
				stringstream error;
				error << "parameter file at line " <<l<<": ";
				error << "syntex error ";
				terminate(error.str());
				}
			}
		}// fit

	/* par */
	if(c[0]=="par") {
		for (i=0;i<A.NP;i++) {
			if (strcmp(c[1].c_str(),A.pid[i])==0 && c.size()==8) {
				understood = true;
				sscanf(c[2].c_str(),"%lf",&A.lb[i]);
				sscanf(c[3].c_str(),"%lf",&A.ub[i]);
				sscanf(c[4].c_str(),"%lf",&A.step[i]);
				sscanf(c[5].c_str(),"%lf",&A.conv[i]);
				sscanf(c[6].c_str(),"%lf",&A.Kl[i]);
				sscanf(c[7].c_str(),"%lf",&A.Ku[i]);
				/*
				string a(INDENT,' ');
				cout << left<< A.pid[i]<<endl;
				cout <<a<<left<<setw(40-INDENT);
				cout << "lower bound"<<": "<<A.lb[i]<<endl;
				cout <<a<<left<<setw(40-INDENT);
				cout << "upper bound"<<": "<<A.ub[i]<<endl;
				cout <<a<<left<<setw(40-INDENT);
				cout << "grids step"<<": "<<A.step[i]<<endl;
				cout <<a<<left<<setw(40-INDENT);
				cout << "convergency limit";
				cout << ": "<<A.conv[i]<<endl;
				cout <<a<<left<<setw(40-INDENT);
				cout << "penalty for lower bound";
				cout << ": "<<A.Kl[i]<<endl;
				cout <<a<<left<<setw(40-INDENT);
				cout << "penalty for upper bound";
				cout << ": "<<A.Ku[i]<<endl;
				*/
				}
			}
		}

	/* fix */
	if(c[0]=="fix") {
		for (i=0;i<A.NP;i++) 
			if (strcmp(c[1].c_str(),A.pid[i]) == 0) {
				understood = true;
				sscanf(c[2].c_str(),"%lf",&A.p[i]);
				cout << left;
				cout << "fix " << c[1];
				cout << " = " << A.p[i] << endl;
				SET_FIXED(A.pstat[i]);
				}
		}

	if (!understood) {
		stringstream error;
		error << "parameter file at line " <<l<<": ";
		error << "syntex error";
		terminate(error.str());
		}

	}// while

	parameter_file.close();
}

void parse (const char *line,vector <string> &c)
{
	const int l=strlen(line);
	int i,p;
	string s=line;
	c.clear();
    for (p=0,i=0;i<=l;i++) {
    	if (line[i]==' ' || line[i]=='\t' || line[i]=='\0') {
            c.push_back(s.substr(p,i-p));
		while (line[i]==' ' || line[i]=='\t') i++;
            p=i;
            }
	if (line[i]=='#') break;
	}
}
