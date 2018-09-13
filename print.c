#include "gfit.h"

using namespace std;

void print_help()
{
/* usage */
cout <<"\n";
cout <<" Global FITting ver. "<<VERSION<<" by Sung-Hun Bae\n\n";
cout <<" gfit [-d data file] [-p parameter file(default:gfit.par)]\n\n";
cout <<left<<setw(40)<<" -h : [H]elp (this message)";
cout <<" -v : [V]erbose data file\n";
cout <<left<<" gfit -d a.dat\n";
cout <<left<<" gfit -d a.dat -p x.par\n";
cout <<"\n";
exit(1);
}

void terminate(const string error)
{
cout << endl << endl;
cout << PROMPT << "ERROR " << error << endl;
cout << PROMPT << "terminated" << endl;
exit(1);
}

void print_selected(CSP &A)
{
int c,r;
string indent(INDENT,' ');
for (c=0,r=0;r<A.NR;r++) {
 	if(A.flag[r]) {
		if (c%10==0) cout << indent;
		printf("%4d ",A.resid[r]);
		if (c%10==9) cout << endl;
		c++;
		}//10 residues per line
	}
if (c%10!=0) cout << endl;
}

void print_grid_search (CSP &A)
{
string indent(INDENT,' ');
cout << indent;
for (int i=0;i<A.NP;i++)
	if (!IS_FIXED(A.pstat[i])) 
		printf("%s %.4f ",A.pid[i],A.p[i]);
cout << endl;
}

void print_iteration(int iter,CSP &A)
{
const extern bool opt_1,opt_3,opt_4,opt_5,opt_n;
string indent(INDENT,' ');
string hline(70,'=');
cout <<indent<<hline<<endl;
cout << PROMPT << "iteration ";
printf("%2d ",iter);
printf("Kd %8.05f",A.p[_Kd_]);
cout <<endl;// in order to empty output buffer
cout <<indent<<hline<<endl;
}

void print_boundary(CSP &A)
{
string indent(INDENT,' ');
// header
printf("%s",indent.c_str());
printf("   %-6s %8s [ %9s ... %9s | %9s ] %s\n",
"Param ","   Value","Lower BND","Upper BND","Grid Step","Conv.Limit");
for (int i=0;i<A.NP;i++) {
	printf("%s",indent.c_str());
	if (!IS_GLOBL(A.pstat[i])) printf("L"); else printf("G");
	if (!IS_FIXED(A.pstat[i])) printf("*"); else printf(" ");
	printf(" ");
	printf("%-5s %8.03e [ %9.03e ... %9.03e | %9.03e ] %.2E\n",
	A.pid[i],A.p[i],A.lb[i],A.ub[i],A.step[i],A.conv[i]);
	}
}

void print_verb(const string s) {cout << PROMPT << s << endl;}
